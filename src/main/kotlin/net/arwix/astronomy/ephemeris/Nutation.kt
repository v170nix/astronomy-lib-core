/*
 * Copyright 2017 Vitaliy Sheyanov vit.onix@gmail.com
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package net.arwix.astronomy.ephemeris

import net.arwix.astronomy.annotation.Ecliptic
import net.arwix.astronomy.annotation.Geocentric
import net.arwix.astronomy.core.ARCSEC_TO_RAD
import net.arwix.astronomy.core.DEG_TO_RAD
import net.arwix.astronomy.core.vector.Matrix
import net.arwix.astronomy.core.vector.Vector
import java.lang.Math.*

sealed class Nutation(val t: Double) : NutationTransformation {

    class IAU1980(t: Double, obliquity: Obliquity) : Nutation(t), NutationTransformation by ImplNutationTransform(calcNutation_IAU1980(t), obliquity.meanEps)
    class IAU2000(t: Double, obliquity: Obliquity) : Nutation(t), NutationTransformation by ImplNutationTransform(calcNutation_IAU2000(t), obliquity.meanEps)
    // Apply precession adjustments, see Wallace & Capitaine, 2006, Eqs.5
    class IAU2006(t: Double, obliquity: Obliquity) : Nutation(t), NutationTransformation by ImplNutationTransform(calcNutation_IAU2000(t).let {
        NutationResult(it.longitude * (1.0 + (0.4697E-6 - 2.7774E-6 * t)), it.obliquity * (1.0 + (2.7774E-6 * t)))
    }, obliquity.meanEps)

    class Fast(t: Double, obliquity: Obliquity) : Nutation(t), NutationTransformation by ImplNutationTransform(getFastNutation(t), obliquity.meanEps)
}

internal data class NutationResult(val longitude: Double, val obliquity: Double)

interface NutationTransformation {
    @Ecliptic val eclipticMatrix: Matrix

    @Ecliptic fun applyNutationToEclipticVector(@Ecliptic vector: Vector): Vector
    @Ecliptic fun removeNutationFromEclipticVector(@Ecliptic vector: Vector): Vector

    @Geocentric val geocentricMatrix: Matrix

    /**
     * Nutates equatorial coordinates from mean dynamical equator and equinox of date to true
     * equator and equinox, or the opposite. See AA Explanatory Supplement, page 114-115.
     * For dates between 1900 and 2100 and the flag to prefer precision in the ephemeris
     * object set to false, the approximate code by Peter Duffet for nutation is used.
     *
     * @param vector equatorial coordinates
     * @return Output equatorial coordinates
     */
    @Geocentric fun applyNutationToGeocentricVector(@Geocentric vector: Vector): Vector
    @Geocentric fun removeNutationFromGeocentricVector(@Geocentric vector: Vector): Vector
}

private class ImplNutationTransform(nutationAngles: NutationResult, meanEps: Double) : NutationTransformation {

    override val eclipticMatrix = getEclipticMatrix(nutationAngles)
    override val geocentricMatrix = getGeocentricMatrix(nutationAngles, meanEps)


    companion object {
        private fun getGeocentricMatrix(nutationAngles: NutationResult, meanEps: Double): Matrix {
            val nut = nutationAngles
            val oblm = meanEps
            val dpsi = nut.longitude
            return Matrix(Matrix.Axis.X, -(oblm + nut.obliquity)) * Matrix(Matrix.Axis.Z, -dpsi) * Matrix(Matrix.Axis.X, oblm)
        }

        private fun getEclipticMatrix(nutationAngles: NutationResult): Matrix {
            val nut = nutationAngles
            val dpsi = nut.longitude
            return Matrix(Matrix.Axis.X, -(nut.obliquity)) * Matrix(Matrix.Axis.Z, -dpsi)
        }
    }

    override fun applyNutationToEclipticVector(vector: Vector) = eclipticMatrix * vector
    override fun removeNutationFromEclipticVector(vector: Vector) = vector * eclipticMatrix
    override fun applyNutationToGeocentricVector(vector: Vector) = geocentricMatrix * vector
    override fun removeNutationFromGeocentricVector(vector: Vector) = vector * geocentricMatrix

}

private fun getFastNutation(t: Double): NutationResult {

    val T2: Double = t * t

    var aq = 100.0021358 * t
    var b = 360.0 * (aq - floor(aq))
    val L1 = 279.6967 + .000303 * T2 + b
    val L2 = 2.0 * L1 * DEG_TO_RAD

    aq = 1336.855231 * t
    b = 360.0 * (aq - floor(aq))
    val D1 = 270.4342 - .001133 * T2 + b
    val D2 = 2.0 * D1 * DEG_TO_RAD

    aq = 99.99736056 * t
    b = 360.0 * (aq - floor(aq))
    val M1 = (358.4758 - .00015 * T2 + b) * DEG_TO_RAD

    aq = 1325.552359 * t
    b = 360.0 * (aq - floor(aq))
    val M2 = (296.1046 + .009192 * T2 + b) * DEG_TO_RAD

    aq = 5.372616667 * t
    b = 360.0 * (aq - floor(aq))
    val N1 = (259.1833 + .002078 * T2 - b) * DEG_TO_RAD
    val N2 = 2.0 * N1 // * Constant.DEG_TO_RAD // A bug in Peter Duffett's code !?

    var DP = (-17.2327 - .01737 * t) * sin(N1)
    DP += (-1.2729 - .00013 * t) * sin(L2) + .2088 * sin(N2)
    DP += -.2037 * sin(D2) + (.1261 - .00031 * t) * sin(M1)
    DP += .0675 * sin(M2) - (.0497 - .00012 * t) * sin(L2 + M1)
    DP += -.0342 * sin(D2 - N1) - .0261 * sin(D2 + M2)
    DP += +.0214 * sin(L2 - M1) - .0149 * sin(L2 - D2 + M2)
    DP += .0124 * sin(L2 - N1) + .0114 * sin(D2 - M2)

    var DOSR = (9.21 + .00091 * t) * cos(N1)
    DOSR += (.5522 - .00029 * t) * cos(L2) - .0904 * cos(N2)
    DOSR += .0884 * cos(D2) + .0216 * cos(L2 + M1)
    DOSR += .0183 * cos(D2 - N1) + .0113 * cos(D2 + M2)
    DOSR += -.0093 * cos(L2 - M1) - .0066 * cos(L2 - N1)

    DP *= ARCSEC_TO_RAD
    DOSR *= ARCSEC_TO_RAD

    return NutationResult(DP, DOSR)
}


