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
import net.arwix.astronomy.annotation.Equatorial
import net.arwix.astronomy.core.ARCSEC_TO_RAD
import net.arwix.astronomy.core.MINUTES_PER_DEGREE
import net.arwix.astronomy.core.PI2
import net.arwix.astronomy.core.SECONDS_PER_DEGREE
import net.arwix.astronomy.core.vector.Matrix
import net.arwix.astronomy.core.vector.Vector
import net.arwix.astronomy.math.polynomialSum
import java.lang.Math.cos
import java.lang.Math.sin

sealed class Obliquity(val t: Double) : ObliquityTransformation {
    class WILLIAMS_1994(t: Double) : Obliquity(t), ObliquityTransformation by ImplObliquityTransform(getEps(t, rvalStart_WIL, coeffs_WIL))
    class SIMON_1994(t: Double) : Obliquity(t), ObliquityTransformation by ImplObliquityTransform(getEps(t, rvalStart_SIM, coeffs_SIM))
    class LASKAR_1996(t: Double) : Obliquity(t), ObliquityTransformation by ImplObliquityTransform(getEps(t, rvalStart_LAS, coeffs_LAS))
    class IAU_1976(t: Double) : Obliquity(t), ObliquityTransformation by ImplObliquityTransform(getEps(t, rvalStart_IAU, coeffs_IAU))
    class IAU_2006(t: Double) : Obliquity(t), ObliquityTransformation by ImplObliquityTransform(getEps(t, rvalStart_CAP, coeffs_CAP))
    class VONDRAK_2011(t: Double) : Obliquity(t), ObliquityTransformation by ImplObliquityTransform(getVondrakEps(t))
}

interface ObliquityTransformation {

    val meanEps: Double
    val eclipticToEquatorialMatrix: Matrix
    val equatorialToEclipticMatrix: Matrix

    @Equatorial fun rotateFromEclipticToEquatorial(@Ecliptic vector: Vector): Vector
    @Ecliptic fun rotateFromEquatorialToEcliptic(@Equatorial vector: Vector): Vector
}

private class ImplObliquityTransform(override val meanEps: Double) : ObliquityTransformation {

    override val eclipticToEquatorialMatrix = Matrix(Matrix.Axis.X, -meanEps)
    override val equatorialToEclipticMatrix = eclipticToEquatorialMatrix.transpose()

    override fun rotateFromEclipticToEquatorial(vector: Vector) = eclipticToEquatorialMatrix * vector
    override fun rotateFromEquatorialToEcliptic(vector: Vector) = equatorialToEclipticMatrix * vector

}

private fun getEps(t: Double, rvalStart: Double, coeffs: DoubleArray): Double {
    return (rvalStart + coeffs.polynomialSum(t / 100.0)) * ARCSEC_TO_RAD
}

private fun getVondrakEps(t: Double): Double {
    return (PI2 * t)
            .let { w ->
                xyper.sumByDouble { doubles -> (w / doubles[0]).let { a -> cos(a) * doubles[1] + sin(a) * doubles[2] } } +
                        xypol.polynomialSum(t)
            } * ARCSEC_TO_RAD
}

// Capitaine et al. 2003, Hilton et al. 2006
val rvalStart_CAP = 23.0 * SECONDS_PER_DEGREE + 26.0 * MINUTES_PER_DEGREE + 21.406
val coeffs_CAP = doubleArrayOf(0.0, -4683.6769, -1.831, 2003.400, -57.6, -434.0)

// Simon et al., 1994
val OBLIQ_COEFFS = 10
val rvalStart_SIM = 23.0 * SECONDS_PER_DEGREE + 26.0 * MINUTES_PER_DEGREE + 21.412
val coeffs_SIM = doubleArrayOf(0.0, -4680.927, -1.52, 1998.9, -51.38, -249.67, -39.05, 7.12, 27.87, 5.79, 2.45)

// Williams et al., DE403 Ephemeris
val rvalStart_WIL = 23.0 * SECONDS_PER_DEGREE + 26.0 * MINUTES_PER_DEGREE + 21.406173
val coeffs_WIL = doubleArrayOf(0.0, -4683.396, -1.75, 1998.9, -51.38, -249.67, -39.05, 7.12, 27.87, 5.79, 2.45)

// Laskar et al.
/*
 * This expansion is from Laskar, cited above. Bretagnon and Simon say,
 * in Planetary Programs and Tables, that it is accurate to 0.1" over a
 * span of 6000 years. Laskar estimates the precision to be 0.01" after
 * 1000 years and a few seconds of arc after 10000 years.
 */
val rvalStart_LAS = 23.0 * SECONDS_PER_DEGREE + 26.0 * MINUTES_PER_DEGREE + 21.448
val coeffs_LAS = doubleArrayOf(0.0, -4680.93, -1.55, 1999.25, -51.38, -249.67, -39.05, 7.12, 27.87, 5.79, 2.45)

// IAU 1976
val rvalStart_IAU = 23.0 * SECONDS_PER_DEGREE + 26.0 * MINUTES_PER_DEGREE + 21.448
val coeffs_IAU = doubleArrayOf(0.0, -4681.5, -5.9, 1813.0)


private val xypol by lazy { doubleArrayOf(84028.206305, 0.3624445, -0.00004039, -110E-9) }

private val xyper by lazy {
    arrayOf(
            doubleArrayOf(409.90, 753.872780, -1704.720302),
            doubleArrayOf(396.15, -247.805823, -862.308358),
            doubleArrayOf(537.22, 379.471484, 447.832178),
            doubleArrayOf(402.90, -53.880558, -889.571909),
            doubleArrayOf(417.15, -90.109153, 190.402846),
            doubleArrayOf(288.92, -353.600190, -56.564991),
            doubleArrayOf(4043.00, -63.115353, -296.222622),
            doubleArrayOf(306.00, -28.248187, -75.859952),
            doubleArrayOf(277.00, 17.703387, 67.473503),
            doubleArrayOf(203.00, 38.911307, 3.014055))
}