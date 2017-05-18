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

package net.arwix.astronomy.ephemeris.precession

import net.arwix.astronomy.annotation.Geocentric
import net.arwix.astronomy.core.ARCSEC_TO_RAD
import net.arwix.astronomy.core.PI2
import net.arwix.astronomy.core.vector.Matrix
import net.arwix.astronomy.core.vector.RectangularVector
import net.arwix.astronomy.core.vector.Vector


object NewPrecession {

    fun getVondrakMatrix(T: Double): Matrix {

        var w = PI2 * T
        val EPS0 = 84381.406 * ARCSEC_TO_RAD
        val (PSIA, OMEGAA, CHIA) = (0..13).fold(RectangularVector(0.0, 0.0, 0.0), { acc, i ->
            var a = w / xyper[i][0]
            var s = Math.sin(a)
            var c = Math.cos(a)
            acc.x += c * xyper[i][1] + s * xyper[i][3]
            acc.y += c * xyper[i][2] + s * xyper[i][4]

            a = w / zper[i][0]
            s = Math.sin(a)
            c = Math.cos(a)
            acc.z += c * zper[i][1] + s * zper[i][2]
            return@fold acc
        }).let {
            w = 1.0
            for (j in 0..3) {
                it.x += xypol[0][j] * w
                it.y += xypol[1][j] * w
                it.z += xypol[2][j] * w
                w *= T
            }
            return@let it
        } * ARCSEC_TO_RAD

        // COMPUTE ELEMENTS OF PRECESSION ROTATION MATRIX
        // EQUIVALENT TO R3(CHI_A)R1(-OMEGA_A)R3(-PSI_A)R1(EPSILON_0)

        return Matrix(Matrix.Axis.Z, CHIA) * Matrix(Matrix.Axis.X, -OMEGAA) * Matrix(Matrix.Axis.Z, -PSIA) * Matrix(Matrix.Axis.X, EPS0)

    }

    fun getIAU2000Matrix(T: Double, isFromJ2000ToApparent: Boolean = true): Matrix {
        val T0: Double = if (isFromJ2000ToApparent) T else 0.0

        var EPS0 = 84381.448
        var PSIA = ((((-0.0 * T + 0.0) * T - 0.001147) * T - 1.07259) * T + 5038.7784) * T - 0.29965 * T0
        var OMEGAA = ((((+0.0 * T - 0.0) * T - 0.007726) * T + 0.05127) * T - 0.0) * T + EPS0 - 0.02524 * T0
        var CHIA = ((((-0.0 * T + 0.0) * T - 0.001125) * T - 2.38064) * T + 10.5526) * T

        EPS0 *= ARCSEC_TO_RAD
        PSIA *= ARCSEC_TO_RAD
        OMEGAA *= ARCSEC_TO_RAD
        CHIA *= ARCSEC_TO_RAD
        return Matrix(Matrix.Axis.Z, CHIA) * Matrix(Matrix.Axis.X, -OMEGAA) * Matrix(Matrix.Axis.Z, -PSIA) * Matrix(Matrix.Axis.X, EPS0)
    }

    fun getIAU2006Matrix(T: Double): Matrix {
        var EPS0 = 84381.406
        var PSIA = ((((-0.0000000951 * T + 0.000132851) * T - 0.00114045) * T - 1.0790069) * T + 5038.481507) * T
        var OMEGAA = ((((+0.0000003337 * T - 0.000000467) * T - 0.00772503) * T + 0.0512623) * T - 0.025754) * T + EPS0
        var CHIA = ((((-0.0000000560 * T + 0.000170663) * T - 0.00121197) * T - 2.3814292) * T + 10.556403) * T

        EPS0 *= ARCSEC_TO_RAD
        PSIA *= ARCSEC_TO_RAD
        OMEGAA *= ARCSEC_TO_RAD
        CHIA *= ARCSEC_TO_RAD
        return Matrix(Matrix.Axis.Z, CHIA) * Matrix(Matrix.Axis.X, -OMEGAA) * Matrix(Matrix.Axis.Z, -PSIA) * Matrix(Matrix.Axis.X, EPS0)
    }

    /**
     * Precession following IAU2000 definitions. From SOFA software library.
     * Reference: Capitaine et al., Astronomy & Astrophysics 400, 1145-1154,
     * 2003. See also Lieske et al. 1977.
     *
     * @param JD0 Julian day of input vector (equatorial rectangular).
     * @param JD Julian day of output. Either JD or JD0 must be equal to Constant.J2000.
     * @param R Input vector.
     * @return Vector refered to mean equinox and equator of JD.
     * @throws JPARSECException If both JD and JD0 are not equal to J2000 epoch.
    </P> */
    fun precessionIAU2000(T: Double, R: Vector, isFromJ2000ToApparent: Boolean = true): Vector {

        return if (isFromJ2000ToApparent) getIAU2000Matrix(T, isFromJ2000ToApparent) * R else R * getIAU2000Matrix(T, isFromJ2000ToApparent)

    }

    /**
     * Precession following Capitaine et al. 2003.
     * <P>
     * Capitaine formula of precession is to be officially adopted by the IAU,
     * see recommendation in the report of the IAU Division I Working Group on
     * Precession and the Ecliptic (Hilton et al. 2006, Celest. Mech., 94,
     * 351-367).
     * Reference: Capitaine et al., Astronomy & Astrophysics 412, 567-586,
     * 2003.
     * @param JD0 Julian day of input vector (equatorial rectangular).
     * @param JD Julian day of output. Either JD or JD0 must be equal to Constant.J2000.
     * @param R Input vector.
     *
     * @return Vector referred to mean equinox and equator of JD.
     * */
    fun precessionIAU2006(T: Double, R: Vector, isFromJ2000ToApparent: Boolean = true): Vector {
        return if (isFromJ2000ToApparent) getIAU2006Matrix(T) * R else R * getIAU2006Matrix(T)
    }


    /**
     * Precession following Vondrak et al. 2011. See A&amp;A 534, A22.

     * @param JD0 Julian day of input vector (equatorial rectangular).
     * @param JD Julian day of output. Either JD or JD0 must be equal to
     *        Constant.J2000.
     * @param R Input vector.
     * @return Vector referred to mean equinox and equator of JD.
     */
    fun precessionVondrak2011(T: Double, @Geocentric R: Vector, isFromJ2000ToApparent: Boolean = true): Vector {

        return if (isFromJ2000ToApparent) {
            getVondrakMatrix(T) * R
        } else {
            getVondrakMatrix(T).transpose() * R
        }
    }


    /**
     * Precess rectangular equatorial coordinates from J2000 epoch.
     *
     * @param JD Equinox of the output in Julian day (TT).
     * @param R Array with x, y, z.
     * @param eph Ephemeris properties.
     * @return Array with corrected x, y, z.
     */
    fun precessFromJ2000(T: Double, @Geocentric R: Vector): Vector {
        return precessionIAU2006(T, R, true)
        // return getVondrakMatrix(T) * R
    }

    /**
     * Precess rectangular equatorial coordinates to J2000 epoch.
     *
     * @param JD Equinox of the input in Julian day (TT).
     * @param R Array with x, y, z.
     * @param eph The ephemeris properties.
     * @return Array with corrected x, y, z.
     */
    fun precessToJ2000(T: Double, @Geocentric R: Vector): Vector {
        //    return R * getVondrakMatrix(T)
        return precessionIAU2000(T, R, false)
    }

    /**
     * Precess rectangular equatorial coordinates between two epochs.
     *
     * @param fromDatejT Equinox of the input in Julian day (TT).
     * @param toDatejT Equinox of the output in Julian day (TT).
     * @param R Array with x, y, z.
     * @return Array with corrected x, y, z.
     */
    fun precess(fromDatejT: Double, toDatejT: Double, R: Vector): Vector {


//        return precessionIAU2000(toDatejT) * (R * (fromDatejT))

//        return getVondrakMatrix(toDatejT) * (R * getVondrakMatrix(fromDatejT))
        val v = precessToJ2000(fromDatejT, R)
        return precessFromJ2000(toDatejT, v)
    }


}

private val xypol = arrayOf(
        doubleArrayOf(8473.343527, 5042.7980307, -0.00740913, 289E-9),
        doubleArrayOf(84283.175915, -0.4436568, 0.00000146, 151E-9),
        doubleArrayOf(-19.657270, 0.0790159, 0.00001472, -61E-9))

private val xyper = arrayOf(
        doubleArrayOf(402.90, -22206.325946, 1267.727824, -3243.236469, -8571.476251),
        doubleArrayOf(256.75, 12236.649447, 1702.324248, -3969.723769, 5309.796459),
        doubleArrayOf(292.00, -1589.008343, -2970.553839, 7099.207893, -610.393953),
        doubleArrayOf(537.22, 2482.103195, 693.790312, -1903.696711, 923.201931),
        doubleArrayOf(241.45, 150.322920, -14.724451, 146.435014, 3.759055),
        doubleArrayOf(375.22, -13.632066, -516.649401, 1300.630106, -40.691114),
        doubleArrayOf(157.87, 389.437420, -356.794454, 1727.498039, 80.437484),
        doubleArrayOf(274.20, 2031.433792, -129.552058, 299.854055, 807.300668),
        doubleArrayOf(203.00, 363.748303, 256.129314, -1217.125982, 83.712326),
        doubleArrayOf(440.00, -896.747562, 190.266114, -471.367487, -368.654854),
        doubleArrayOf(170.72, -926.995700, 95.103991, -441.682145, -191.881064),
        doubleArrayOf(713.37, 37.070667, -332.907067, -86.169171, -4.263770),
        doubleArrayOf(313.00, -597.682468, 131.337633, -308.320429, -270.353691),
        doubleArrayOf(128.38, 66.282812, 82.731919, -422.815629, 11.602861))

private val zper = arrayOf(
        doubleArrayOf(402.90, -13765.924050, -2206.967126),
        doubleArrayOf(256.75, 13511.858383, -4186.752711),
        doubleArrayOf(292.00, -1455.229106, 6737.949677),
        doubleArrayOf(537.22, 1054.394467, -856.922846),
        doubleArrayOf(375.22, -112.300144, 957.149088),
        doubleArrayOf(157.87, 202.769908, 1709.440735),
        doubleArrayOf(274.20, 1936.050095, 154.425505),
        doubleArrayOf(202.00, 327.517465, -1049.071786),
        doubleArrayOf(440.00, -655.484214, -243.520976),
        doubleArrayOf(170.72, -891.898637, -406.539008),
        doubleArrayOf(315.00, -494.780332, -301.504189),
        doubleArrayOf(136.32, 585.492621, 41.348740),
        doubleArrayOf(128.38, -333.322021, -446.656435),
        doubleArrayOf(490.00, 110.512834, 142.525186))