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

import net.arwix.astronomy.annotation.Equatorial
import net.arwix.astronomy.annotation.Geocentric
import net.arwix.astronomy.core.ARCSEC_TO_RAD
import net.arwix.astronomy.core.PI2
import net.arwix.astronomy.core.vector.Matrix
import net.arwix.astronomy.core.vector.RectangularVector
import net.arwix.astronomy.core.vector.Vector


enum class PrecessionMethod {

    /**
     * Constant ID for selecting IAU 1976 formulae of precession, nutation (IAU
     * 1980), and Greenwich mean sidereal time. This will use
     * old formulae that will match results from IMCCE ephemeris server. You
     * may consider using this formulae for VSOP theory. See
     * J. H. Lieske, T. Lederle, W. Fricke, and B. Morando, "Expressions for the
     * Precession Quantities Based upon the IAU (1976) System of Astronomical
     * Constants," Astronomy and Astrophysics 58, 1-16 (1977).
     */
    IAU_1976,

    /**
     * Constant ID for selecting IAU 2000 formulae of precession, obliquity,
     * nutation, and Greenwich mean sidereal time. Pole movement is considered
     * when possible. Planetary rotation models are set to IAU 2000 resolutions.
     * @deprecated It is not recommended since the precession is inconsistent.
     */
    IAU_2000,

    /**
     * Constant ID for selecting IAU 2006 formulae of obliquity, nutation, and
     * Greenwich mean sidereal time. Note that the precession algorithm is from
     * Capitaine et al. 2003 (Astronomy &amp; Astrophysics 412, 567-586, 2003), officially
     * adopted by the IAU as a replacement of the precession part of IAU2000A
     * precession-nutation model. See also Hilton et al. 2006. Pole movement is
     * considered when possible. Should be used by default for better precission
     * in modern theories like JPL DE405 and above. Planetary rotation models are
     * set to IAU 2006 resolutions.
     */
    IAU_2006,

    /**
     * Same as IAU2006, but planetary rotation models are those recommended by
     * the IAU working group on carthographic coordinates, in 2009.
     */
    IAU_2009,

    /**
     * Precession following Vondrak et al. 2011. See A&amp;A 534, A22.
     */
    VONDRAK_2011,

    /**
     * Constant ID for selecting JPL DE403/404/405/406 formulae for precession,
     * obliquity, nutation (IAU 1980), and Greenwich mean sidereal time. Quite
     * similar to Williams formulae. Adequate for planets using Moshier
     * method, Series96, or JPL DE40x ephemerides.
     */
    JPL_DE4xx,

    /**
     * Constant ID for selecting Laskar formulae of precession, nutation (IAU
     * 1980), and Greenwich mean sidereal time. See J. Laskar,
     * "Secular terms of classical planetary theories using the results of
     * general theory," Astronomy and Astrophysics 157, 59070 (1986).
     */
    LASKAR_1986,

    /**
     * Constant ID for selecting Williams formulae of precession (DE403 JPL
     * Ephemeris), nutation (IAU 1980), obliquity, and Greenwich mean sidereal
     * time. See James G. Williams, "Contributions to the Earth's obliquity rate,
     * precession, and nutation," Astron. J. 108, 711-724 (1994). It is convenient
     * to use this when obtaining ephemeris of the Moon using Moshier method.
     */
    WILLIAMS_1994,

    /**
     * Constant ID for selecting SIMON formulae of precession, obliquity,
     * nutation (IAU 1980), and Greenwich mean sidereal time. See
     * J. L. Simon, P. Bretagnon, J. Chapront, M. Chapront-Touze', G. Francou,
     * and J. Laskar, "Numerical Expressions for precession formulae and mean
     * elements for the Moon and the planets," Astronomy and Astrophysics 282,
     * 663-683 (1994).
     */
    SIMON_1994;

    companion object {

        private val WILLIAMS_1994_Matrices by lazy {
            arrayOf(
                    doubleArrayOf(-8.66e-10, -4.759e-8, 2.424e-7, 1.3095e-5, 1.7451e-4, -1.8055e-3, -0.235316, 0.076, 110.5407, 50287.70000),
                    /* Pi from Williams' 1994 paper, in radians. */
                    doubleArrayOf(6.6402e-16, -2.69151e-15, -1.547021e-12, 7.521313e-12, 1.9e-10, -3.54e-9, -1.8103e-7, 1.26e-7, 7.436169e-5, -0.04207794833, 3.052115282424),
                    doubleArrayOf(1.2147e-16, 7.3759e-17, -8.26287e-14, 2.503410e-13, 2.4650839e-11, -5.4000441e-11, 1.32115526e-9, -6.012e-7, -1.62442e-5, 0.00227850649, 0.0))
        }
        private val JPL_DE4xx_Matrices by lazy {
            arrayOf(
                    doubleArrayOf(-8.66e-10, -4.759e-8, 2.424e-7, 1.3095e-5, 1.7451e-4, -1.8055e-3, -0.235316, 0.076, 110.5414, 50287.91959),
                    /* Pi from Williams' 1994 paper, in radians. No change in DE403. */
                    doubleArrayOf(6.6402e-16, -2.69151e-15, -1.547021e-12, 7.521313e-12, 1.9e-10, -3.54e-9, -1.8103e-7, 1.26e-7, 7.436169e-5, -0.04207794833, 3.052115282424),
                    doubleArrayOf(1.2147e-16, 7.3759e-17, -8.26287e-14, 2.503410e-13, 2.4650839e-11, -5.4000441e-11, 1.32115526e-9, -6.012e-7, -1.62442e-5, 0.00227850649, 0.0))
        }
        private val SIMON_1994_Matrices by lazy {
            arrayOf(
                    doubleArrayOf(-8.66e-10, -4.759e-8, 2.424e-7, 1.3095e-5, 1.7451e-4, -1.8055e-3, -0.235316, 0.07732, 111.2022, 50288.200),
                    doubleArrayOf(6.6402e-16, -2.69151e-15, -1.547021e-12, 7.521313e-12, 1.9e-10, -3.54e-9, -1.8103e-7, 2.579e-8, 7.4379679e-5, -0.0420782900, 3.0521126906),
                    doubleArrayOf(1.2147e-16, 7.3759e-17, -8.26287e-14, 2.503410e-13, 2.4650839e-11, -5.4000441e-11, 1.32115526e-9, -5.99908e-7, -1.624383e-5, 0.002278492868, 0.0))
        }

        private val LASKAR_1986_Matrices by lazy {
            arrayOf(
                    doubleArrayOf(-8.66e-10, -4.759e-8, 2.424e-7, 1.3095e-5, 1.7451e-4, -1.8055e-3, -0.235316, 0.07732, 111.1971, 50290.966),
                    doubleArrayOf(6.6402e-16, -2.69151e-15, -1.547021e-12, 7.521313e-12, 6.3190131e-10, -3.48388152e-9, -1.813065896e-7, 2.75036225e-8, 7.4394531426e-5, -0.042078604317, 3.052112654975),
                    doubleArrayOf(1.2147e-16, 7.3759e-17, -8.26287e-14, 2.503410e-13, 2.4650839e-11, -5.4000441e-11, 1.32115526e-9, -5.998737027e-7, -1.6242797091e-5, 0.002278495537, 0.0))
        }

        private val IAU_1976_Matrices by lazy {
            arrayOf(
                    doubleArrayOf(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.006, 111.113, 50290.966),
                    doubleArrayOf(6.6402e-16, -2.69151e-15, -1.547021e-12, 7.521313e-12, 6.3190131e-10, -3.48388152e-9, -1.813065896e-7, 2.75036225e-8, 7.4394531426e-5, -0.042078604317, 3.052112654975),
                    doubleArrayOf(1.2147e-16, 7.3759e-17, -8.26287e-14, 2.503410e-13, 2.4650839e-11, -5.4000441e-11, 1.32115526e-9, -5.998737027e-7, -1.6242797091e-5, 0.002278495537, 0.0))
        }

        private val xypol by lazy {
            arrayOf(
                    doubleArrayOf(8473.343527, 5042.7980307, -0.00740913, 289E-9),
                    doubleArrayOf(84283.175915, -0.4436568, 0.00000146, 151E-9),
                    doubleArrayOf(-19.657270, 0.0790159, 0.00001472, -61E-9))
        }

        private val xyper by lazy {
            arrayOf(
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
        }

        private val zper by lazy {
            arrayOf(
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
        }

        /**
         * Precession following IAU2000 definitions. From SOFA software library.
         * Reference: Capitaine et al., Astronomy & Astrophysics 400, 1145-1154,
         * 2003. See also Lieske et al. 1977.
         */
        @Equatorial @Geocentric
        private fun getIAU2000Matrix(T: Double, isFromJ2000ToApparent: Boolean = true): Matrix {
            val T0: Double = if (isFromJ2000ToApparent) T else 0.0
            var EPS0 = 84381.448
            val (PSIA, OMEGAA, CHIA) = RectangularVector(
                    ((((-0.0 * T + 0.0) * T - 0.001147) * T - 1.07259) * T + 5038.7784) * T - 0.29965 * T0,
                    ((((+0.0 * T - 0.0) * T - 0.007726) * T + 0.05127) * T - 0.0) * T + EPS0 - 0.02524 * T0,
                    ((((-0.0 * T + 0.0) * T - 0.001125) * T - 2.38064) * T + 10.5526) * T) * ARCSEC_TO_RAD

            EPS0 *= ARCSEC_TO_RAD

            return Matrix(Matrix.Axis.Z, CHIA) * Matrix(Matrix.Axis.X, -OMEGAA) * Matrix(Matrix.Axis.Z, -PSIA) * Matrix(Matrix.Axis.X, EPS0)
        }

        /**
         * Precession following Capitaine et al. 2003.
         * Capitaine formula of precession is to be officially adopted by the IAU,
         * see recommendation in the report of the IAU Division I Working Group on
         * Precession and the Ecliptic (Hilton et al. 2006, Celest. Mech., 94, 351-367).
         *
         * Reference: Capitaine et al., Astronomy & Astrophysics 412, 567-586, 2003.
         */
        @Equatorial
        private fun getIAU2006Matrix(T: Double): Matrix {
            var EPS0 = 84381.406
            val (PSIA, OMEGAA, CHIA) = RectangularVector(
                    ((((-0.0000000951 * T + 0.000132851) * T - 0.00114045) * T - 1.0790069) * T + 5038.481507) * T,
                    ((((+0.0000003337 * T - 0.000000467) * T - 0.00772503) * T + 0.0512623) * T - 0.025754) * T + EPS0,
                    ((((-0.0000000560 * T + 0.000170663) * T - 0.00121197) * T - 2.3814292) * T + 10.556403) * T) * ARCSEC_TO_RAD

            EPS0 *= ARCSEC_TO_RAD

            return Matrix(Matrix.Axis.Z, CHIA) * Matrix(Matrix.Axis.X, -OMEGAA) * Matrix(Matrix.Axis.Z, -PSIA) * Matrix(Matrix.Axis.X, EPS0)
        }

        /**
         * Precession following Vondrak et al. 2011. See A&amp;A 534, A22.
         */
        @Equatorial @Geocentric
        private fun getVondrakMatrix(T: Double): Matrix {
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
    }

    @Equatorial
    val d = 1


    @Equatorial
    fun getEquatorialPrecessionMatrix(T: Double, isFromJ2000ToApparent: Boolean = true): Matrix {
        return when (this) {
            PrecessionMethod.IAU_1976 -> TODO()
            PrecessionMethod.IAU_2000 -> getIAU2000Matrix(T, isFromJ2000ToApparent)
            PrecessionMethod.IAU_2006 -> getIAU2006Matrix(T)
            PrecessionMethod.IAU_2009 -> getIAU2006Matrix(T)
            PrecessionMethod.VONDRAK_2011 -> getVondrakMatrix(T)
            PrecessionMethod.JPL_DE4xx -> TODO()
            PrecessionMethod.LASKAR_1986 -> TODO()
            PrecessionMethod.WILLIAMS_1994 -> TODO()
            PrecessionMethod.SIMON_1994 -> TODO()
        }.let {
            if (!isFromJ2000ToApparent) it.transpose() else it
        }
    }

    @Equatorial
    fun getEquatorialPrecessionMatrix(fromT: Double, toT: Double) = getEquatorialPrecessionMatrix(toT, true) * getEquatorialPrecessionMatrix(fromT, false)


    fun getEclipticPrecessionVector(T: Double, isFromJ2000ToApparent: Boolean = true): Matrix {

        val v = when (this) {
            PrecessionMethod.IAU_1976 -> Companion.IAU_1976_Matrices
            PrecessionMethod.IAU_2000 -> TODO()
            PrecessionMethod.IAU_2006 -> TODO()
            PrecessionMethod.IAU_2009 -> TODO()
            PrecessionMethod.VONDRAK_2011 -> TODO()
            PrecessionMethod.JPL_DE4xx -> Companion.JPL_DE4xx_Matrices
            PrecessionMethod.LASKAR_1986 -> Companion.LASKAR_1986_Matrices
            PrecessionMethod.WILLIAMS_1994 -> Companion.WILLIAMS_1994_Matrices
            PrecessionMethod.SIMON_1994 -> Companion.SIMON_1994_Matrices
        }
        var pA: Double
        var W: Double
        var z: Double
        val element1: DoubleArray
        val element2: DoubleArray
        val element3: DoubleArray
        var p1 = -1
        var p2 = -1
        var p3 = -1
        var i: Int

        /*
		 * Precession in longitude
		 */
        val T10 = T / 10.0 /* thousands of years */
        p1++
        element1 = v[0]
        pA = element1[p1]
        i = 0
        while (i < 9) {
            p1++
            pA = pA * T10 + element1[p1]
            i++
        }
        pA *= ARCSEC_TO_RAD * T10

        /*
		 * Node of the moving ecliptic on the JD0 ecliptic.
		 */
        p2++
        element2 = v[1]
        W = element2[p2]
        i = 0
        while (i < 10) {
            p2++
            W = W * T10 + element2[p2]
            i++
        }

        /*
		 * Rotate about new x axis by the inclination of the moving ecliptic on
		 * the JD0 ecliptic.
		 */
        p3++
        element3 = v[2]
        z = element3[p3]
        i = 0
        while (i < 10) {
            p3++
            z = z * T10 + element3[p3]
            i++
        }

        if (!isFromJ2000ToApparent) z = -z

        //   return Matrix(Matrix.Axis.Z, W) * Matrix(Matrix.Axis.X, z) * Matrix(Matrix.Axis.Z, -(W + pA))
        return Matrix(Matrix.Axis.Z, -(W + pA)) * Matrix(Matrix.Axis.X, z) * Matrix(Matrix.Axis.Z, W)

        // z = W // z = z

        //            Matrix(Matrix.Axis.Z, -(Pi + p_a)) * Matrix(Matrix.Axis.X, pi) * Matrix(Matrix.Axis.Z, Pi)

        //     return RectangularVector(pA, W, z)
    }

//    fun transformPrec(T: Double, R: Vector): Vector {
//        val angles = getEclipticPrecessionVector(T)
//        var eps = Obliquity.meanObliquity(0.0)
//        val x = doubleArrayOf(0.0, 0.0, 0.0)
//        var z = 0.0;
//
//        x[0] = R[0]
//        x[1] = R[1]
//        x[2] = R[2]
//
////        x[0] = R[0];
////        var z = Math.cos(eps) * R[1] + Math.sin(eps) * R[2];
////        x[2] = -Math.sin(eps) * R[1] + Math.cos(eps) * R[2];
////        x[1] = z
//
//        /*
//		 * Precession in longitude
//		 */
//        val T10 = T / 10.0; /* thousands of years */
//        val pA = angles[0];
//        /*
//		 * Node of the moving ecliptic on the JD0 ecliptic.
//		 */
//        var W = angles[1];
//
//        /*
//         * Rotate about z axis to the node.
//         */
//        z = W;
//        var B = Math.cos(z);
//        var A = Math.sin(z);
//        z = B * x[0] + A * x[1];
//        x[1] = -A * x[0] + B * x[1];
//        x[0] = z;
//
//        /*
//         * Rotate about new x axis by the inclination of the moving ecliptic on
//         * the JD0 ecliptic.
//         */
//        z = angles[2];
//
//        B = Math.cos(z);
//        A = Math.sin(z);
//        z = B * x[1] + A * x[2];
//        x[2] = -A * x[1] + B * x[2];
//        x[1] = z;
//
//        /*
//         * Rotate about new z axis back from the node.
//         */
//        z = -W - pA;
//        B = Math.cos(z);
//        A = Math.sin(z);
//        z = B * x[0] + A * x[1];
//        x[1] = -A * x[0] + B * x[1];
//        x[0] = z;
//
//        /*
//         * Rotate about x axis to final equator.
//         */
////        eps = Obliquity.meanObliquity(T);
////
////        z = Math.cos(eps) * x[1] - Math.sin(eps) * x[2];
////        x[2] = Math.sin(eps) * x[1] + Math.cos(eps) * x[2];
////        x[1] = z;
//        return RectangularVector(x)
//    }

    /**
     * Precess rectangular equatorial coordinates from J2000 epoch.
     * @param JD Equinox of the output in Julian day (TT).
     * @param R Array with x, y, z.
     * @param eph Ephemeris properties.
     * @return Array with corrected x, y, z.
     */
    @Equatorial
    fun precessFromJ2000(T: Double, @Equatorial R: Vector): Vector {
        return getEquatorialPrecessionMatrix(T, true) * R
    }

    /**
     * Precess rectangular equatorial coordinates to J2000 epoch.
     *
     * @param JD Equinox of the input in Julian day (TT).
     * @param R Array with x, y, z.
     * @param eph The ephemeris properties.
     * @return Array with corrected x, y, z.
     */
    @Equatorial
    fun precessToJ2000(T: Double, @Equatorial R: Vector): Vector {
        return getEquatorialPrecessionMatrix(T, false) * R
    }


    /**
     * Precess rectangular equatorial coordinates between two epochs.
     *
     * @param fromDatejT Equinox of the input in Julian day (TT).
     * @param toDatejT Equinox of the output in Julian day (TT).
     * @param R Array with x, y, z.
     * @return Array with corrected x, y, z.
     */
    @Equatorial
    fun precess(fromDatejT: Double, toDatejT: Double, R: Vector): Vector {
        return getEquatorialPrecessionMatrix(fromDatejT, toDatejT) * R
    }


}