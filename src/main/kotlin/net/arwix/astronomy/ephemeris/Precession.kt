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
import net.arwix.astronomy.core.PI2
import net.arwix.astronomy.core.vector.Matrix
import net.arwix.astronomy.core.vector.RectangularVector
import net.arwix.astronomy.core.vector.Vector


sealed class Precession(val t: Double) : PrecessionTransformation {

    /**
     * Precession following Vondrak et al. 2011. See A&amp;A 534, A22.
     *
     * @param t Julian centuries of input vector (equatorial).
     * @return PrecessionTransformation referred to mean equinox and equator of JD.
     */
    @Equatorial class Vondrak2011(t: Double) : Precession(t), PrecessionTransformation by Impl(false, VonrakMatrix.getMatrix(t))

    /**
     * Same as IAU2006, but planetary rotation models are those recommended by
     * the IAU working group on carthographic coordinates, in 2009.
     *
     * @param t Julian centuries of input vector (equatorial).
     * @return PrecessionTransformation referred to mean equinox and equator of JD.
     */
    @Equatorial class IAU2009(t: Double) : Precession(t), PrecessionTransformation by Impl(false, getIAU2006Matrix(t))

    /**
     * Precession following Capitaine et al. 2003.
     *
     * Capitaine formula of precession is to be officially adopted by the IAU,
     * see recommendation in the report of the IAU Division I Working Group on
     * Precession and the Ecliptic (Hilton et al. 2006, Celest. Mech., 94,
     * 351-367).
     * Reference: Capitaine et al., Astronomy & Astrophysics 412, 567-586,
     * 2003.
     *
     * @param t Julian centuries of input vector (equatorial).
     * @return PrecessionTransformation referred to mean equinox and equator of JD.
     */
    @Equatorial class IAU2006(t: Double) : Precession(t), PrecessionTransformation by Impl(false, getIAU2006Matrix(t))

    /**
     * Precession following IAU2000 definitions. From SOFA software library.
     * Reference: Capitaine et al., Astronomy & Astrophysics 400, 1145-1154,
     * 2003. See also Lieske et al. 1977.
     *
     * @param t Julian centuries of input vector (equatorial).
     * @return PrecessionTransformation referred to mean equinox and equator of JD.
     */
    @Equatorial class IAU2000(t: Double) : Precession(t), PrecessionTransformation by Impl(false, getIAU2000Matrix(t, true), getIAU2000Matrix(t, false))

    /**
     * Precession for selecting JPL DE403/404/405/406 formulae for precession,
     * obliquity, nutation (IAU 1980), and Greenwich mean sidereal time. Quite
     * similar to Williams formulae. Adequate for planets using Moshier
     * method, Series96, or JPL DE40x ephemerides.
     *
     * @param t Julian centuries of input vector (ecliptic).
     * @return PrecessionTransformation referred to mean equinox and ecliptic of JD.
     */
    @Ecliptic class DE4xx(val tt: Double) : Precession(tt), PrecessionTransformation by Impl(true, getEclipticMatrix(JPL_DE4xx_Matrices, tt)) {

    }

    /**
     * Precession for selecting SIMON formulae of precession, obliquity,
     * nutation (IAU 1980), and Greenwich mean sidereal time. See
     * J. L. Simon, P. Bretagnon, J. Chapront, M. Chapront-Touze', G. Francou,
     * and J. Laskar, "Numerical Expressions for precession formulae and mean
     * elements for the Moon and the planets," Astronomy and Astrophysics 282,
     * 663-683 (1994).
     *
     * @param t Julian centuries of input vector (ecliptic).
     * @return PrecessionTransformation referred to mean equinox and ecliptic of JD.
     */
    @Ecliptic class Simon1994(t: Double) : Precession(t), PrecessionTransformation by Impl(true, getEclipticMatrix(SIMON_1994_Matrices, t))

    /**
     * Precession for selecting Williams formulae of precession (DE403 JPL
     * Ephemeris), nutation (IAU 1980), obliquity, and Greenwich mean sidereal
     * time. See James G. Williams, "Contributions to the Earth's obliquity rate,
     * precession, and nutation," Astron. J. 108, 711-724 (1994). It is convenient
     * to use this when obtaining ephemeris of the Moon using Moshier method.
     *
     * @param t Julian centuries of input vector (ecliptic).
     * @return PrecessionTransformation referred to mean equinox and ecliptic of JD.
     */
    @Ecliptic class Williams1994(t: Double) : Precession(t), PrecessionTransformation by Impl(true, getEclipticMatrix(WILLIAMS_1994_Matrices, t))

    /**
     * Precession for selecting Laskar formulae of precession, nutation (IAU
     * 1980), and Greenwich mean sidereal time. See J. Laskar,
     * "Secular terms of classical planetary theories using the results of
     * general theory," Astronomy and Astrophysics 157, 59070 (1986).
     *
     * @param t Julian centuries of input vector (ecliptic).
     * @return PrecessionTransformation referred to mean equinox and ecliptic of JD.
     */
    @Ecliptic class Laskar1986(t: Double) : Precession(t), PrecessionTransformation by Impl(true, getEclipticMatrix(LASKAR_1986_Matrices, t))

    /**
     * Precession for selecting IAU 1976 formulae of precession, nutation (IAU
     * 1980), and Greenwich mean sidereal time. This will use
     * old formulae that will match results from IMCCE ephemeris server. You
     * may consider using this formulae for VSOP82 theory. See
     * J. H. Lieske, T. Lederle, W. Fricke, and B. Morando, "Expressions for the
     * Precession Quantities Based upon the IAU (1976) System of Astronomical
     * Constants," Astronomy and Astrophysics 58, 1-16 (1977).
     *
     * @param t Julian centuries of input vector (ecliptic).
     * @return PrecessionTransformation referred to mean equinox and ecliptic of JD.
     */
    @Ecliptic class IAU1976(t: Double) : Precession(t), PrecessionTransformation by Impl(true, getEclipticMatrix(IAU_1976_Matrices, t))


    fun getNearestObliquityModel(t: Double = this.t) =
            when (this) {
                is Precession.Vondrak2011 -> Obliquity.Vondrak2011(t)
                is Precession.IAU2009,
                is Precession.IAU2006,
                is Precession.IAU2000 -> Obliquity.IAU2006(t)
                is Precession.Williams1994,
                is Precession.DE4xx -> Obliquity.Williams1994(t)
                is Precession.Simon1994 -> Obliquity.Simon1994(t)
                is Precession.Laskar1986 -> Obliquity.Laskar1996(t)
                is Precession.IAU1976 -> Obliquity.IAU1976(t)

            }

    fun getNearsetNutationModel(obliquity: Obliquity, t: Double = this.t) = when (this) {
        is Precession.IAU2009,
        is Precession.IAU2006 -> Nutation.IAU2006(t, obliquity)
        is Precession.Williams1994,
        is Precession.DE4xx,
        is Precession.Simon1994,
        is Precession.Laskar1986,
        is Precession.IAU1976 -> Nutation.IAU1980(t, obliquity)
        else -> Nutation.IAU2000(t, obliquity)
    }

}

interface PrecessionTransformation {
    val isEcliptic: Boolean
    val fromJ2000Matrix: Matrix
    val toJ2000Matrix: Matrix
    fun transformFromJ2000(vector: Vector): Vector
    fun transformToJ2000(vector: Vector): Vector
}

private class Impl(
        override val isEcliptic: Boolean,
        override val fromJ2000Matrix: Matrix,
        override val toJ2000Matrix: Matrix = fromJ2000Matrix.transpose()) : PrecessionTransformation {

    override fun transformFromJ2000(vector: Vector) = fromJ2000Matrix * vector
    override fun transformToJ2000(vector: Vector) = toJ2000Matrix * vector
}

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

private fun getEclipticMatrix(list: Array<DoubleArray>, T: Double, isFromJ2000ToApparent: Boolean = true): Matrix {
    val T10 = T / 10.0 /* thousands of years */
    val pA = list[0].fold(0.0, { acc, d -> acc * T10 + d }) * ARCSEC_TO_RAD * T10
    val W = list[1].fold(0.0, { acc, d -> acc * T10 + d })
    val z = list[2].fold(0.0, { acc, d -> acc * T10 + d }).let { if (!isFromJ2000ToApparent) -it else it }

    return Matrix(Matrix.Axis.Z, -(W + pA)) * Matrix(Matrix.Axis.X, z) * Matrix(Matrix.Axis.Z, W)
}


private fun getIAU2000Matrix(T: Double, isFromJ2000ToApparent: Boolean): Matrix {
    val T0: Double = if (isFromJ2000ToApparent) T else 0.0
    var EPS0 = 84381.448
    val (PSIA, OMEGAA, CHIA) = RectangularVector(
            ((((-0.0 * T + 0.0) * T - 0.001147) * T - 1.07259) * T + 5038.7784) * T - 0.29965 * T0,
            ((((+0.0 * T - 0.0) * T - 0.007726) * T + 0.05127) * T - 0.0) * T + EPS0 - 0.02524 * T0,
            ((((-0.0 * T + 0.0) * T - 0.001125) * T - 2.38064) * T + 10.5526) * T) * ARCSEC_TO_RAD

    EPS0 *= ARCSEC_TO_RAD

    return Matrix(Matrix.Axis.Z, CHIA) * Matrix(Matrix.Axis.X, -OMEGAA) * Matrix(Matrix.Axis.Z, -PSIA) * Matrix(Matrix.Axis.X, EPS0)
}


private fun getIAU2006Matrix(T: Double): Matrix {
    var EPS0 = 84381.406
    val (PSIA, OMEGAA, CHIA) = RectangularVector(
            ((((-0.0000000951 * T + 0.000132851) * T - 0.00114045) * T - 1.0790069) * T + 5038.481507) * T,
            ((((+0.0000003337 * T - 0.000000467) * T - 0.00772503) * T + 0.0512623) * T - 0.025754) * T + EPS0,
            ((((-0.0000000560 * T + 0.000170663) * T - 0.00121197) * T - 2.3814292) * T + 10.556403) * T) * ARCSEC_TO_RAD

    EPS0 *= ARCSEC_TO_RAD

    return Matrix(Matrix.Axis.Z, CHIA) * Matrix(Matrix.Axis.X, -OMEGAA) * Matrix(Matrix.Axis.Z, -PSIA) * Matrix(Matrix.Axis.X, EPS0)
}


private object VonrakMatrix {

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

    fun getMatrix(T: Double): Matrix {
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