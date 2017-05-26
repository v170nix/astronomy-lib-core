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

package net.arwix.astronomy.core.kepler

import net.arwix.astronomy.core.ARCSEC_TO_DEG
import net.arwix.astronomy.core.DEG_TO_RAD
import net.arwix.astronomy.core.GM_Sun
import net.arwix.astronomy.math.polynomialSum
import net.arwix.astronomy.math.radians.Radian
import java.lang.Math.cos
import java.lang.Math.sin

private fun Double.deg() = this * ARCSEC_TO_DEG

sealed class KeplerBodySimonJ2000 : KeplerElements {

    object Mercury : KeplerBodySimonJ2000(),
            KeplerElements by ImplKeplerElements(
                    aCoefficients = doubleArrayOf(0.3870983098),
                    LCoefficients = doubleArrayOf(252.25090552 - 0.047.deg(), 5381016286.88982.deg(), -1.92789.deg(), +0.00639.deg()),
                    eCoefficients = doubleArrayOf(0.2056317526, 0.0002040653, -28349e-10, -1805e-10, 23e-10, -2e-10),
                    WCoefficients = doubleArrayOf(77.45611904, 5719.1159.deg(), -4.83016.deg(), -0.02464.deg(), -0.00016.deg(), 0.00004.deg()),
                    iCoefficients = doubleArrayOf(7.00498625, -214.25629.deg(), 0.28977.deg(), 0.15421.deg(), -0.00169.deg(), -0.00002.deg()),
                    OCoefficients = doubleArrayOf(48.33089304, -4515.21727.deg(), -31.79892.deg(), -0.71933.deg(), 0.01242.deg()), k = 10.0)

    object Venus : KeplerBodySimonJ2000(),
            KeplerElements by ImplKeplerElements(
                    aCoefficients = doubleArrayOf(0.7233298200),
                    LCoefficients = doubleArrayOf(181.97980085, 2106641364.33548.deg(), 0.59381.deg(), -0.00627.deg()),
                    eCoefficients = doubleArrayOf(0.0067719164, -0.0004776521, 98127e-10, 4639e-10, 123e-10, -3e-10),
                    WCoefficients = doubleArrayOf(131.56370300, 175.48640.deg(), -498.48184.deg(), -20.50042.deg(), -0.72432.deg(), 0.00224.deg()),
                    iCoefficients = doubleArrayOf(3.39466189, -30.84437.deg(), -11.67836.deg(), 0.03338.deg(), 0.00269.deg(), 0.00004.deg()),
                    OCoefficients = doubleArrayOf(76.67992019, -10008.48154.deg(), -51.32614.deg(), -0.5891.deg(), -0.004665.deg()), k = 10.0)


    object Earth : KeplerBodySimonJ2000(),
            KeplerElements by ImplKeplerElements(
                    aCoefficients = doubleArrayOf(1.0000010178),
                    LCoefficients = doubleArrayOf(100.46645683, 1295977422.83429.deg(), -2.04411.deg(), -0.00523.deg()),
                    eCoefficients = doubleArrayOf(0.0167086342, -0.0004203654, -0.0000126734, 1444e-10, -2e-10, 3e-10),
                    WCoefficients = doubleArrayOf(102.93734808, 11612.35290.deg(), 53.27577.deg(), -0.14095.deg(), 0.11440.deg(), 0.00478.deg()),
                    iCoefficients = doubleArrayOf(0.0, 469.97289.deg(), -3.35053.deg(), -0.12374.deg(), 0.00027.deg(), -0.00001.deg(), 0.00001.deg()),
                    OCoefficients = doubleArrayOf(174.87317577, -8679.27034.deg(), 15.34191.deg(), 0.00532.deg(), -0.03734.deg(), -0.00073.deg(), 0.00004.deg()), k = 10.0)


    object Mars : KeplerBodySimonJ2000(),
            KeplerElements by ImplKeplerElements(

                    aCoefficients = doubleArrayOf(1.5236793419, 3e-10),
                    LCoefficients = doubleArrayOf(355.43299958, 689050774.93988.deg(), 0.94264.deg(), -0.01043.deg()),
                    eCoefficients = doubleArrayOf(0.0934006477, 0.0009048438, -80641e-10, -2519e-10, 124e-10, -10e-10),
                    WCoefficients = doubleArrayOf(336.06023395, 15980.45908.deg(), -62.32800.deg(), 1.86464.deg(), -0.04603.deg(), -0.00164.deg()),
                    iCoefficients = doubleArrayOf(1.84972648, -293.31722.deg(), -8.11830.deg(), -0.10326.deg(), -0.00153.deg(), 0.00048.deg()),
                    OCoefficients = doubleArrayOf(49.55809321, -10620.90088.deg(), -230.57416.deg(), -7.06942.deg(), -0.6892.deg(), -0.05829.deg()), k = 10.0)


    object Jupiter : KeplerBodySimonJ2000(),
            KeplerElements by ImplKeplerElements(

                    aCoefficients = doubleArrayOf(5.2026032092, 19132e-10, -39e-10, -60e-10, -10e-10, 1e-10),
                    LCoefficients = doubleArrayOf(34.35151874, 109256603.77991.deg(), -30.60378.deg(), 0.05706.deg(), 0.04667.deg(), 0.00591.deg(), -0.00034.deg()),
                    eCoefficients = doubleArrayOf(0.0484979255, 0.0016322542, -0.0000471366, -20063e-10.deg(), 1018e-10.deg(), -21e-10.deg(), 1e10.deg()),
                    WCoefficients = doubleArrayOf(14.33120687, 7758.75163.deg(), 259.95938.deg(), -16.14731.deg(), 0.74704.deg(), -0.02087.deg(), -0.00016.deg()),
                    iCoefficients = doubleArrayOf(1.30326698, -71.55890.deg(), 11.95297.deg(), 0.340909.deg(), -0.02710.deg(), -0.00124.deg(), 0.00003.deg()),
                    OCoefficients = doubleArrayOf(100.46440702, 6362.03561.deg(), 326.52178.deg(), -26.18091.deg(), -2.10322.deg(), 0.04453.deg(), 0.01154.deg()), k = 10.0)


    object Saturn : KeplerBodySimonJ2000(),
            KeplerElements by ImplKeplerElements(
                    aCoefficients = doubleArrayOf(9.5549091915, -0.0000213896, 444e-10, 670e-10, 110e-10, -7e-10, -1e-10),
                    LCoefficients = doubleArrayOf(50.07744430, 43996098.55732.deg(), 75.61614.deg(), -0.16618.deg(), -0.11484.deg(), -0.01452.deg(), 0.00083.deg()),
                    eCoefficients = doubleArrayOf(0.0555481426, -0.0034664062, -0.0000643639, 33956e-10, -219e-10, -3e-10, 6e-10),
                    WCoefficients = doubleArrayOf(93.05723748, 20395.49439.deg(), 190.25952.deg(), 17.68303.deg(), 1.23148.deg(), 0.10310.deg(), 0.00702.deg()),
                    iCoefficients = doubleArrayOf(2.48887878, 91.85195.deg(), -17.66225.deg(), 0.06105.deg(), 0.02638.deg(), -0.00152.deg(), -0.00012.deg()),
                    OCoefficients = doubleArrayOf(113.66550252, -9240.19942.deg(), -66.23743.deg(), 1.72778.deg(), 0.2699.deg(), 0.03610.deg(), -0.00248.deg()), k = 10.0)


    object Uranus : KeplerBodySimonJ2000(),
            KeplerElements by ImplKeplerElements(
                    aCoefficients = doubleArrayOf(19.2184460618, -3716e-10, 979e-10),
                    LCoefficients = doubleArrayOf(314.05500511, 15424811.93933.deg(), -1.75083.deg(), 0.02156.deg()),
                    eCoefficients = doubleArrayOf(0.0463812221, -0.0002729293, 0.0000078913, 2447e-10.deg(), -171e-10.deg()),
                    WCoefficients = doubleArrayOf(173.00529106, 3215.56238.deg(), -34.09288.deg(), 1.48909.deg(), 0.066.deg()),
                    iCoefficients = doubleArrayOf(0.77319689, -60.72723.deg(), 1.25759.deg(), 0.05808.deg(), 0.00031.deg()),
                    OCoefficients = doubleArrayOf(74.00595701, 2669.15033.deg(), 145.93964.deg(), 0.42917.deg(), -0.0912.deg()), k = 10.0)


    object Neptune : KeplerBodySimonJ2000(),
            KeplerElements by ImplKeplerElements(

                    aCoefficients = doubleArrayOf(30.1103868694, -16635e-10, 686e-10),
                    LCoefficients = doubleArrayOf(304.34866548, 7865503.20744.deg(), 0.21103.deg(), -0.00895.deg()),
                    eCoefficients = doubleArrayOf(0.0094557470, 0.0000603263, 0.0, -483e-10.deg()),
                    WCoefficients = doubleArrayOf(48.12027554, 1050.71912.deg(), 27.39717.deg()),
                    iCoefficients = doubleArrayOf(1.76995259, 8.12333.deg(), 0.08135.deg(), -0.00046.deg()),
                    OCoefficients = doubleArrayOf(131.78405702, -221.94322.deg(), -0.78728.deg(), -0.28070.deg(), 0.00049.deg()), k = 10.0)

}

/**
 * Keplerian elements and their rates, with respect to the mean ecliptic and equinox of J2000,
 * valid for the time-interval 3000 BC -- 3000 AD.
 */
sealed class KeplerBodyJPL : KeplerElements {

    //https://ssd.jpl.nasa.gov/?planets#ephem
    //https://ssd.jpl.nasa.gov/txt/p_elem_t2.txt
    //https://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf
    object Mercury : KeplerBodyJPL(),
            KeplerElements by ImplKeplerElements(

                    aCoefficients = doubleArrayOf(0.38709843),
                    eCoefficients = doubleArrayOf(0.20563661, 0.00002123),
                    iCoefficients = doubleArrayOf(7.00559432, -0.00590158),
                    LCoefficients = doubleArrayOf(252.25166724, 149472.67486623),
                    WCoefficients = doubleArrayOf(77.45771895, 0.15940013),
                    OCoefficients = doubleArrayOf(48.33961819, -0.12214182))

    object Venus : KeplerBodyJPL(),
            KeplerElements by ImplKeplerElements(

                    aCoefficients = doubleArrayOf(0.72332102, -0.00000026),
                    eCoefficients = doubleArrayOf(0.00676399, -0.00005107),
                    iCoefficients = doubleArrayOf(3.39777545, 0.00043494),
                    LCoefficients = doubleArrayOf(181.97970850, 58517.81560260),
                    WCoefficients = doubleArrayOf(131.76755713, 0.05679648),
                    OCoefficients = doubleArrayOf(76.67261496, -0.27274174))

    object EarthMoonBarycenter : KeplerBodyJPL(),
            KeplerElements by ImplKeplerElements(

                    aCoefficients = doubleArrayOf(1.00000018, -0.00000003),
                    eCoefficients = doubleArrayOf(0.01673163, -0.00003661),
                    iCoefficients = doubleArrayOf(-0.00054346, -0.01337178),
                    LCoefficients = doubleArrayOf(100.46691572, 35999.37306329),
                    WCoefficients = doubleArrayOf(102.93005885, 0.31795260),
                    OCoefficients = doubleArrayOf(-5.11260389, -0.24123856))

    object Mars : KeplerBodyJPL(),
            KeplerElements by ImplKeplerElements(

                    aCoefficients = doubleArrayOf(1.52371243, 0.00000097),
                    eCoefficients = doubleArrayOf(0.09336511, 0.00009149),
                    iCoefficients = doubleArrayOf(1.85181869, -0.00724757),
                    LCoefficients = doubleArrayOf(-4.56813164, 19140.29934243),
                    WCoefficients = doubleArrayOf(-23.91744784, 0.45223625),
                    OCoefficients = doubleArrayOf(49.71320984, -0.26852431))

    object Jupiter : KeplerBodyJPL(),
            KeplerElements by ImplKeplerElements(


                    aCoefficients = doubleArrayOf(5.20248019, -0.00002864),
                    eCoefficients = doubleArrayOf(0.04853590, 0.00018026),
                    iCoefficients = doubleArrayOf(1.29861416, -0.00322699),
                    LCoefficients = doubleArrayOf(34.33479152, 3034.90371757),
                    WCoefficients = doubleArrayOf(14.27495244, 0.18199196),
                    OCoefficients = doubleArrayOf(100.29282654, 0.13024619)) {
        private val f = 38.35125000 * DEG_TO_RAD
        private val c = 0.06064060
        private val s = -0.35635438


        override fun getMeanAnomaly(T: Double) = getLongitude(T) - getPerihelionLongitude(T) - 0.00012452 * T * T + c * cos(f * T) + s * sin(f * T)
    }

    object Saturn : KeplerBodyJPL(),
            KeplerElements by ImplKeplerElements(


                    doubleArrayOf(9.54149883, -0.00003065),
                    doubleArrayOf(0.05550825, -0.00032044),
                    doubleArrayOf(2.49424102, 0.00451969),
                    doubleArrayOf(50.07571329, 1222.11494724),
                    doubleArrayOf(92.86136063, 0.54179478),
                    doubleArrayOf(113.63998702, -0.25015002)) {
        private val f = 38.35125000 * DEG_TO_RAD
        private val c = -0.13434469
        private val s = 0.87320147


        override fun getMeanAnomaly(T: Double) = getLongitude(T) - getPerihelionLongitude(T) + 0.00025899 * T * T + c * cos(f * T) + s * sin(f * T)
    }

    object Uranus : KeplerBodyJPL(), KeplerElements by ImplKeplerElements(
            doubleArrayOf(19.18797948, -0.00020455),
            doubleArrayOf(0.04685740, -0.00001550),
            doubleArrayOf(0.77298127, -0.00180155),
            doubleArrayOf(314.20276625, 428.49512595),
            doubleArrayOf(172.43404441, 0.09266985),
            doubleArrayOf(73.96250215, 0.05739699)) {
        private val f = 7.67025000 * DEG_TO_RAD
        private val c = -0.97731848
        private val s = 0.17689245

        override fun getMeanAnomaly(T: Double) = getLongitude(T) - getPerihelionLongitude(T) + 0.00058331 * T * T + c * cos(f * T) + s * sin(f * T)
    }

    object Neptune : KeplerBodyJPL(), KeplerElements by ImplKeplerElements(
            doubleArrayOf(30.06952752, 0.00006447),
            doubleArrayOf(0.00895439, 0.00000818),
            doubleArrayOf(1.77005520, 0.00022400),
            doubleArrayOf(304.22289287, 218.46515314),
            doubleArrayOf(46.68158724, 0.01009938),
            doubleArrayOf(131.78635853, -0.00606302)) {
        private val f = 7.67025000 * DEG_TO_RAD
        private val c = 0.68346318
        private val s = -0.10162547

        override fun getMeanAnomaly(T: Double) = getLongitude(T) - getPerihelionLongitude(T) - 0.00041348 * T * T + c * cos(f * T) + s * sin(f * T)
    }

    object Pluto : KeplerBodyJPL(), KeplerElements by ImplKeplerElements(
            doubleArrayOf(39.48686035, 0.00449751),
            doubleArrayOf(0.24885238, 0.00006016),
            doubleArrayOf(17.14104260, 0.00000501),
            doubleArrayOf(238.96535011, 145.18042903),
            doubleArrayOf(224.09702598, -0.00968827),
            doubleArrayOf(110.30167986, -0.00809981)) {

        override fun getMeanAnomaly(T: Double) = getLongitude(T) - getPerihelionLongitude(T) - 0.01262724 * T * T
    }

}

private class ImplKeplerElements(
        //    override val T: Double,
        private val aCoefficients: DoubleArray,
        private val eCoefficients: DoubleArray,
        private val iCoefficients: DoubleArray,
        private val LCoefficients: DoubleArray,
        private val WCoefficients: DoubleArray,
        private val OCoefficients: DoubleArray,
        private val k: Double = 1.0

) : KeplerElements {
    override fun getSemiMajorAxis(T: Double) = aCoefficients.polynomialSum(T / k)
    override fun getEccentricity(T: Double) = eCoefficients.polynomialSum(T / k)
    override fun getInclination(T: Double) = iCoefficients.polynomialSum(T / k) * DEG_TO_RAD
    override fun getPerihelionLongitude(T: Double) = WCoefficients.polynomialSum(T / k) * DEG_TO_RAD
    override fun getAscendingNodeLongitude(T: Radian) = OCoefficients.polynomialSum(T / k) * DEG_TO_RAD
    override fun getLongitude(T: Double) = LCoefficients.polynomialSum(T / k) * DEG_TO_RAD
}


interface KeplerElements {
    //   val T: Double
//    val a: Double //a
//    val e: Radian // e
//    val inclination: Radian // i
    //  val Longitude: Radian // L
//    val perihelionLongitude: Radian // w
//    val ascendingNodeLongitude: Radian // O

    fun getSemiMajorAxis(T: Double): Double //a
    fun getEccentricity(T: Double): Radian // e
    fun getInclination(T: Double): Radian //i

    fun getLongitude(T: Double): Radian //L
    fun getPerihelionLongitude(T: Double): Radian //w
    fun getAscendingNodeLongitude(T: Radian): Radian // O

    fun getMeanAnomaly(T: Double) = getLongitude(T) - getPerihelionLongitude(T)

    fun getOrbitalPlane(T: Double): OrbitalPlane {

        val p = EarthMoonElements(T).precession
        val M = getMeanAnomaly(T)
        return EllipticOrbit.getOrbitalPlane(GM_Sun, M, getSemiMajorAxis(T), getEccentricity(T)).let {
            val PQR = KeplerOrbit.createGaussianMatrix(getAscendingNodeLongitude(T) + p, getInclination(T), getPerihelionLongitude(T) - getAscendingNodeLongitude(T))
            OrbitalPlane(PQR * it.position, PQR * it.velocity)
        }
    }
}