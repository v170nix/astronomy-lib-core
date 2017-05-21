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

import net.arwix.astronomy.core.DEG_TO_RAD
import net.arwix.astronomy.core.GM_Sun
import net.arwix.astronomy.math.polynomialSum
import net.arwix.astronomy.math.radians.Radian
import java.lang.Math.cos
import java.lang.Math.sin

/**
 * Keplerian elements and their rates, with respect to the mean ecliptic and equinox of J2000,
 * valid for the time-interval 3000 BC -- 3000 AD.
 */
sealed class JPLKeplerObject : KeplerElements {

    //https://ssd.jpl.nasa.gov/?planets#ephem
    //https://ssd.jpl.nasa.gov/txt/p_elem_t2.txt
    //https://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf
    class Mercury(override val T: Double) : JPLKeplerObject(),
            KeplerElements by ImplKeplerElements(T, aCoefficients, eCoefficients, iCoefficients, LCoefficients, WCoefficients, OCoefficients) {
        companion object {
            private val aCoefficients = doubleArrayOf(0.38709843)
            private val eCoefficients = doubleArrayOf(0.20563661, 0.00002123)
            private val iCoefficients = doubleArrayOf(7.00559432, -0.00590158)
            private val LCoefficients = doubleArrayOf(252.25166724, 149472.67486623)
            private val WCoefficients = doubleArrayOf(77.45771895, 0.15940013)
            private val OCoefficients = doubleArrayOf(48.33961819, -0.12214182)
        }
    }

    class Venus(override val T: Double) : JPLKeplerObject(),
            KeplerElements by ImplKeplerElements(T, aCoefficients, eCoefficients, iCoefficients, LCoefficients, WCoefficients, OCoefficients) {
        companion object {
            private val aCoefficients = doubleArrayOf(0.72332102, -0.00000026)
            private val eCoefficients = doubleArrayOf(0.00676399, -0.00005107)
            private val iCoefficients = doubleArrayOf(3.39777545, 0.00043494)
            private val LCoefficients = doubleArrayOf(181.97970850, 58517.81560260)
            private val WCoefficients = doubleArrayOf(131.76755713, 0.05679648)
            private val OCoefficients = doubleArrayOf(76.67261496, -0.27274174)
        }
    }

    class EarthMoonBarycenter(override val T: Double) : JPLKeplerObject(),
            KeplerElements by ImplKeplerElements(T, aCoefficients, eCoefficients, iCoefficients, LCoefficients, WCoefficients, OCoefficients) {
        companion object {
            private val aCoefficients = doubleArrayOf(1.00000018, -0.00000003)
            private val eCoefficients = doubleArrayOf(0.01673163, -0.00003661)
            private val iCoefficients = doubleArrayOf(-0.00054346, -0.01337178)
            private val LCoefficients = doubleArrayOf(100.46691572, 35999.37306329)
            private val WCoefficients = doubleArrayOf(102.93005885, 0.31795260)
            private val OCoefficients = doubleArrayOf(-5.11260389, -0.24123856)
        }
    }

    class Mars(override val T: Double) : JPLKeplerObject(),
            KeplerElements by ImplKeplerElements(T, aCoefficients, eCoefficients, iCoefficients, LCoefficients, WCoefficients, OCoefficients) {
        companion object {
            private val aCoefficients = doubleArrayOf(1.52371243, 0.00000097)
            private val eCoefficients = doubleArrayOf(0.09336511, 0.00009149)
            private val iCoefficients = doubleArrayOf(1.85181869, -0.00724757)
            private val LCoefficients = doubleArrayOf(-4.56813164, 19140.29934243)
            private val WCoefficients = doubleArrayOf(-23.91744784, 0.45223625)
            private val OCoefficients = doubleArrayOf(49.71320984, -0.26852431)
        }
    }

    class Jupiter(override val T: Double) : JPLKeplerObject(),
            KeplerElements by ImplKeplerElements(T, aCoefficients, eCoefficients, iCoefficients, LCoefficients, WCoefficients, OCoefficients) {

        companion object {
            private val aCoefficients = doubleArrayOf(5.20248019, -0.00002864)
            private val eCoefficients = doubleArrayOf(0.04853590, 0.00018026)
            private val iCoefficients = doubleArrayOf(1.29861416, -0.00322699)
            private val LCoefficients = doubleArrayOf(34.33479152, 3034.90371757)
            private val WCoefficients = doubleArrayOf(14.27495244, 0.18199196)
            private val OCoefficients = doubleArrayOf(100.29282654, 0.13024619)
            private val f = 38.35125000 * DEG_TO_RAD
            private val c = 0.06064060
            private val s = -0.35635438
        }

        override fun getMeanAnomaly() = Longitude - perihelionLongitude - 0.00012452 * T * T + c * cos(f * T) + s * sin(f * T)
    }

    class Saturn(override val T: Double) : JPLKeplerObject(),
            KeplerElements by ImplKeplerElements(T, aCoefficients, eCoefficients, iCoefficients, LCoefficients, WCoefficients, OCoefficients) {

        companion object {
            private val aCoefficients = doubleArrayOf(9.54149883, -0.00003065)
            private val eCoefficients = doubleArrayOf(0.05550825, -0.00032044)
            private val iCoefficients = doubleArrayOf(2.49424102, 0.00451969)
            private val LCoefficients = doubleArrayOf(50.07571329, 1222.11494724)
            private val WCoefficients = doubleArrayOf(92.86136063, 0.54179478)
            private val OCoefficients = doubleArrayOf(113.63998702, -0.25015002)
            private val f = 38.35125000 * DEG_TO_RAD
            private val c = -0.13434469
            private val s = 0.87320147
        }

        override fun getMeanAnomaly() = Longitude - perihelionLongitude + 0.00025899 * T * T + c * cos(f * T) + s * sin(f * T)
    }

    class Uranus(override val T: Double) : JPLKeplerObject(),
            KeplerElements by ImplKeplerElements(T, aCoefficients, eCoefficients, iCoefficients, LCoefficients, WCoefficients, OCoefficients) {

        companion object {
            private val aCoefficients = doubleArrayOf(19.18797948, -0.00020455)
            private val eCoefficients = doubleArrayOf(0.04685740, -0.00001550)
            private val iCoefficients = doubleArrayOf(0.77298127, -0.00180155)
            private val LCoefficients = doubleArrayOf(314.20276625, 428.49512595)
            private val WCoefficients = doubleArrayOf(172.43404441, 0.09266985)
            private val OCoefficients = doubleArrayOf(73.96250215, 0.05739699)
            private val f = 7.67025000 * DEG_TO_RAD
            private val c = -0.97731848
            private val s = 0.17689245
        }

        override fun getMeanAnomaly() = Longitude - perihelionLongitude + 0.00058331 * T * T + c * cos(f * T) + s * sin(f * T)
    }

    class Neptune(override val T: Double) : JPLKeplerObject(),
            KeplerElements by ImplKeplerElements(T, aCoefficients, eCoefficients, iCoefficients, LCoefficients, WCoefficients, OCoefficients) {

        companion object {
            private val aCoefficients = doubleArrayOf(30.06952752, 0.00006447)
            private val eCoefficients = doubleArrayOf(0.00895439, 0.00000818)
            private val iCoefficients = doubleArrayOf(1.77005520, 0.00022400)
            private val LCoefficients = doubleArrayOf(304.22289287, 218.46515314)
            private val WCoefficients = doubleArrayOf(46.68158724, 0.01009938)
            private val OCoefficients = doubleArrayOf(131.78635853, -0.00606302)
            private val f = 7.67025000 * DEG_TO_RAD
            private val c = 0.68346318
            private val s = -0.10162547
        }

        override fun getMeanAnomaly() = Longitude - perihelionLongitude - 0.00041348 * T * T + c * cos(f * T) + s * sin(f * T)
    }

    class Pluto(override val T: Double) : JPLKeplerObject(),
            KeplerElements by ImplKeplerElements(T, aCoefficients, eCoefficients, iCoefficients, LCoefficients, WCoefficients, OCoefficients) {

        companion object {
            private val aCoefficients = doubleArrayOf(39.48686035, 0.00449751)
            private val eCoefficients = doubleArrayOf(0.24885238, 0.00006016)
            private val iCoefficients = doubleArrayOf(17.14104260, 0.00000501)
            private val LCoefficients = doubleArrayOf(238.96535011, 145.18042903)
            private val WCoefficients = doubleArrayOf(224.09702598, -0.00968827)
            private val OCoefficients = doubleArrayOf(110.30167986, -0.00809981)
            private val f = 7.67025000 * DEG_TO_RAD
            private val c = 0.68346318
            private val s = -0.10162547
        }

        override fun getMeanAnomaly() = Longitude - perihelionLongitude - 0.01262724 * T * T
    }

}

private class ImplKeplerElements(
        override val T: Double,
        aCoefficients: DoubleArray,
        eCoefficients: DoubleArray,
        iCoefficients: DoubleArray,
        LCoefficients: DoubleArray,
        WCoefficients: DoubleArray,
        OCoefficients: DoubleArray

) : KeplerElements {
    override val a = aCoefficients.polynomialSum(T)
    override val e = eCoefficients.polynomialSum(T)
    override val inclination = iCoefficients.polynomialSum(T) * DEG_TO_RAD
    override val Longitude = LCoefficients.polynomialSum(T) * DEG_TO_RAD
    override val perihelionLongitude = WCoefficients.polynomialSum(T) * DEG_TO_RAD
    override val ascendingNodeLongitude = OCoefficients.polynomialSum(T) * DEG_TO_RAD
}


interface KeplerElements {
    val T: Double
    val a: Double //a
    val e: Radian // e
    val inclination: Radian // i
    val Longitude: Radian // L
    val perihelionLongitude: Radian // w
    val ascendingNodeLongitude: Radian // O

    fun getMeanAnomaly() = Longitude - perihelionLongitude

    fun getOrbitalPlane(): OrbitalPlane {
        val p = 1.3970 * DEG_TO_RAD // прецессия Земли за сто лет
        val M = getMeanAnomaly()
        return EllipticOrbit.getOrbitalPlane(GM_Sun, M, a, e).let {
            val PQR = KeplerianOrbit.createGaussianMatrix(ascendingNodeLongitude + p * T, inclination, perihelionLongitude - ascendingNodeLongitude)
            OrbitalPlane(PQR * it.position, PQR * it.velocity)
        }
    }
}