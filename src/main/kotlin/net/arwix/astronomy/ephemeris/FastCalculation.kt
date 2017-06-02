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
import net.arwix.astronomy.core.AU
import net.arwix.astronomy.core.DEG_TO_RAD
import net.arwix.astronomy.core.coordinates.FunGetGeocentricEclipticCoordinates
import net.arwix.astronomy.core.kepler.KeplerBodySimonJ2000
import net.arwix.astronomy.core.vector.SphericalVector
import net.arwix.astronomy.math.Degrees.Degree
import net.arwix.astronomy.math.Degrees.normalizeDegree
import net.arwix.astronomy.math.Degrees.toRad
import net.arwix.astronomy.math.polynomialSum
import net.arwix.astronomy.math.radians.Radian
import net.arwix.astronomy.math.radians.normalize
import net.arwix.astronomy.math.radians.toDeg
import java.lang.Math.*

object FastCalculation {


    private fun getNutation(t: Double): Radian {
        val A = (124.90 - 1934.134 * t + 0.002063 * t * t).normalizeDegree().toRad()
        val B = (201.11 + 72001.5377 * t + 0.00057 * t * t).normalizeDegree().toRad()
        return (-.004778 * Math.sin(A) - .0003667 * Math.sin(B)).toRad()
    }

    private fun getSolarAnomaly(t: Double): Degree =
            (357.5291092 + 35999.0502909 * t - .0001536 * t * t + 1.0 / 24490000.0 * t * t * t).normalizeDegree()


    // MOON PARAMETERS (Formulae from "Calendrical Calculations") //D
    private fun getLunarElongation(t: Double): Degree =
            doubleArrayOf(297.8501921, 445267.1114034, -0.0018819, 1.0 / 545868.0, -1.0 / 113065000.0)
                    .polynomialSum(t).toRad().normalize().toDeg()

    // Anomalistic phase /ML
    private fun getLunarAnomaly(t: Double): Degree =
            doubleArrayOf(134.9633964, 477198.8675055, .0087414, 1.0 / 69699.0, -1.0 / 14712000.0)
                    .polynomialSum(t).normalizeDegree()


    // ascending node F
    private fun getAscendingNode(t: Double): Degree =
            doubleArrayOf(93.2720950, 483202.0175233, -0.0036539, -1.0 / 3526000.0, 1.0 / 863310000.0)
                    .polynomialSum(t).normalizeDegree()

    private fun getMeanLunarLongitude(t: Double): Degree =
            doubleArrayOf(218.3164477, 481267.88123421, -0.0015786, 1.0 / 538841.0, -1.0 / 65194000)
                    .polynomialSum(t).normalizeDegree()


    private var args = arrayOf(
            doubleArrayOf(403406.0, 270.54861, 0.9287892),
            doubleArrayOf(195207.0, 340.19128, 35999.1376958),
            doubleArrayOf(119433.0, 63.91854, 35999.4089666),
            doubleArrayOf(112392.0, 331.2622, 35998.7287385),
            doubleArrayOf(3891.0, 317.843, 71998.20261),
            doubleArrayOf(2819.0, 86.631, 71998.4403),
            doubleArrayOf(1721.0, 240.052, 36000.35726),
            doubleArrayOf(660.0, 310.26, 71997.4812),
            doubleArrayOf(350.0, 247.23, 32964.4678),
            doubleArrayOf(334.0, 260.87, -19.441),
            doubleArrayOf(314.0, 297.82, 445267.1117),
            doubleArrayOf(268.0, 343.14, 45036.884),
            doubleArrayOf(242.0, 166.79, 3.1008),
            doubleArrayOf(234.0, 81.53, 22518.4434),
            doubleArrayOf(158.0, 3.5, -19.9739),
            doubleArrayOf(132.0, 132.75, 65928.9345),
            doubleArrayOf(129.0, 182.95, 9038.0293),
            doubleArrayOf(114.0, 162.03, 3034.7684),
            doubleArrayOf(99.0, 29.8, 33718.148),
            doubleArrayOf(93.0, 266.4, 3034.448),
            doubleArrayOf(86.0, 249.2, -2280.773),
            doubleArrayOf(78.0, 157.6, 29929.992),
            doubleArrayOf(72.0, 257.8, 31556.493),
            doubleArrayOf(68.0, 185.1, 149.588),
            doubleArrayOf(64.0, 69.9, 9037.75),
            doubleArrayOf(46.0, 8.0, 107997.405),
            doubleArrayOf(38.0, 197.1, -4444.176),
            doubleArrayOf(37.0, 250.4, 151.771),
            doubleArrayOf(32.0, 65.3, 67555.316),
            doubleArrayOf(29.0, 162.7, 31556.08),
            doubleArrayOf(28.0, 341.5, -4561.54),
            doubleArrayOf(27.0, 291.6, 107996.706),
            doubleArrayOf(27.0, 98.5, 1221.655),
            doubleArrayOf(25.0, 146.7, 62894.167),
            doubleArrayOf(24.0, 110.0, 31437.369),
            doubleArrayOf(21.0, 5.2, 14578.298),
            doubleArrayOf(21.0, 342.6, -31931.757),
            doubleArrayOf(20.0, 230.9, 34777.243),
            doubleArrayOf(18.0, 256.1, 1221.999),
            doubleArrayOf(17.0, 45.3, 62894.511),
            doubleArrayOf(14.0, 242.9, -4442.039),
            doubleArrayOf(13.0, 115.2, 107997.909),
            doubleArrayOf(13.0, 151.8, 119.066),
            doubleArrayOf(13.0, 285.3, 16859.071),
            doubleArrayOf(12.0, 53.3, -4.578),
            doubleArrayOf(10.0, 126.6, 26895.292),
            doubleArrayOf(10.0, 205.7, -39.127),
            doubleArrayOf(10.0, 85.9, 12297.536),
            doubleArrayOf(10.0, 146.1, 90073.778))

    private val dLng = doubleArrayOf(
            0.0, 2.0, 2.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0,
            0.0, 1.0, 0.0, 2.0, 0.0, 0.0, 4.0, 0.0, 4.0, 2.0,
            2.0, 1.0, 1.0, 2.0, 2.0, 4.0, 2.0, 0.0, 2.0, 2.0,
            1.0, 2.0, 0.0, 0.0, 2.0, 2.0, 2.0, 4.0, 0.0, 3.0,
            2.0, 4.0, 0.0, 2.0, 2.0, 2.0, 4.0, 0.0, 4.0, 1.0,
            2.0, 0.0, 1.0, 3.0, 4.0, 2.0, 0.0, 1.0, 2.0, 2.0)

    private val mLng = doubleArrayOf(
            0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, -1.0,
            1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            1.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0,
            0.0, -2.0, 1.0, 2.0, -2.0, 0.0, 0.0, -1.0, 0.0, 0.0,
            1.0, -1.0, 2.0, 2.0, 1.0, -1.0, 0.0, 0.0, -1.0, 0.0,
            1.0, 0.0, 1.0, 0.0, 0.0, -1.0, 2.0, 1.0, 0.0, 0.0)

    private val mpLng = doubleArrayOf(
            1.0, -1.0, 0.0, 2.0, 0.0, 0.0, -2.0, -1.0, 1.0, 0.0,
            -1.0, 0.0, 1.0, 0.0, 1.0, 1.0, -1.0, 3.0, -2.0, -1.0,
            0.0, -1.0, 0.0, 1.0, 2.0, 0.0, -3.0, -2.0, -1.0, -2.0,
            1.0, 0.0, 2.0, 0.0, -1.0, 1.0, 0.0, -1.0, 2.0, -1.0,
            1.0, -2.0, -1.0, -1.0, -2.0, 0.0, 1.0, 4.0, 0.0, -2.0,
            0.0, 2.0, 1.0, -2.0, -3.0, 2.0, 1.0, -1.0, 3.0, -1.0)

    private val fLng = doubleArrayOf(
            0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, -2.0, 2.0, -2.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 2.0, 0.0, 2.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0,
            -2.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0)

    private val sinLng = doubleArrayOf(
            6288774.0, 1274027.0, 658314.0, 213618.0, -185116.0, -114332.0, 58793.0, 57066.0, 53322.0, 45758.0,
            -40923.0, -34720.0, -30383.0, 15327.0, -12528.0, 10980.0, 10675.0, 10034.0, 8548.0, -7888.0,
            -6766.0, -5163.0, 4987.0, 4036.0, 3994.0, 3861.0, 3665.0, -2689.0, -2602.0, 2390.0,
            -2348.0, 2236.0, -2120.0, -2069.0, 2048.0, -1773.0, -1595.0, 1215.0, -1110.0, -892.0,
            -810.0, 759.0, -713.0, -700.0, 691.0, 596.0, 549.0, 537.0, 520.0, -487.0,
            -399.0, -381.0, 351.0, -340.0, 330.0, 327.0, -323.0, 299.0, 294.0, 0.0)

    private val dLat = doubleArrayOf(
            0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0, 0.0, 2.0, 0.0,
            2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0, 4.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 4.0, 4.0,
            0.0, 4.0, 2.0, 2.0, 2.0, 2.0, 0.0, 2.0, 2.0, 2.0,
            2.0, 4.0, 2.0, 2.0, 0.0, 2.0, 1.0, 1.0, 0.0, 2.0,
            1.0, 2.0, 0.0, 4.0, 4.0, 1.0, 4.0, 1.0, 4.0, 2.0)

    private val mLat = doubleArrayOf(
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            -1.0, 0.0, 0.0, 1.0, -1.0, -1.0, -1.0, 1.0, 0.0, 1.0,
            0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            1.0, 0.0, -1.0, -2.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            0.0, -1.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, -2.0)

    private val mpLat = doubleArrayOf(
            0.0, 1.0, 1.0, 0.0, -1.0, -1.0, 0.0, 2.0, 1.0, 2.0,
            0.0, -2.0, 1.0, 0.0, -1.0, 0.0, -1.0, -1.0, -1.0, 0.0,
            0.0, -1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 3.0, 0.0, -1.0,
            1.0, -2.0, 0.0, 2.0, 1.0, -2.0, 3.0, 2.0, -3.0, -1.0,
            0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, -2.0, -1.0,
            1.0, -2.0, 2.0, -2.0, -1.0, 1.0, 1.0, -1.0, 0.0, 0.0)

    private val fLat = doubleArrayOf(
            1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, -1.0,
            -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, -1.0, -1.0, 1.0,
            3.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 1.0,
            -3.0, 1.0, -3.0, -1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0,
            1.0, 1.0, 1.0, -1.0, 3.0, -1.0, -1.0, 1.0, -1.0, -1.0,
            1.0, -1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0)

    private val sinLat = doubleArrayOf(
            5128122.0, 280602.0, 277693.0, 173237.0, 55413.0, 46271.0, 32573.0, 17198.0, 9266.0, 8822.0,
            8216.0, 4324.0, 4200.0, -3359.0, 2463.0, 2211.0, 2065.0, -1870.0, 1828.0, -1794.0,
            -1749.0, -1565.0, -1491.0, -1475.0, -1410.0, -1344.0, -1335.0, 1107.0, 1021.0, 833.0,
            777.0, 671.0, 607.0, 596.0, 491.0, -451.0, 439.0, 422.0, 421.0, -366.0,
            -351.0, 331.0, 315.0, 302.0, -283.0, -229.0, 223.0, 223.0, -220.0, -220.0,
            -185.0, 181.0, -177.0, 176.0, 166.0, -164.0, 132.0, -119.0, 115.0, 107.0)

    private val cosLng = doubleArrayOf(
            -20905355.0, -3699111.0, -2955968.0, -569925.0, 48888.0, -3149.0, 246158.0, -152138.0, -170733.0, -204586.0,
            -129620.0, 108743.0, 104755.0, 10321.0, 0.0, 79661.0, -34782.0, -23210.0, -21636.0, 24208.0,
            30824.0, -8379.0, -16675.0, -12831.0, -10445.0, -11650.0, 14403.0, -7003.0, 0.0, 10056.0,
            6322.0, -9884.0, 5751.0, 0.0, -4950.0, 4130.0, 0.0, -3958.0, 0.0, 3258.0,
            2616.0, -1897.0, -2117.0, 2354.0, 0.0, 0.0, -1423.0, -1117.0, -1571.0, -1739.0,
            0.0, -4421.0, 0.0, 0.0, 0.0, 0.0, 1165.0, 0.0, 0.0, 8752.0)


    @Geocentric @Ecliptic
    val getSunGeocentricEclipticApparentCoordinates: FunGetGeocentricEclipticCoordinates by lazy {
        { t: Double ->
            val aberration = (0.0000974 * cos((177.63 + 35999.01848 * t).toRad()) - 0.005575).normalizeDegree().toRad()
            val longitude = (282.7771834
                    + 36000.76952744 * t
                    + 0.000005729577951308232
                    * args.fold(0.0, { acc, (x, y, z) -> acc + x * sin(y.toRad() + t * z.toRad()) }))
                    .normalizeDegree().toRad()

            val solarAnomaly = (357.5291092 + 35999.0502909 * t - .0001536 * t * t + t * t * t / 24490000.0)
                    .normalizeDegree().toRad()


            var c = (1.9146 - .004817 * t - .000014 * t * t) * Math.sin(solarAnomaly)
            c += (.019993 - .000101 * t) * Math.sin(2 * solarAnomaly)
            c += .00029 * Math.sin(3.0 * solarAnomaly) // Correction to the mean ecliptic longitude
            val nutation: Radian = getNutation(t)

            val ecc = KeplerBodySimonJ2000.Earth.getEccentricity(t)
            val v = solarAnomaly + c * DEG_TO_RAD
            val distance = 1.000001018 * (1.0 - ecc * ecc) / (1.0 + ecc * Math.cos(v)) // In UA

            SphericalVector(longitude + nutation + aberration, 0.0, distance)
        }
    }

    fun getMoonGeocentricEclipticApparentLongitude(t: Double): Radian {

        val solarAnomaly = getSolarAnomaly(t) //M
//        val meanSynodicMonth = 29.530588853
        val lunarElongation = getLunarElongation(t) //D
//
//        val age = meanSynodicMonth * lunarElongation / PI2

        val lunarAnomaly = getLunarAnomaly(t) //ML
        val ascendingNode = getAscendingNode(t) //F

        val E = 1.0 - 0.002516 * t - 0.0000074 * t * t
        val E2 = E * E

        val meanLunarLongitude = getMeanLunarLongitude(t) // L`
        val venus: Degree = 3958.0 / 1000000.0 * sin((119.75 + t * 131.849).toRad())
        val jupiter: Degree = 318.0 / 1000000.0 * sin((53.09 + t * 479264.29).toRad())
        val nutation: Radian = getNutation(t)

        val correctionLongitude: Degree = sinLng.foldIndexed(0.0, {
            index: Int, acc: Double, v: Double ->
            val w = dLng[index]
            val x = mLng[index]
            val y = mpLng[index]
            val z = fLng[index]
            acc + v * (if (abs(x) == 1.0) E else if (abs(x) == 2.0) E2 else 1.0) * sin((w * lunarElongation + x * solarAnomaly + y * lunarAnomaly + z * ascendingNode).toRad())
        }) / 1000000.0

        val flatEarth: Degree = (1962.0 / 1000000.0 * sin((meanLunarLongitude - ascendingNode).toRad()))

        return (meanLunarLongitude + correctionLongitude + venus + jupiter + flatEarth + nutation.toDeg()).toRad()
    }

    fun getMoonGeocentricEclipticApparentLatitude(t: Double): Radian {

        val solarAnomaly = getSolarAnomaly(t) //M
        val lunarElongation = getLunarElongation(t) //D
        val lunarAnomaly = getLunarAnomaly(t) //ML
        val ascendingNode = getAscendingNode(t) //F
        val meanLunarLongitude = getMeanLunarLongitude(t) // L`

        val E = 1.0 - 0.002516 * t - 0.0000074 * t * t
        val E2 = E * E

        val venus: Degree = 175.0 / 1000000.0 * sin((119.75 + t * 131.849).toRad())
        val flatEarth: Degree = (1962.0 / 1000000.0 * sin((meanLunarLongitude - ascendingNode).toRad()))

        val lat0: Degree = sinLat.foldIndexed(0.0, {
            index: Int, acc: Double, v: Double ->
            val w = dLat[index]
            val x = mLat[index]
            val y = mpLat[index]
            val z = fLat[index]
            acc + v * (if (abs(x) == 1.0) E else if (abs(x) == 2.0) E2 else 1.0) * sin((w * lunarElongation + x * solarAnomaly + y * lunarAnomaly + z * ascendingNode).toRad())
        }) / 1000000.0

        val lat = lat0 + venus + flatEarth + 382.0 / 1000000.0 * sin((313.45 + t * 481266.484).toRad())
        return lat.toRad()
    }

    fun getMoonGeocentricApparentDistance(t: Double): Double {

        val solarAnomaly = getSolarAnomaly(t) //M
        val lunarElongation = getLunarElongation(t) //D
        val lunarAnomaly = getLunarAnomaly(t) //ML
        val ascendingNode = getAscendingNode(t) //F

        val E = 1.0 - 0.002516 * t - 0.0000074 * t * t
        val E2 = E * E

        val correction: Double = cosLng.foldIndexed(0.0, {
            index: Int, acc: Double, v: Double ->
            val w = dLng[index]
            val x = mLng[index]
            val y = mpLng[index]
            val z = fLng[index]
            acc + v * (if (abs(x) == 1.0) E else if (abs(x) == 2.0) E2 else 1.0) * cos((w * lunarElongation + x * solarAnomaly + y * lunarAnomaly + z * ascendingNode).toRad())
        })
        return (385000560 + correction) / 1000.0 / AU
    }


    @Geocentric @Ecliptic
    val getMoonGeocentricEclipticApparentCoordinates: FunGetGeocentricEclipticCoordinates by lazy {
        { t: Double ->

            val solarAnomaly = getSolarAnomaly(t) //M
            val lunarElongation = getLunarElongation(t) //D
            val lunarAnomaly = getLunarAnomaly(t) //ML
            val ascendingNode = getAscendingNode(t) //F

            val E = 1.0 - 0.002516 * t - 0.0000074 * t * t
            val E2 = E * E

            val meanLunarLongitude = getMeanLunarLongitude(t) // L`
            val venus: Degree = 3958.0 / 1000000.0 * sin((119.75 + t * 131.849).toRad())
            val jupiter: Degree = 318.0 / 1000000.0 * sin((53.09 + t * 479264.29).toRad())
            val nutation: Radian = getNutation(t)

            val correctionLongitude: Degree = sinLng.foldIndexed(0.0, {
                index: Int, acc: Double, v: Double ->
                val w = dLng[index]
                val x = mLng[index]
                val y = mpLng[index]
                val z = fLng[index]
                acc + v * (if (abs(x) == 1.0) E else if (abs(x) == 2.0) E2 else 1.0) * sin((w * lunarElongation + x * solarAnomaly + y * lunarAnomaly + z * ascendingNode).toRad())
            }) / 1000000.0

            val flatEarth: Degree = (1962.0 / 1000000.0 * sin((meanLunarLongitude - ascendingNode).toRad()))

            val longitude = (meanLunarLongitude + correctionLongitude + venus + jupiter + flatEarth + nutation.toDeg()).toRad()


            val venusLat: Degree = 175.0 / 1000000.0 * sin((119.75 + t * 131.849).toRad())
            val flatEarthLat: Degree = (1962.0 / 1000000.0 * sin((meanLunarLongitude - ascendingNode).toRad()))

            val lat0: Degree = sinLat.foldIndexed(0.0, {
                index: Int, acc: Double, v: Double ->
                val w = dLat[index]
                val x = mLat[index]
                val y = mpLat[index]
                val z = fLat[index]
                acc + v * (if (abs(x) == 1.0) E else if (abs(x) == 2.0) E2 else 1.0) * sin((w * lunarElongation + x * solarAnomaly + y * lunarAnomaly + z * ascendingNode).toRad())
            }) / 1000000.0

            val lat = lat0 + venusLat + flatEarthLat + 382.0 / 1000000.0 * sin((313.45 + t * 481266.484).toRad())
            val latitude = lat.toRad()

            val correction: Double = cosLng.foldIndexed(0.0, {
                index: Int, acc: Double, v: Double ->
                val w = dLng[index]
                val x = mLng[index]
                val y = mpLng[index]
                val z = fLng[index]
                acc + v * (if (abs(x) == 1.0) E else if (abs(x) == 2.0) E2 else 1.0) * cos((w * lunarElongation + x * solarAnomaly + y * lunarAnomaly + z * ascendingNode).toRad())
            })
            val distance = (385000560 + correction) / 1000.0 / AU


            SphericalVector(longitude, latitude, distance)
        }
    }


}