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

package net.arwix.astronomy.swiss

import net.arwix.astronomy.core.ARCSEC_TO_RAD
import net.arwix.astronomy.core.JULIAN_DAYS_PER_CENTURY
import net.arwix.astronomy.core.kepler.EarthMoonElements
import net.arwix.astronomy.core.kepler.KeplerBodySimonJ2000
import net.arwix.astronomy.core.vector.RectangularVector
import net.arwix.astronomy.core.vector.SphericalVector
import net.arwix.astronomy.core.vector.Vector
import net.arwix.astronomy.math.mod3600
import net.arwix.astronomy.math.radians.normalize


/**
 * Prepare lookup table of sin and cos ( i*Lj ) for required multiple angles.
 * @param k
 * @param arg
 * @param n
 */
private fun sscc(k: Int, arg: Double, n: Int, ssArray: Array<DoubleArray>, ccArray: Array<DoubleArray>) {
    var cv: Double
    var sv: Double
    var oldCv: Double
    val su = Math.sin(arg)
    val cu = Math.cos(arg)
    ssArray[k][0] = su /* sin(L) */
    ccArray[k][0] = cu /* cos(L) */
    sv = 2.0 * su * cu
    cv = cu * cu - su * su
    ssArray[k][1] = sv /* sin(2L) */
    ccArray[k][1] = cv
    var i = 2
    while (i < n) {
        oldCv = cv
        cv = cu * cv - su * sv
        sv = su * oldCv + cu * sv
        ssArray[k][i] = sv /* sin( i+1 L ) */
        ccArray[k][i] = cv
        i++
    }
}

/**
 * Generic program to accumulate sum of trigonometric series in three
 * variables (e.g., longitude, latitude, radius) of the same list of
 * arguments.

 * @param tt Julian centuries
 * @param arg_tbl
 * @param distance
 * @param lat_tbl
 * @param lon_tbl
 * @param rad_tbl
 * @param max_harmonic
 * @param max_power_of_t
 * @param maxargs
 * @param timescale
 * @param trunclvl
 * @return An array with x, y, z (AU).
 */
internal fun gplan(tt: Double, arg_tbl: IntArray, distance: Double, lat_tbl: DoubleArray, lon_tbl: DoubleArray,
                   rad_tbl: DoubleArray, max_harmonic: IntArray, max_power_of_t: Int, maxargs: Int, timescale: Double, trunclvl: Double): Vector {


    /* Compute mean elements at Julian date J. */
    val ss = Array(20) { DoubleArray(41) }
    val cc = Array(20) { DoubleArray(41) }
    val T = tt * JULIAN_DAYS_PER_CENTURY / timescale

    /* From Simon et al (1994) */
    val freqs = doubleArrayOf(
            /* Arc sec per 10000 Julian years. */
            53810162868.8982, 21066413643.3548, 12959774228.3429, 6890507749.3988, 1092566037.7991, 439960985.5372, 154248119.3933, 78655032.0744, 52272245.1795)
    val phases = doubleArrayOf(
            /* Arc sec. */
            252.25090552 * 3600.0, 181.97980085 * 3600.0, 100.46645683 * 3600.0, 355.43299958 * 3600.0, 34.35151874 * 3600.0, 50.07744430 * 3600.0, 314.05500511 * 3600.0, 304.34866548 * 3600.0, 860492.1546)

    max_harmonic.forEachIndexed { index, harmonic -> if (harmonic > 0) sscc(index, ((freqs[index] * T).mod3600() + phases[index]) * ARCSEC_TO_RAD, harmonic, ss, cc) }

    val vector = accumulate(T, arg_tbl, lon_tbl, lat_tbl, rad_tbl, ss, cc, true)

    //TODO isRectangular?
    if (distance == 0.0) return RectangularVector((ARCSEC_TO_RAD * vector.phi).normalize(), (ARCSEC_TO_RAD * vector.theta).normalize(), (ARCSEC_TO_RAD * vector.r).normalize())

    vector.set(vector.phi * ARCSEC_TO_RAD, vector.theta * ARCSEC_TO_RAD, distance * (1.0 + ARCSEC_TO_RAD * vector.r))
    return vector //.getVectorOfType(VectorType.RECTANGULAR).toArray() // doubleArrayOf(x, y, z)
}


/**
 * Generic program to accumulate sum of trigonometric series in three
 * variables (e.g., longitude, latitude, radius) of the same list of
 * arguments.

 * @param tt Julian centuries
 * @param arg_tbl
 * @param distance
 * @param lat_tbl
 * @param lon_tbl
 * @param rad_tbl
 * @param max_harmonic
 * @param max_power_of_t
 * @param maxargs
 * @param timescale
 * @param trunclvl
 * @return An array with x, y, z (AU).
 */
internal fun g3plan(tt: Double, arg_tbl: IntArray, distance: Double, lat_tbl: DoubleArray, lon_tbl: DoubleArray,
                    rad_tbl: DoubleArray, max_harmonic: IntArray, max_power_of_t: Int, maxargs: Int, timescale: Double, trunclvl: Double,
                    libration: Boolean): Vector {

    val ss = Array(20) { DoubleArray(41) }
    val cc = Array(20) { DoubleArray(41) }

    val T100 = tt * JULIAN_DAYS_PER_CENTURY / timescale

    val args = meanElements(tt)
    val pA_precession = EarthMoonElements(tt).precession
    val Ea_arcsec = args[2]
    if (libration) args[13] -= pA_precession // Only librations
    /* Calculate sin( i*MM ), etc. for needed multiple angles. */

    max_harmonic.forEachIndexed { i, harmonic -> if (harmonic > 0) sscc(i, args[i], max_harmonic[i], ss, cc) }

    val vector = accumulate(T100, arg_tbl, lon_tbl, lat_tbl, rad_tbl, ss, cc)
    vector.set(vector.phi * 0.0001, vector.theta * 0.0001, vector.r * 0.0001)
    //TODO isRectangular?
    if (distance == 0.0) return RectangularVector(ARCSEC_TO_RAD * vector.phi + Ea_arcsec, ARCSEC_TO_RAD * vector.theta, ARCSEC_TO_RAD * vector.r)

    vector.set(vector.phi * ARCSEC_TO_RAD + Ea_arcsec, vector.theta * ARCSEC_TO_RAD, distance * (1.0 + ARCSEC_TO_RAD * vector.r))
    return vector; //.getVectorOfType(VectorType.RECTANGULAR) as RectangularVector // doubleArrayOf(x, y, z)
}


/**
 * Generic program to accumulate sum of trigonometric series in one
 * variables (e.g., latitude) of the same list of arguments.

 * @param tt Julian centuries
 * @param arg_tbl
 * @param lat_tbl
 * @param max_harmonic
 * @param timescale
 * @return Latitude (rad).
 */
internal fun g1plan(tt: Double, arg_tbl: IntArray, lat_tbl: DoubleArray, max_harmonic: IntArray, timescale: Double): Double {

    var m: Int
    var k: Int
    var j: Int
    var buffer: Double

    var su: Double
    var cu: Double

    var np: Int

    val ss = Array(20) { DoubleArray(41) }
    val cc = Array(20) { DoubleArray(41) }

    val T100 = tt * JULIAN_DAYS_PER_CENTURY / timescale

    val args = meanElements(tt)

    /* Calculate sin( i*MM ), etc. for needed multiple angles. */
    max_harmonic.forEachIndexed { i, harmonic -> if (harmonic > 0) sscc(i, args[i], max_harmonic[i], ss, cc) }

    var sb = 0.0

    var p_index = -1
    var pb_index = -1

    while (true) {
        /* argument of sine and cosine */
        /* Number of periodic arguments. */
        np = arg_tbl[++p_index]

        if (np < 0) break
        if (np == 0) { /* It is a polynomial term. */
            sb += (0..arg_tbl[++p_index] - 1).fold(lat_tbl[++pb_index], { acc, _ -> acc * T100 + lat_tbl[++pb_index] })
            continue
        }

        val (sv, cv) = (0..np - 1).fold(doubleArrayOf(0.0, 1.0), { acc, i ->
            // What harmonic
            j = arg_tbl[++p_index]
            if (j == 0) return@fold if (i == np - 1 && acc[0] == 0.0 && acc[1] == 1.0) doubleArrayOf(0.0, 0.0) else acc
            // Which planet
            m = arg_tbl[++p_index] - 1
            k = Math.abs(j) - 1
            // sin(k*angle)
            su = if (j < 0) -ss[m][k] else ss[m][k]
            cu = cc[m][k]
            buffer = su * acc[1] + cu * acc[0]
            acc[1] = cu * acc[1] - su * acc[0]
            acc[0] = buffer
            acc
        })

        /* Highest power of T. */
        /* Latitude. */
        val range = (0..arg_tbl[++p_index] - 1)

        sb += range
                .fold(doubleArrayOf(lat_tbl[++pb_index], lat_tbl[++pb_index]), { acc, _ ->
                    acc[0] = acc[0] * T100 + lat_tbl[++pb_index]
                    acc[1] = acc[1] * T100 + lat_tbl[++pb_index]
                    acc
                })
                .let { it[0] * cv + it[1] * sv }
    }

    return ARCSEC_TO_RAD * sb * 0.0001
}

/**
 * Generic program to accumulate sum of trigonometric series in two
 * variables (e.g., longitude, radius) of the same list of arguments.

 * @param tt Julian centuries
 * @param arg_tbl
 * @param distance
 * @param lon_tbl
 * @param rad_tbl
 * @param max_harmonic
 * @param timescale
 * @return An array with x, y, z (AU).
 */
internal fun g2plan(tt: Double, arg_tbl: IntArray, distance: Double, lon_tbl: DoubleArray, rad_tbl: DoubleArray,
                    max_harmonic: IntArray, timescale: Double, lat: Double): Vector {

    var np: Int

    var m: Int
    var k: Int
    var j: Int
    var buffer: Double

    var su: Double
    var cu: Double

    val ss = Array(20) { DoubleArray(41) }
    val cc = Array(20) { DoubleArray(41) }

    val T100 = tt * JULIAN_DAYS_PER_CENTURY / timescale
    val args = meanElements(tt)
    val LP_equinox = args[13]

    /* Calculate sin( i*MM ), etc. for needed multiple angles. */

    max_harmonic.forEachIndexed { i, harmonic -> if (harmonic > 0) sscc(i, args[i], max_harmonic[i], ss, cc) }

    val vector = SphericalVector(0.0, 0.0, 0.0)

    var p_index = -1
    var pl_index = -1
    var pr_index = -1

    while (true) {
        // argument of sine and cosine
        // Number of periodic arguments
        np = arg_tbl[++p_index]
        if (np < 0) break
        if (np == 0) {// It is a polynomial term
            val range = (0..arg_tbl[++p_index] - 1)
            // "Longitude" polynomial (phi)
            vector.phi += range.fold(lon_tbl[++pl_index], { acc, _ -> acc * T100 + lon_tbl[++pl_index] })
            // Radius polynomial (psi)
            vector.r += range.fold(rad_tbl[++pr_index], { acc, _ -> acc * T100 + rad_tbl[++pr_index] })
            continue
        }

        val (sv, cv) = (0..np - 1).fold(doubleArrayOf(0.0, 1.0), { acc, i ->
            // What harmonic
            j = arg_tbl[++p_index]
            if (j == 0) return@fold if (i == np - 1 && acc[0] == 0.0 && acc[1] == 1.0) doubleArrayOf(0.0, 0.0) else acc
            // Which planet
            m = arg_tbl[++p_index] - 1
            k = Math.abs(j) - 1
            // sin(k*angle)
            su = if (j < 0) -ss[m][k] else ss[m][k]
            cu = cc[m][k]
            buffer = su * acc[1] + cu * acc[0]
            acc[1] = cu * acc[1] - su * acc[0]
            acc[0] = buffer
            acc
        })

        /* Highest power of T. */
        /* Longitude. */
        val range = (0..arg_tbl[++p_index] - 1)

        vector.phi += range
                .fold(doubleArrayOf(lon_tbl[++pl_index], lon_tbl[++pl_index]), { acc, _ ->
                    acc[0] = acc[0] * T100 + lon_tbl[++pl_index]
                    acc[1] = acc[1] * T100 + lon_tbl[++pl_index]
                    acc
                })
                .let { it[0] * cv + it[1] * sv }

        /* Radius. */
        vector.r += range
                .fold(doubleArrayOf(rad_tbl[++pr_index], rad_tbl[++pr_index]), { acc, _ ->
                    acc[0] = acc[0] * T100 + rad_tbl[++pr_index]
                    acc[1] = acc[1] * T100 + rad_tbl[++pr_index]
                    acc
                })
                .let { it[0] * cv + it[1] * sv }
    }

    vector.phi *= 0.0001
    vector.r *= 0.0001

    //TODO Spherical?
    if (distance == 0.0)
        return RectangularVector(ARCSEC_TO_RAD * vector.phi + LP_equinox, lat, ARCSEC_TO_RAD * vector.r)

    vector.phi = ARCSEC_TO_RAD * vector.phi + LP_equinox
    vector.theta = lat
    vector.r = distance * (1.0 + ARCSEC_TO_RAD * vector.r)

    return vector.toType<RectangularVector>()
}


/**
 * Obtain mean elements of the planets.
 * @param T Julian centuries
 * @return An array with the mean longitudes.
 */
private fun meanElements(T: Double): DoubleArray {


    val Args = DoubleArray(20)

    /** Mean longitudes of planets (Simon et al, 1994) .047" subtracted from
     * constant term for offset to DE403 origin.
     */

    val delta = -.047 * ARCSEC_TO_RAD

    Args[0] = KeplerBodySimonJ2000.Mercury(T).Longitude + delta
    Args[1] = KeplerBodySimonJ2000.Venus(T).Longitude + delta
    Args[2] = KeplerBodySimonJ2000.Earth(T).Longitude + delta
    Args[3] = KeplerBodySimonJ2000.Mars(T).Longitude + delta
    Args[4] = KeplerBodySimonJ2000.Jupiter(T).Longitude + delta
    Args[5] = KeplerBodySimonJ2000.Saturn(T).Longitude + delta
    Args[6] = KeplerBodySimonJ2000.Uranus(T).Longitude + delta
    Args[7] = KeplerBodySimonJ2000.Neptune(T).Longitude + delta

    /* Copied from cmoon.c, DE404 version. */
    /* Mean elongation of moon = elongation */
    val elements = EarthMoonElements(T)
    Args[9] = elements.elongation
    Args[10] = elements.ascendingNode
    Args[11] = elements.sunAnomaly
    Args[12] = elements.moonAnomaly
    Args[13] = elements.longitude
    Args[14] = elements.LA
    Args[15] = elements.LB
    Args[16] = elements.LC
    Args[17] = Args[13] - Args[10]
    Args[18] = elements.NB
    return Args
}

private fun accumulate(T: Double, arg_tbl: IntArray, lon_tbl: DoubleArray, lat_tbl: DoubleArray, rad_tbl: DoubleArray,
                       ss: Array<DoubleArray>, cc: Array<DoubleArray>, isMod: Boolean = false): SphericalVector {
    var p_index = -1
    var pl_index = -1
    var pb_index = -1
    var pr_index = -1

    var m: Int
    var k: Int
    var j: Int
    var buffer: Double

    var su: Double
    var cu: Double

    var np: Int
    val vector = SphericalVector(0.0, 0.0, 0.0)
    while (true) {
        // argument of sine and cosine
        // Number of periodic arguments
        np = arg_tbl[++p_index]
        if (np < 0) break
        if (np == 0) {// It is a polynomial term
            val range = (0..arg_tbl[++p_index] - 1)
            // "Longitude" polynomial (phi)
            vector.phi += range.fold(lon_tbl[++pl_index], { acc, _ -> acc * T + lon_tbl[++pl_index] }).let { if (isMod) it.mod3600() else it }

            // "Latitude" polynomial (theta)
            vector.theta += range.fold(lat_tbl[++pb_index], { acc, _ -> acc * T + lat_tbl[++pb_index] })

            // Radius polynomial (psi)
            vector.r += range.fold(rad_tbl[++pr_index], { acc, _ -> acc * T + rad_tbl[++pr_index] })
            continue
        }


        val (sv, cv) = (0..np - 1).fold(doubleArrayOf(0.0, 1.0), { acc, i ->
            // What harmonic
            j = arg_tbl[++p_index]
            if (j == 0) return@fold if (i == np - 1 && acc[0] == 0.0 && acc[1] == 1.0) doubleArrayOf(0.0, 0.0) else acc
            // Which planet
            m = arg_tbl[++p_index] - 1
            k = Math.abs(j) - 1
            // sin(k*angle)
            su = if (j < 0) -ss[m][k] else ss[m][k]
            cu = cc[m][k]
            buffer = su * acc[1] + cu * acc[0]
            acc[1] = cu * acc[1] - su * acc[0]
            acc[0] = buffer
            acc
        })

        /* Highest power of T. */
        /* Longitude. */
        val range = (0..arg_tbl[++p_index] - 1)

        vector.phi += range
                .fold(doubleArrayOf(lon_tbl[++pl_index], lon_tbl[++pl_index]), { acc, _ ->
                    acc[0] = acc[0] * T + lon_tbl[++pl_index]
                    acc[1] = acc[1] * T + lon_tbl[++pl_index]
                    acc
                })
                .let { it[0] * cv + it[1] * sv }

        /* Latitude. */
        vector.theta += range
                .fold(doubleArrayOf(lat_tbl[++pb_index], lat_tbl[++pb_index]), { acc, _ ->
                    acc[0] = acc[0] * T + lat_tbl[++pb_index]
                    acc[1] = acc[1] * T + lat_tbl[++pb_index]
                    acc
                })
                .let { it[0] * cv + it[1] * sv }

        /* Radius. */
        vector.r += range
                .fold(doubleArrayOf(rad_tbl[++pr_index], rad_tbl[++pr_index]), { acc, _ ->
                    acc[0] = acc[0] * T + rad_tbl[++pr_index]
                    acc[1] = acc[1] * T + rad_tbl[++pr_index]
                    acc
                })
                .let { it[0] * cv + it[1] * sv }
    }

    return vector
}
