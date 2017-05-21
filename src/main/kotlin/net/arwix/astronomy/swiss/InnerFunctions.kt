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
import net.arwix.astronomy.core.vector.SphericalVector
import net.arwix.astronomy.core.vector.VectorType
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

 * @param J Julian day.
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
                   rad_tbl: DoubleArray, max_harmonic: IntArray, max_power_of_t: Int, maxargs: Int, timescale: Double, trunclvl: Double): DoubleArray {


    /* Compute mean elements at Julian date J. */
    val ss = Array(20) { DoubleArray(41) }
    val cc = Array(20) { DoubleArray(41) }

    var j: Int
    var k: Int
    var m: Int
    var su: Double
    var cu: Double
    var buffer: Double

    val T = tt * JULIAN_DAYS_PER_CENTURY / timescale

    /* From Simon et al (1994) */
    val freqs = doubleArrayOf(
            /* Arc sec per 10000 Julian years. */
            53810162868.8982, 21066413643.3548, 12959774228.3429, 6890507749.3988, 1092566037.7991, 439960985.5372, 154248119.3933, 78655032.0744, 52272245.1795)
    val phases = doubleArrayOf(
            /* Arc sec. */
            252.25090552 * 3600.0, 181.97980085 * 3600.0, 100.46645683 * 3600.0, 355.43299958 * 3600.0, 34.35151874 * 3600.0, 50.07744430 * 3600.0, 314.05500511 * 3600.0, 304.34866548 * 3600.0, 860492.1546)

    max_harmonic.forEachIndexed { index, harmonic ->
        if (harmonic > 0) sscc(index, ((freqs[index] * T).mod3600() + phases[index]) * ARCSEC_TO_RAD, harmonic, ss, cc)
    }

    var np: Int
    val vector = SphericalVector(0.0, 0.0, 0.0)
    var p_index = -1
    var pl_index = -1
    var pb_index = -1
    var pr_index = -1

    while (true) {
        // argument of sine and cosine
        // Number of periodic arguments
        np = arg_tbl[++p_index]
        if (np < 0) break
        if (np == 0) {// It is a polynomial term
            val range = (0..arg_tbl[++p_index] - 1)
            // "Longitude" polynomial (phi)
            vector.phi += range.fold(lon_tbl[++pl_index], { acc, _ -> acc * T + lon_tbl[++pl_index] }).mod3600()

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

    if (distance == 0.0)
        return doubleArrayOf((ARCSEC_TO_RAD * vector.phi).normalize(), (ARCSEC_TO_RAD * vector.theta).normalize(), (ARCSEC_TO_RAD * vector.r).normalize())


    vector.set(vector.phi * ARCSEC_TO_RAD, vector.theta * ARCSEC_TO_RAD, distance * (1.0 + ARCSEC_TO_RAD * vector.r))


    return vector.getVectorOfType(VectorType.RECTANGULAR).toArray() // doubleArrayOf(x, y, z)
}


/**
 * Generic program to accumulate sum of trigonometric series in three
 * variables (e.g., longitude, latitude, radius) of the same list of
 * arguments.

 * @param J Julian day.
 *
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
//private fun g3plan(tt: Double, arg_tbl: IntArray, distance: Double, lat_tbl: DoubleArray, lon_tbl: DoubleArray,
//                   rad_tbl: DoubleArray, max_harmonic: IntArray, max_power_of_t: Int, maxargs: Int, timescale: Double, trunclvl: Double,
//                   libration: Boolean): DoubleArray {
//
//    /* Compute mean elements at Julian date J. */
//    val ss = Array(20) { DoubleArray(41) }
//    val cc = Array(20) { DoubleArray(41) }
//
//    var i: Int
//    var j: Int
//    var k: Int
//    var m: Int
//    var k1: Int
//    var ip: Int
//    var np: Int
//    var nt: Int
//    val p: IntArray
//    val pl: DoubleArray
//    val pb: DoubleArray
//    val pr: DoubleArray
//    var su: Double
//    var cu: Double
//    var sv: Double
//    var cv: Double
//    var t: Double
//    var sl: Double
//    var sb: Double
//    var sr: Double
//
//    val args = PlanetEphem.meanElements(J)
//    if (libration) args[13] -= PlanetEphem.pA_precession // Only librations
//    val T = tt * JULIAN_DAYS_PER_CENTURY / timescale
//
//    /* Calculate sin( i*MM ), etc. for needed multiple angles. */
//    i = 0
//    while (i < maxargs) {
//        if (max_harmonic[i] > 0) {
//            sscc(i, args[i], max_harmonic[i])
//        }
//        i++
//    }
//
//    /* Point to start of table of arguments. */
//    p = arg_tbl
//
//    /* Point to tabulated cosine and sine amplitudes. */
//    pl = lon_tbl
//    pb = lat_tbl
//    pr = rad_tbl
//
//    sl = 0.0
//    sb = 0.0
//    sr = 0.0
//
//    np = 0
//    nt = 0
//    cu = 0.0
//
//    var p_index = -1
//    var pl_index = -1
//    var pb_index = -1
//    var pr_index = -1
//
//    while (true) {
//        /* argument of sine and cosine */
//        /* Number of periodic arguments. */
//        p_index++
//        np = p[p_index]
//        if (np < 0)
//            break
//        if (np == 0) { /* It is a polynomial term. */
//            p_index++
//            nt = p[p_index]
//            /* "Longitude" polynomial (phi). */
//            pl_index++
//            cu = pl[pl_index]
//            ip = 0
//            while (ip < nt) {
//                pl_index++
//                cu = cu * T + pl[pl_index]
//                ip++
//            }
//            sl += cu
//            /* "Latitude" polynomial (theta). */
//            pb_index++
//            cu = pb[pb_index]
//            ip = 0
//            while (ip < nt) {
//                pb_index++
//                cu = cu * T + pb[pb_index]
//                ip++
//            }
//            sb += cu
//            /* Radius polynomial (psi). */
//            pr_index++
//            cu = pr[pr_index]
//            ip = 0
//            while (ip < nt) {
//                pr_index++
//                cu = cu * T + pr[pr_index]
//                ip++
//            }
//            sr += cu
//            continue
//        }
//
//        k1 = 0
//        cv = 0.0
//        sv = 0.0
//        ip = 0
//        while (ip < np) {
//            /* What harmonic. */
//            p_index++
//            j = p[p_index]
//            /* Which planet. */
//            p_index++
//            m = p[p_index] - 1
//            if (j != 0) {
//                k = Math.abs(j) - 1
//
//                su = ss[m][k] /* sin(k*angle) */
//                if (j < 0)
//                    su = -su
//
//                cu = cc[m][k]
//                if (k1 == 0) { /* set first angle */
//                    sv = su
//                    cv = cu
//                    k1 = 1
//                } else { /* combine angles */
//                    t = su * cv + cu * sv
//                    cv = cu * cv - su * sv
//                    sv = t
//                }
//            }
//            ip++
//        }
//
//        /* Highest power of T. */
//        p_index++
//        nt = p[p_index]
//        /* Longitude. */
//        pl_index++
//        cu = pl[pl_index]
//        pl_index++
//        su = pl[pl_index]
//        ip = 0
//        while (ip < nt) {
//            pl_index++
//            cu = cu * T + pl[pl_index]
//            pl_index++
//            su = su * T + pl[pl_index]
//            ip++
//        }
//        sl += cu * cv + su * sv
//        /* Latitude. */
//        pb_index++
//        cu = pb[pb_index]
//        pb_index++
//        su = pb[pb_index]
//        ip = 0
//        while (ip < nt) {
//            pb_index++
//            cu = cu * T + pb[pb_index]
//            pb_index++
//            su = su * T + pb[pb_index]
//            ip++
//        }
//        sb += cu * cv + su * sv
//        /* Radius. */
//        pr_index++
//        cu = pr[pr_index]
//        pr_index++
//        su = pr[pr_index]
//        ip = 0
//        while (ip < nt) {
//            pr_index++
//            cu = cu * T + pr[pr_index]
//            pr_index++
//            su = su * T + pr[pr_index]
//            ip++
//        }
//        sr += cu * cv + su * sv
//    }
//
//    sl = sl * 0.0001
//    sb = sb * 0.0001
//    sr = sr * 0.0001
//
//    if (distance == 0.0)
//        return doubleArrayOf(ARCSEC_TO_RAD * sl + PlanetEphem.Ea_arcsec, ARCSEC_TO_RAD * sb, ARCSEC_TO_RAD * sr)
//
//    val pobj = DoubleArray(3)
//    pobj[0] = ARCSEC_TO_RAD * sl + PlanetEphem.Ea_arcsec
//    pobj[1] = ARCSEC_TO_RAD * sb
//    pobj[2] = distance * (1.0 + ARCSEC_TO_RAD * sr)
//
//    val x = pobj[2] * Math.cos(pobj[0]) * Math.cos(pobj[1])
//    val y = pobj[2] * Math.sin(pobj[0]) * Math.cos(pobj[1])
//    val z = pobj[2] * Math.sin(pobj[1])
//
//    return doubleArrayOf(x, y, z)
//}
