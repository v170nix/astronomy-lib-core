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

import net.arwix.astronomy.core.ARCSEC_TO_RAD
import net.arwix.astronomy.math.mod3600


/**
 * IAU1980 Theory of Nutation, with optional support for free core nutation.
 * Results are set in fields [Nutation.nutationInLongitude] and
 * [Nutation.nutationInObliquity].
 * <P>
 * Each term in the expansion has a trigonometric argument given by
 * W = i*MM + j*MS + k*FF + l*DD + m*OM
 * where the variables are defined below.
 * The nutation in longitude is a sum of terms of the form (a + bT) *
 * sin(W). The terms for nutation in obliquity are of the form (c + dT) *
 * cos(W). The coefficients are arranged in the tabulation as follows:
 * Coefficient:

</P> * <pre>
 * i  j  k  l  m      a      b      c     d
 * <BR></BR>
 * 0, 0, 0, 0, 1, -171996, -1742, 92025, 89
</pre> *

 * The first line of the table, above, is done separately since two of the
 * values do not fit into 16 bit integers. The values a and c are arc
 * seconds times 10000. b and d are arc seconds per Julian century times
 * 100000. i through m are integers. See the program for interpretation of
 * MM, MS, etc., which are mean orbital elements of the Sun and Moon.
 *
 * If terms with coefficient less than X are omitted, the peak errors will
 * be:
 * <pre>
 * omit a &lt;    error in longitude,     omit c &lt;    error in obliquity
 * <BR></BR>
 * .0005&quot;      .0100&quot;                  .0008&quot;      .0094&quot;
 * <BR></BR>
 * .0046       .0492                   .0095       .0481
 * <BR></BR>
 * .0123       .0880                   .0224       .0905
 * <BR></BR>
 * .0386       .1808                   .0895       .1129
</pre></P> * <P>
 * References:
 * "Summary of 1980 IAU Theory of Nutation (Final Report of the IAU Working
 * Group on Nutation)", P. K. Seidelmann et al., in Transactions of the IAU
 * Vol. XVIII A, Reports on Astronomy, P. A. Wayman, ed.; elongation. Reidel Pub.
 * Co., 1982.
 * "Nutation and the Earth's Rotation", I.A.U. Symposium No. 78, May, 1977,
 * page 256. I.A.U., 1980.
 * Woolard, E.W., "A redevelopment of the theory of nutation", The
 * Astronomical Journal, 58, 1-3 (1953).
 * This program implements all of the 1980 IAU nutation series. Results
 * checked at 100 points against the 1986 Astronomical Almanac, all agreed.
 * Translated to Java from code by S. L. Moshier, November 1987
 * October, 1992 - typo fixed in nutation matrix
 * October, 1995 - fixed typo in node argument, tested against JPL DE403
 * ephemeris file.
 * @param T Time in Julian centuries from J2000.
 * */
internal fun calcNutation_IAU1980(T: Double): NutationResult {

    /**
     * Array to hold sines of multiple angles
     */
    val ss = Array(5) { DoubleArray(8) }

    /**
     * Array to hold cosines of multiple angles
     */
    val cc = Array(5) { DoubleArray(8) }


    /*
	 * Prepare lookup table of sin and cos ( i*Lj ) for required multiple angles
	 */
    fun sscc(k: Int, arg: Double, n: Int) {
        var s: Double
        val su = Math.sin(arg)
        val cu = Math.cos(arg)
        ss[k][0] = su /* sin(L) */
        cc[k][0] = cu /* cos(L) */
        var sv = 2.0 * su * cu
        var cv = cu * cu - su * su
        ss[k][1] = sv /* sin(2L) */
        cc[k][1] = cv

        var i = 2
        while (i < n) {
            s = su * cv + cu * sv
            cv = cu * cv - su * sv
            sv = s
            ss[k][i] = sv /* sin( i+1 L ) */
            cc[k][i] = cv
            i++
        }
    }


    var f: Double
    var g: Double
    var cu: Double
    var su: Double
    var cv: Double
    var sv: Double
    var sw: Double

    var j: Int
    var k: Int
    var k1: Int
    var m: Int

    val T2 = T * T
    val T10 = T / 10.0

    /* Fundamental arguments in the FK5 reference system. */

    /*
		 * longitude of the mean ascending node of the lunar orbit on the
		 * ecliptic, measured from the mean equinox of date
		 */
    val OM = ((-6962890.539 * T + 450160.280).mod3600() + (0.008 * T + 7.455) * T2) * ARCSEC_TO_RAD

    /*
		 * mean longitude of the Sun minus the mean longitude of the Sun's
		 * perigee
		 */
    val MS = ((129596581.224 * T + 1287099.804).mod3600() - (0.012 * T + 0.577) * T2) * ARCSEC_TO_RAD

    /*
		 * mean longitude of the Moon minus the mean longitude of the Moon's
		 * perigee
		 */
    val MM = ((1717915922.633 * T + 485866.733).mod3600() + (0.064 * T + 31.310) * T2) * ARCSEC_TO_RAD

    /*
		 * mean longitude of the Moon minus the mean longitude of the Moon's
		 * node
		 */
    val FF = ((1739527263.137 * T + 335778.877).mod3600() + (0.011 * T - 13.257) * T2) * ARCSEC_TO_RAD

    /*
		 * mean elongation of the Moon from the Sun.
		 */
    val DD = ((1602961601.328 * T + 1072261.307).mod3600() + (0.019 * T - 6.891) * T2) * ARCSEC_TO_RAD

    /*
		 * Calculate sin( i*MM ), etc. for needed multiple angles
		 */
    sscc(0, MM, 3)
    sscc(1, MS, 2)
    sscc(2, FF, 4)
    sscc(3, DD, 4)
    sscc(4, OM, 2)

    var C = 0.0
    var D = 0.0
    var p = -1 /* point to start of table */

    var i = 0
    while (i < 105) {
        /* argument of sine and cosine */
        k1 = 0
        cv = 0.0
        sv = 0.0
        m = 0
        while (m < 5) {
            p++
            j = net.arwix.astronomy.ephemeris.IAU1980_NT[p]
            if (j != 0) {
                k = j
                if (j < 0)
                    k = -k
                su = ss[m][k - 1] /* sin(k*angle) */
                if (j < 0)
                    su = -su
                cu = cc[m][k - 1]
                if (k1 == 0) { /* set first angle */
                    sv = su
                    cv = cu
                    k1 = 1
                } else { /* combine angles */
                    sw = su * cv + cu * sv
                    cv = cu * cv - su * sv
                    sv = sw
                }
            }
            m++
        }
        /* longitude coefficient */
        p++
        f = IAU1980_NT[p].toDouble()
        p++
        k = IAU1980_NT[p]
        if (k != 0)
            f += T10 * k

        /* obliquity coefficient */
        p++
        g = IAU1980_NT[p].toDouble()
        p++
        k = IAU1980_NT[p]
        if (k != 0)
            g += T10 * k

        /* accumulate the terms */
        C += f * sv
        D += g * cv
        i++
    }

    /* first terms, not in table: */
    C += (-1742.0 * T10 - 171996.0) * ss[4][0] /* sin(OM) */
    D += (89.0 * T10 + 92025.0) * cc[4][0] /* cos(OM) */

    /* Save answers, expressed in radians */
    val nutationInLongitude = 0.0001 * ARCSEC_TO_RAD * C
    val nutationInObliquity = 0.0001 * ARCSEC_TO_RAD * D


    return NutationResult(nutationInLongitude, nutationInObliquity)
}


/**
 * IAU1980 model
 */
private val IAU1980_NT by lazy {
    intArrayOf(0, 0, 0, 0, 2, 2062, 2, -895, 5, -2, 0, 2, 0, 1, 46, 0, -24, 0, 2, 0, -2, 0, 0, 11, 0, 0, 0, -2, 0, 2, 0, 2, -3, 0, 1, 0, 1, -1, 0, -1, 0, -3, 0, 0, 0, 0, -2, 2, -2, 1, -2, 0, 1, 0, 2, 0, -2, 0, 1, 1, 0, 0, 0, 0, 0, 2, -2, 2, -13187, -16, 5736, -31, 0, 1, 0, 0, 0, 1426, -34, 54, -1, 0, 1, 2, -2, 2, -517, 12, 224, -6, 0, -1, 2, -2, 2, 217, -5, -95, 3, 0, 0, 2, -2, 1, 129, 1, -70, 0, 2, 0, 0, -2, 0, 48, 0, 1, 0, 0, 0, 2, -2, 0, -22, 0, 0, 0, 0, 2, 0, 0, 0, 17, -1, 0, 0, 0, 1, 0, 0, 1, -15, 0, 9, 0, 0, 2, 2, -2, 2, -16, 1, 7, 0, 0, -1, 0, 0, 1, -12, 0, 6, 0, -2, 0, 0, 2, 1, -6, 0, 3, 0, 0, -1, 2, -2, 1, -5, 0, 3, 0, 2, 0, 0, -2, 1, 4, 0, -2, 0, 0, 1, 2, -2, 1, 4, 0, -2, 0, 1, 0, 0, -1, 0, -4, 0, 0, 0, 2, 1, 0, -2, 0, 1, 0, 0, 0, 0, 0, -2, 2, 1, 1, 0, 0, 0, 0, 1, -2, 2, 0, -1, 0, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, -1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 2, -2, 0, -1, 0, 0, 0, 0, 0, 2, 0, 2, -2274, -2, 977, -5, 1, 0, 0, 0, 0, 712, 1, -7, 0, 0, 0, 2, 0, 1, -386, -4, 200, 0, 1, 0, 2, 0, 2, -301, 0, 129, -1, 1, 0, 0, -2, 0, -158, 0, -1, 0, -1, 0, 2, 0, 2, 123, 0, -53, 0, 0, 0, 0, 2, 0, 63, 0, -2, 0, 1, 0, 0, 0, 1, 63, 1, -33, 0, -1, 0, 0, 0, 1, -58, -1, 32, 0, -1, 0, 2, 2, 2, -59, 0, 26, 0, 1, 0, 2, 0, 1, -51, 0, 27, 0, 0, 0, 2, 2, 2, -38, 0, 16, 0, 2, 0, 0, 0, 0, 29, 0, -1, 0, 1, 0, 2, -2, 2, 29, 0, -12, 0, 2, 0, 2, 0, 2, -31, 0, 13, 0, 0, 0, 2, 0, 0, 26, 0, -1, 0, -1, 0, 2, 0, 1, 21, 0, -10, 0, -1, 0, 0, 2, 1, 16, 0, -8, 0, 1, 0, 0, -2, 1, -13, 0, 7, 0, -1, 0, 2, 2, 1, -10, 0, 5, 0, 1, 1, 0, -2, 0, -7, 0, 0, 0, 0, 1, 2, 0, 2, 7, 0, -3, 0, 0, -1, 2, 0, 2, -7, 0, 3, 0, 1, 0, 2, 2, 2, -8, 0, 3, 0, 1, 0, 0, 2, 0, 6, 0, 0, 0, 2, 0, 2, -2, 2, 6, 0, -3, 0, 0, 0, 0, 2, 1, -6, 0, 3, 0, 0, 0, 2, 2, 1, -7, 0, 3, 0, 1, 0, 2, -2, 1, 6, 0, -3, 0, 0, 0, 0, -2, 1, -5, 0, 3, 0, 1, -1, 0, 0, 0, 5, 0, 0, 0, 2, 0, 2, 0, 1, -5, 0, 3, 0, 0, 1, 0, -2, 0, -4, 0, 0, 0, 1, 0, -2, 0, 0, 4, 0, 0, 0, 0, 0, 0, 1, 0, -4, 0, 0, 0, 1, 1, 0, 0, 0, -3, 0, 0, 0, 1, 0, 2, 0, 0, 3, 0, 0, 0, 1, -1, 2, 0, 2, -3, 0, 1, 0, -1, -1, 2, 2, 2, -3, 0, 1, 0, -2, 0, 0, 0, 1, -2, 0, 1, 0, 3, 0, 2, 0, 2, -3, 0, 1, 0, 0, -1, 2, 2, 2, -3, 0, 1, 0, 1, 1, 2, 0, 2, 2, 0, -1, 0, -1, 0, 2, -2, 1, -2, 0, 1, 0, 2, 0, 0, 0, 1, 2, 0, -1, 0, 1, 0, 0, 0, 2, -2, 0, 1, 0, 3, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 1, 2, 2, 0, -1, 0, -1, 0, 0, 0, 2, 1, 0, -1, 0, 1, 0, 0, -4, 0, -1, 0, 0, 0, -2, 0, 2, 2, 2, 1, 0, -1, 0, -1, 0, 2, 4, 2, -2, 0, 1, 0, 2, 0, 0, -4, 0, -1, 0, 0, 0, 1, 1, 2, -2, 2, 1, 0, -1, 0, 1, 0, 2, 2, 1, -1, 0, 1, 0, -2, 0, 2, 4, 2, -1, 0, 1, 0, -1, 0, 4, 0, 2, 1, 0, 0, 0, 1, -1, 0, -2, 0, 1, 0, 0, 0, 2, 0, 2, -2, 1, 1, 0, -1, 0, 2, 0, 2, 2, 2, -1, 0, 0, 0, 1, 0, 0, 2, 1, -1, 0, 0, 0, 0, 0, 4, -2, 2, 1, 0, 0, 0, 3, 0, 2, -2, 2, 1, 0, 0, 0, 1, 0, 2, -2, 0, -1, 0, 0, 0, 0, 1, 2, 0, 1, 1, 0, 0, 0, -1, -1, 0, 2, 1, 1, 0, 0, 0, 0, 0, -2, 0, 1, -1, 0, 0, 0, 0, 0, 2, -1, 2, -1, 0, 0, 0, 0, 1, 0, 2, 0, -1, 0, 0, 0, 1, 0, -2, -2, 0, -1, 0, 0, 0, 0, -1, 2, 0, 1, -1, 0, 0, 0, 1, 1, 0, -2, 1, -1, 0, 0, 0, 1, 0, -2, 2, 0, -1, 0, 0, 0, 2, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 2, 4, 2, -1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0)
}
