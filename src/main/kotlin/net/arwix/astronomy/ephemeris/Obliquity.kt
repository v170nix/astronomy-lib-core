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
import net.arwix.astronomy.core.MINUTES_PER_DEGREE
import net.arwix.astronomy.core.PI2
import net.arwix.astronomy.core.SECONDS_PER_DEGREE
import net.arwix.astronomy.ephemeris.precession.PrecessionMethod
import java.lang.Math.cos
import java.lang.Math.sin


object Obliquity {

    private val xypol = doubleArrayOf(84028.206305, 0.3624445, -0.00004039, -110E-9)

    private val xyper = arrayOf(doubleArrayOf(409.90, 753.872780, -1704.720302),
            doubleArrayOf(396.15, -247.805823, -862.308358),
            doubleArrayOf(537.22, 379.471484, 447.832178),
            doubleArrayOf(402.90, -53.880558, -889.571909),
            doubleArrayOf(417.15, -90.109153, 190.402846),
            doubleArrayOf(288.92, -353.600190, -56.564991),
            doubleArrayOf(4043.00, -63.115353, -296.222622),
            doubleArrayOf(306.00, -28.248187, -75.859952),
            doubleArrayOf(277.00, 17.703387, 67.473503),
            doubleArrayOf(203.00, 38.911307, 3.014055))


    /**
     * Calculates the mean obliquity at a given time. The code for this method
     * comes from S. L. Moshier. In case IAU2006 or IAU2009 algorithms are used
     * and the Vondrak precession flag is enabled in the configuration, the obliquity
     * will be calculated using the formulae by Vondrak. Time span in Julian centuries
     * from J2000 can be expanded in this case from +/- 100 to +/- 2000.<P>
     * References:
     * Vondrak et al. 2011, A&amp;A 534, A22.
     *
     * Capitaine et al., Astronomy and Astrophysics 412, 567-586, (2003).
     *
     * James G. Williams, "Contributions to the Earth's obliquity rate,
     * precession, and nutation," Astron. J. 108, 711-724 (1994).
     *
     * J. L. Simon, P. Bretagnon, J. Chapront, M. Chapront-Touze', G. Francou,
     * and J. Laskar, "Numerical Expressions for precession formulae and mean
     * elements for the Moon and the planets," Astronomy and Astrophysics 282,
     * 663-683 (1994).
     *
     * J. H. Lieske, T. Lederle, W. Fricke, and B. Morando, "Expressions for the
     * Precession Quantities Based upon the IAU (1976) System of Astronomical
     * Constants," Astronomy and Astrophysics 58, 1-16 (1977).
     *
     * J. Laskar's expansion comes from "Secular terms of classical planetary
     * theories using the results of general theory," Astronomy and Astrophysics
     * 157, 59070 (1986).
     *
     * @param t Time in Julian centuries from J2000 in TT.
     *          Valid range is the years -8000 to +12000 (t = -100 to 100). In case of using
     *          the Vondrak formulae validity is t between +/- 2000.
     * @return The mean obliquity (epsilon sub 0) in radians.
     */
    fun meanObliquity(t: Double, method: PrecessionMethod = PrecessionMethod.VONDRAK_2011): Double {

        // The obliquity formula come from Meeus, Astro Algorithms, 2ed.
        var rval = 0.0
        var u: Double
        val u0: Double

        // Capitaine et al. 2003, Hilton et al. 2006
        val rvalStart_CAP = 23.0 * SECONDS_PER_DEGREE + 26.0 * MINUTES_PER_DEGREE + 21.406
        val coeffs_CAP = doubleArrayOf(-468367.69, -183.1, 200340.0, -5760.0, -43400.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        // Simon et al., 1994
        val OBLIQ_COEFFS = 10
        val rvalStart_SIM = 23.0 * SECONDS_PER_DEGREE + 26.0 * MINUTES_PER_DEGREE + 21.412
        val coeffs_SIM = doubleArrayOf(-468092.7, -152.0, 199890.0, -5138.0, -24967.0, -3905.0, 712.0, 2787.0, 579.0, 245.0)

        // Williams et al., DE403 Ephemeris
        val rvalStart_WIL = 23.0 * SECONDS_PER_DEGREE + 26.0 * MINUTES_PER_DEGREE + 21.406173
        val coeffs_WIL = doubleArrayOf(-468339.6, -175.0, 199890.0, -5138.0, -24967.0, -3905.0, 712.0, 2787.0, 579.0, 245.0)

        // Laskar et al.
        /*
		 * This expansion is from Laskar, cited above. Bretagnon and Simon say,
		 * in Planetary Programs and Tables, that it is accurate to 0.1" over a
		 * span of 6000 years. Laskar estimates the precision to be 0.01" after
		 * 1000 years and a few seconds of arc after 10000 years.
		 */
        val rvalStart_LAS = 23.0 * SECONDS_PER_DEGREE + 26.0 * MINUTES_PER_DEGREE + 21.448
        val coeffs_LAS = doubleArrayOf(-468093.0, -155.0, 199925.0, -5138.0, -24967.0, -3905.0, 712.0, 2787.0, 579.0, 245.0)

        // IAU 1976
        val rvalStart_IAU = 23.0 * SECONDS_PER_DEGREE + 26.0 * MINUTES_PER_DEGREE + 21.448
        val coeffs_IAU = doubleArrayOf(-468150.0, -590.0, 181300.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        // Select the desired formula
        val (rvalStart, coeffs) = when (method) {
            PrecessionMethod.WILLIAMS_1994, PrecessionMethod.JPL_DE4xx -> Pair(rvalStart_WIL, coeffs_WIL)
            PrecessionMethod.SIMON_1994 -> Pair(rvalStart_SIM, coeffs_SIM)
            PrecessionMethod.LASKAR_1986 -> Pair(rvalStart_LAS, coeffs_LAS)
            PrecessionMethod.IAU_1976 -> Pair(rvalStart_IAU, coeffs_IAU)
            else -> Pair(rvalStart_CAP, coeffs_CAP)
        }


        if (Math.abs(t) > 100 || method == PrecessionMethod.VONDRAK_2011) {
            var y = 0.0
            var w = PI2 * t

            for (i in 0..9) {
                val a = w / xyper[i][0]
                y += cos(a) * xyper[i][1] + sin(a) * xyper[i][2]
            }

            w = 1.0
            for (j in 0..3) {
                y += xypol[j] * w
                w *= t
            }
            rval = y * ARCSEC_TO_RAD

            if (Math.abs(t) > 100) TODO("Date is too far from J2000, obliquity forced to Vondrk et al. 2011 model.")
        } else {
            u0 = t / 100.0
            u = u0 // u is in julian 10000's of years
            rval = rvalStart

            for (i in 0..OBLIQ_COEFFS - 1) {
                rval += u * coeffs[i] / 100.0
                u *= u0
            }

            // convert from seconds to radians
            rval *= ARCSEC_TO_RAD

            if (Math.abs(t) > 100.0) TODO("This date is too far from J2000 epoch. Obliquity is probably incorrect.")
        }

        return rval
    }
}