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

import net.arwix.astronomy.core.ARCSEC_TO_RAD
import net.arwix.astronomy.math.polynomialSum
import net.arwix.astronomy.math.radians.Radian


class EarthMoonElements(val T: Double) {
    private val T2 = T * T

    /**
     * Mean elongation of moon = D
     */
    val elongation: Radian by lazy {
        val x = doubleArrayOf(1.0722612202445078e+06, 1.6029616009939659e+09).polynomialSum(T) +
                doubleArrayOf(-6.7352202374457519e+00, 6.9492746836058421e-03, -3.702060118571e-005, +2.560078201452e-009, 2.555243317839e-011, -3.207663637426e-013).polynomialSum(T) * T2
        ARCSEC_TO_RAD * x
    }

    /**
     * Mean distance of moon from its ascending node = F
     */
    val ascendingNode: Radian by lazy {
        val x = doubleArrayOf(3.3577951412884740e+05, 1.7395272628437717e+09).polynomialSum(T) +
                doubleArrayOf(-1.3117809789650071e+01, -7.5311878482337989e-04, -2.165750777942e-006, -2.790392351314e-009, 4.189032191814e-011, 4.474984866301e-013).polynomialSum(T) * T2
        ARCSEC_TO_RAD * x
    }
    /**
     * Mean anomaly of sun = l' (J. Laskar)
     */
    val sunAnomaly: Radian by lazy {
        val x = doubleArrayOf(1.2871027407441526e+06, 1.2959658102304320e+08).polynomialSum(T) +
                doubleArrayOf(-5.5281306421783094e-01, 8.7473717367324703e-05, -1.1297037031e-5, -4.77258489e-8, 8.8555011e-11, 4.237343e-13, -3.83508e-15, -1.0390e-17, 1.62e-20).polynomialSum(T) * T2
        ARCSEC_TO_RAD * x
    }

    /**
     * Mean anomaly of moon = l
     */
    val anomaly: Radian by lazy {
        val x = doubleArrayOf(4.8586817465825332e+05, 1.7179159228846793e+09).polynomialSum(T) +
                doubleArrayOf(3.1501359071894147e+01, 5.2099641302735818e-02, -2.536291235258e-004, -2.506365935364e-008, 3.452144225877e-011, -1.755312760154e-012).polynomialSum(T) * T2
        ARCSEC_TO_RAD * x
    }

    /**
     * Mean longitude of moon, re mean ecliptic and equinox of date = L
     */
    val longitude: Radian by lazy {
        val x = doubleArrayOf(7.8593980921052420e+05, 1.7325643720442266e+09).polynomialSum(T) +
                doubleArrayOf(-5.6550460027471399e+00, 6.9017248528380490e-03, -6.073960534117e-005, -1.024222633731e-008, 2.235210987108e-010, 7.200592540556e-014).polynomialSum(T) * T2
        ARCSEC_TO_RAD * x
    }

    /**
     *  Lunar free librations
     *  74.7 years. Denoted W or LA
     */
    val lA: Radian by lazy {
        val x = (-0.112 * T + 1.73655499e6) * T - 389552.81;
        ARCSEC_TO_RAD * (x);
    }

    /**
     * 2.891725 years. Denoted LB
     */
    val LB: Radian by lazy {
        ARCSEC_TO_RAD * (4.48175409e7 * T + 806045.7)
    }

    /**
     * 24.2 years. Denoted P or LC
     */
    val LC: Radian by lazy {
        ARCSEC_TO_RAD * (5.36486787e6 * T - 391702.8)
    }


    /**
     * Precession of the equinox pA
     */
    val precession: Radian by lazy {
        doubleArrayOf(0.0, 5028.791959, 1.105414, 0.000076, -0.0000235316, -1.8055e-8, 1.7451e-10, 1.3095e-12, 2.424e-15, -4.759e-17, -8.66e-20).polynomialSum(T) * ARCSEC_TO_RAD
    }

    /**
     * Usual node term re equinox of date, denoted NA
     */
    val NA: Radian by lazy { longitude - ascendingNode }


    /**
     * Fancy node term, denoted NB.
     * Capital Pi of ecliptic motion (Williams 1994)
     */
    val NB: Radian by lazy {
        val x = (((-0.000004 * T + 0.000026) * T + 0.153382) * T - 867.919986) * T + 629543.967373
        NA + ARCSEC_TO_RAD * (3.24e5 - x) - precession;
    }


}