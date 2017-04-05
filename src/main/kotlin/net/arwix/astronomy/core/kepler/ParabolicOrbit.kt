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

import net.arwix.astronomy.core.vector.RectangularVector
import java.lang.Math.*

object ParabolicOrbit {


    //------------------------------------------------------------------------------
//
// Stumpff: computes values for the Stumpff functions C1, C2 and C3
//
// Input:
//
//   E2       Square of eccentric anomaly (E2=E*E) in [rad^2]
//
// Output:
//
//   c1       Value of C1 = sin(E)/E
//   c2       Value of C2 = (1-cos(E))/(E*E)
//   c3       Value of C3 = (E-sin(E))/(E^3)
//
//------------------------------------------------------------------------------
    private fun Stumpff(E2: Double): Triple<Double, Double, Double> {
        val eps = Math.ulp(100.0)
        var c1 = 0.0;
        var c2 = 0.0;
        var c3 = 0.0
        var add = 1.0;
        var n = 1.0
        do {
            c1 += add; add /= (2.0 * n)
            c2 += add; add /= (2.0 * n + 1.0)
            c3 += add; add *= -E2
            n += 1.0
        } while (abs(add) >= eps)
        return Triple(c1, c2, c3)
    }

    //------------------------------------------------------------------------------
//
// Parab: computes position and velocity vectors for parabolic and near
//        parabolic orbits
//
// Input:
//
//   GM       Product of gravitational constant and centre mass [AU^3*d^-2]
//   t0       Time of perihelion passage
//   t        Time for calculation
//   q        Perihelion distance in [AU]
//   e        Eccentricity of the orbit (~1)
//
// Output:
//
//   r        Position w.r.t. orbital plane in [AU]
//   v        Velocity w.r.t. orbital plane in [AU/d]
//
// Note: t0 and t in Julian centuries since J2000
//
//------------------------------------------------------------------------------
    fun getOrbitalPlane(GM: Double, t0: Double, t: Double, q: Double, e: Double): OrbitalPlane {
        //
        // Constants
        //
        val maxit = 15
        val eps = Math.ulp(100.0)

        //
        // Variables
        //
        var i = 0
        var E20: Double
        var E2 = 0.0
        var u: Double
        var u2: Double
        var c: Triple<Double, Double, Double>
        var fac = 0.5 * e

        val k = sqrt(GM / (q * (1.0 + e)));
        val tau = sqrt(GM) * (t - t0);

        do {
            ++i;
            E20 = E2
            val A = 1.5 * sqrt(fac / (q * q * q)) * tau
            val B = pow(sqrt(A * A + 1.0) + A, 1.0 / 3.0)
            u = B - 1.0 / B
            u2 = u * u
            E2 = u2 * (1.0 - e) / fac
            c = Stumpff(E2)
            fac = 3.0 * e * c.third
            if (i == maxit) {
                throw IndexOutOfBoundsException(" Convergence problems in Parab")
            }
        } while (abs(E2 - E20) >= eps)

        val R = q * (1.0 + u2 * c.second * e / fac)
        val r = RectangularVector(q * (1.0 - u2 * c.second / fac), q * sqrt((1.0 + e) / fac) * u * c.first, 0.0)
        return OrbitalPlane(r, RectangularVector(-k * r.y / R, k * (r.x / R + e), 0.0))
    }

}