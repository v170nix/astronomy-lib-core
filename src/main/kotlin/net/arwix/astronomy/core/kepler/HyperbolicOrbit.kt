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


object HyperbolicOrbit {

    //------------------------------------------------------------------------------
//
// HypAnom: computes the eccentric anomaly for hyperbolic orbits
//
// Input:
//
//   Mh       Mean anomaly in [rad]
//   e        Eccentricity of the orbit (>1)
//
// <return>:  Eccentric anomaly in [rad]
//
//------------------------------------------------------------------------------
    fun getAnomaly(Mh: Double, e: Double): Double {
        //
        // Constants
        //
        val maxit = 15
        val eps = Math.ulp(100.0)
        //
        // Variables
        //
        var i = 0
        var f: Double
        // Starting value
        var H = log(2.0 * abs(Mh) / e + 1.8)
        if (Mh < 0.0) H = -H
        // Iteration
        do {
            f = e * sinh(H) - H - Mh
            H -= f / (e * cosh(H) - 1.0)
            ++i
            if (i == maxit) {
                throw IndexOutOfBoundsException(" Convergence problems in HypAnom")
            }
        } while (abs(f) > eps * (1.0 + abs(H + Mh)))
        return H
    }

    //------------------------------------------------------------------------------
//
// Hyperb: computes position and velocity vectors for elliptic orbits
//
// Input:
//
//   GM       Product of gravitational constant and centre mass [AU^3*d^-2]
//   t0       Time of perihelion passage
//   t        Time for calculation
//   a        Semimajor axis of the orbit in [AU]
//   e        Eccentricity of the orbit (>1)
//
// Output:
//
//   r        Position w.r.t. orbital plane in [AU]
//   v        Velocity w.r.t. orbital plane in [AU/d]
//
// Note: t0 and t in Julian centuries since J2000
//
//------------------------------------------------------------------------------
    fun getOrbitalPlane(GM: Double, t0: Double, t: Double, a: Double, e: Double): OrbitalPlane {
        //
        // Variables
        //
        val a = abs(a)
        val k = sqrt(GM / a)

        val Mh = k * (t - t0) / a
        val H = getAnomaly(Mh, e)
        val coshH = cosh(H)
        val sinhH = sinh(H)
        val fac = sqrt((e + 1.0) * (e - 1.0))
        val rho = e * coshH - 1.0


        return OrbitalPlane(RectangularVector(a * (e - coshH), a * fac * sinhH, 0.0), RectangularVector(-k * sinhH / rho, k * fac * coshH / rho, 0.0))
    }

}