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

import net.arwix.astronomy.core.PI2
import net.arwix.astronomy.core.modulo
import net.arwix.astronomy.core.vector.SphericalVector
import java.lang.Math.*


object EllipticOrbit {

    //------------------------------------------------------------------------------
//
// EccAnom: computes the eccentric anomaly for elliptic orbits
//
// Input:
//
//   M        Mean anomaly in [rad]
//   e        Eccentricity of the orbit [0,1[
//
// <return>:  Eccentric anomaly in [rad]
//
//------------------------------------------------------------------------------
    fun getEccentricAnomaly(M: Double, e: Double): Double {
        //
        // Constants
        //
        val maxit = 15
        val eps = Math.ulp(100.0);


        //
        // Variables
        //
        var i = 0
        var E: Double;
        var f: Double
        // Starting value
        val M = M.modulo(PI2)
        if (e < 0.8) E = M; else E = Math.PI

        // Iteration
        do {
            f = E - e * sin(E) - M
            E -= f / (1.0 - e * cos(E))
            ++i
            if (i == maxit) {
                throw IndexOutOfBoundsException(" Convergence problems in EccAnom")
            }
        } while (abs(f) > eps)

        return E
    }

    //------------------------------------------------------------------------------
//
// Ellip: computes position and velocity vectors for elliptic orbits
//
// Input:
//
//   GM       Product of gravitational constant and centre mass [AU^3*d^-2]
//   M        Mean anomaly in [rad]
//   a        Semi-major axis of the orbit in [AU]
//   e        Eccentricity of the orbit (<1)
//
// Output:
//
//   r        Position w.r.t. orbital plane in [AU]
//   v        Velocity w.r.t. orbital plane in [AU/d]
//
//------------------------------------------------------------------------------
    fun getOrbitalPlane(GM: Double, M: Double, a: Double, e: Double): OrbitalPlane {
        val k = sqrt(GM / a)
        val E = getEccentricAnomaly(M, e)
        val cosE = cos(E)
        val sinE = sin(E)
        val fac = sqrt((1.0 - e) * (1.0 + e))
        val rho = 1.0 - e * cosE
        return OrbitalPlane(SphericalVector(a * (cosE - e), a * fac * sinE, 0.0), SphericalVector(-k * sinE / rho, k * fac * cosE / rho, 0.0))
    }
}