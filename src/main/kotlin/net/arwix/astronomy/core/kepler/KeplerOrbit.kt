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

import net.arwix.astronomy.core.vector.Matrix

object KeplerOrbit {
//------------------------------------------------------------------------------
//
// GaussVec: computes the transformation matrix from the orbital plane
//           coordinate system to the ecliptic
//
// Input:
//
//   Omega    Longitude of the ascending node of the orbit in [rad]
//   i        Inclination of the orbit to the ecliptic in [rad]
//   omega    Argument of perihelion in [rad]
//
// <return>:  Transformation matrix containing the Gaussian vectors P, Q and R
//
//------------------------------------------------------------------------------

    fun createGaussianMatrix(Omega: Double, i: Double, omega: Double): Matrix
            = Matrix(Matrix.Axis.Z, -Omega) * Matrix(Matrix.Axis.X, -i) * Matrix(Matrix.Axis.Z, -omega)

    //------------------------------------------------------------------------------
//
// Kepler: computes position and velocity vectors for Keplerian orbits w.r.t.
//         the ecliptic
//
// Input:
//
//   GM       Product of gravitational constant and centre mass [AU^3*d^-2]
//   t0       Time of perihelion passage
//   t        Time for calculation
//   q        Perihelion distance in [AU]
//   e        Eccentricity of the orbit
//   PQR      Transformation orbital plane -> ecliptic (Gaussian vectors)
//
// Output:
//
//   r        Heliocentric ecliptical position in [AU]
//   v        Heliocentric ecliptical velocity vector in [AU/d]
//
// Note: t0 and t in Julian centuries since J2000
//
//------------------------------------------------------------------------------
    fun getOrbitalPlane(GM: Double, t0: Double, t: Double, q: Double, e: Double, PQR: Matrix): OrbitalPlane {
        val M0 = 0.1
        val eps = 0.1
        //
        // Variables
        //
        val delta = Math.abs(1.0 - e)
        val invax = delta / q;
        val tau = Math.sqrt(GM) * (t - t0);
        val M = tau * Math.sqrt(invax * invax * invax);
        val orbit = if ((M < M0) && (delta < eps)) ParabolicOrbit.getOrbitalPlane(GM, t0, t, q, e)
        else if (e < 1.0) EllipticOrbit.getOrbitalPlane(GM, M, 1.0 / invax, e)
        else HyperbolicOrbit.getOrbitalPlane(GM, t0, t, 1.0 / invax, e)
        return OrbitalPlane(PQR * orbit.position, PQR * orbit.velocity)
    }

}