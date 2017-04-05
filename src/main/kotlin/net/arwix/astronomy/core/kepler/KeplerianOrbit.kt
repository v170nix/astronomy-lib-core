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

import net.arwix.astronomy.core.GM_Sun
import net.arwix.astronomy.core.vector.Matrix
import net.arwix.astronomy.core.vector.RectangularVector

object KeplerianOrbit {

    enum class Planet {SUN, MERCURY, VENUS, EARTH, MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO;

        fun getOrbitalPlane(T: Double): OrbitalPlane {
            val p = 1.3970
            var a = 0.0;
            var e = 0.0;
            var M0 = 0.0;
            var n = 0.0
            var O = 0.0;
            var i = 0.0;
            var w = 0.0;
            var T0 = 0.0
            when (this) {
                KeplerianOrbit.Planet.SUN -> return OrbitalPlane(RectangularVector(), RectangularVector())
                MERCURY -> {
                    a = 0.387099; e = 0.205634; M0 = 174.7947; n = 149472.6738;
                    O = 48.331; i = 7.0048; w = 77.4552; T0 = 0.0;
                }
                KeplerianOrbit.Planet.VENUS -> {
                    a = 0.723332; e = 0.006773; M0 = 50.4071; n = 58517.8149;
                    O = 76.680; i = 3.3946; w = 131.5718; T0 = 0.0;
                }
                KeplerianOrbit.Planet.EARTH -> {
                    a = 1.000000; e = 0.016709; M0 = 357.5256; n = 35999.3720;
                    O = 174.876; i = 0.0000; w = 102.9400; T0 = 0.0;
                }
                KeplerianOrbit.Planet.MARS -> {
                    a = 1.523692; e = 0.093405; M0 = 19.3879; n = 19140.3023;
                    O = 49.557; i = 1.8496; w = 336.0590; T0 = 0.0;
                }
                KeplerianOrbit.Planet.JUPITER -> {
                    a = 5.204267; e = 0.048775; M0 = 18.8185; n = 3033.6272;
                    O = 100.4908; i = 1.3046; w = 15.5576; T0 = 0.0;
                }
                KeplerianOrbit.Planet.SATURN -> {
                    a = 9.582018; e = 0.055723; M0 = 320.3477; n = 1213.8664;
                    O = 113.6427; i = 2.4852; w = 89.6567; T0 = 0.0;
                }
                KeplerianOrbit.Planet.URANUS -> {
                    a = 19.229412; e = 0.044406; M0 = 142.9559; n = 426.9282;
                    O = 73.9893; i = 0.7726; w = 170.5310; T0 = 0.0;
                }
                KeplerianOrbit.Planet.NEPTUNE -> {
                    a = 30.103658; e = 0.011214; M0 = 267.7649; n = 217.9599;
                    O = 131.7942; i = 1.7680; w = 37.4435; T0 = 0.0;
                }
                KeplerianOrbit.Planet.PLUTO -> {
                    a = 39.264230; e = 0.244672; M0 = 15.0233; n = 146.3183;
                    O = 110.2867; i = 17.1514; w = 224.0499; T0 = 0.0;
                }
            }

            return EllipticOrbit.getOrbitalPlane(GM_Sun, Math.toRadians(M0 + n * (T - T0)), a, e).let {
                val PQR = createGaussianMatrix(Math.toRadians(O + p * T), Math.toRadians(i), Math.toRadians(w - O))
                OrbitalPlane(PQR * it.position, PQR * it.velocity)
            }
        }

    }


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