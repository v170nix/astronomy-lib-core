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
import net.arwix.astronomy.core.coordinates.Coordinates
import net.arwix.astronomy.core.coordinates.FunGetGeocentricEclipticCoordinates
import net.arwix.astronomy.core.coordinates.FunGetHeliocentricEclipticCoordinates
import net.arwix.astronomy.core.kepler.EarthMoonElements
import net.arwix.astronomy.core.kepler.KeplerObjectSimonJ2000
import net.arwix.astronomy.core.vector.Matrix
import net.arwix.astronomy.core.vector.RectangularVector
import net.arwix.astronomy.core.vector.VectorType
import net.arwix.astronomy.ephemeris.Precession


sealed class SwissObject : Coordinates {

    class Libration : SwissObject(), Coordinates by SwissLibrationImpl
    class MERCURY : SwissObject(), Coordinates by SwissBaseImpl(MercurySwissData)
    class VENUS : SwissObject(), Coordinates by SwissBaseImpl(VenusSwissData)
    class Earth_Moon_Barycenter : SwissObject(), Coordinates by SwissEarthMoonBarycenterImpl
    class EARTH : SwissObject(), Coordinates by SwissBaseImpl(EarthMoonBarycenterSwissData)
    class Moon(precession: Precession.WILLIAMS_1994) : SwissObject(), Coordinates by SwissMoonImpl(precession)
    class MARS : SwissObject(), Coordinates by SwissBaseImpl(MarsSwissData)
    class JUPITER : SwissObject(), Coordinates by SwissBaseImpl(JupiterSwissData)
    class SATURN : SwissObject(), Coordinates by SwissBaseImpl(SaturnSwissData)
    class URANUS : SwissObject(), Coordinates by SwissBaseImpl(UranusSwissData)
    class NEPTUNE : SwissObject(), Coordinates by SwissBaseImpl(NeptuneSwissData)
    class Pluto : SwissObject(), Coordinates by SwissBaseImpl(PlutoSwissData)
}

private class SwissBaseImpl(val swissData: SwissData) : Coordinates {
    override val getGeocentricEclipticCoordinates: FunGetGeocentricEclipticCoordinates get() = TODO("not implemented")
    override val getGeocentricEquatorialCoordinates: FunGetGeocentricEclipticCoordinates get() = TODO("not implemented")

    override val getHeliocentricEclipticCoordinates: FunGetHeliocentricEclipticCoordinates
        get() = { t ->
            RectangularVector(
                    gplan(t, swissData.args, swissData.distance, swissData.tabb, swissData.tabl,
                            swissData.tabr, swissData.max_harmonic, swissData.max_power_of_t,
                            swissData.maxargs, swissData.timescale, swissData.trunclvl)
            )
        }
}

private object SwissEarthMoonBarycenterImpl : Coordinates {
    override val getGeocentricEclipticCoordinates: FunGetGeocentricEclipticCoordinates get() = TODO("not implemented")
    override val getGeocentricEquatorialCoordinates: FunGetGeocentricEclipticCoordinates get() = TODO("not implemented")

    override val getHeliocentricEclipticCoordinates: FunGetHeliocentricEclipticCoordinates
        get() = { t ->
            RectangularVector(g3plan(t, EarthMoonBarycenterSwissData.args, EarthMoonBarycenterSwissData.distance,
                    EarthMoonBarycenterSwissData.tabb, EarthMoonBarycenterSwissData.tabl,
                    EarthMoonBarycenterSwissData.tabr, EarthMoonBarycenterSwissData.max_harmonic,
                    EarthMoonBarycenterSwissData.max_power_of_t, EarthMoonBarycenterSwissData.maxargs,
                    EarthMoonBarycenterSwissData.timescale, EarthMoonBarycenterSwissData.trunclvl, false))

        }
}

private class SwissMoonImpl(val precession: Precession.WILLIAMS_1994) : Coordinates {
    override val getGeocentricEclipticCoordinates: FunGetGeocentricEclipticCoordinates get() = TODO("not implemented")
    override val getGeocentricEquatorialCoordinates: FunGetGeocentricEclipticCoordinates get() = TODO("not implemented")

    override val getHeliocentricEclipticCoordinates: FunGetHeliocentricEclipticCoordinates
        get() = { t ->
            val moon_lat = g1plan(t, MoonLatSwissData.args, MoonLatSwissData.tabl, MoonLatSwissData.max_harmonic, MoonLatSwissData.timescale)
            val vector = g2plan(t, MoonLonSwissData.args, MoonLonSwissData.distance,
                    MoonLonSwissData.tabl, MoonLonSwissData.tabr, MoonLonSwissData.max_harmonic, MoonLonSwissData.timescale, moon_lat)
            // Here we apply Williams formula to pass to J2000, since this is
            // the one chosen by Moshier
//            var vector: Vector = RectangularVector(p)
//            val precession = Precession.WILLIAMS_1994(t)
            precession.transformToJ2000(vector)
            //           var velocityVector = RectangularVector(0.0, 0.0, 0.0)
        }

}

private object SwissLibrationImpl : Coordinates {
    override val getGeocentricEclipticCoordinates: FunGetGeocentricEclipticCoordinates get() = TODO("not implemented")
    override val getGeocentricEquatorialCoordinates: FunGetGeocentricEclipticCoordinates get() = TODO("not implemented")

    override val getHeliocentricEclipticCoordinates: FunGetHeliocentricEclipticCoordinates
        get() = { t ->
            var p = g3plan(t, LibrationSwissData.args, LibrationSwissData.distance,
                    LibrationSwissData.tabb, LibrationSwissData.tabl,
                    LibrationSwissData.tabr, LibrationSwissData.max_harmonic,
                    LibrationSwissData.max_power_of_t, LibrationSwissData.maxargs,
                    LibrationSwissData.timescale, LibrationSwissData.trunclvl, true).getVectorOfType(VectorType.RECTANGULAR) as RectangularVector

            p[0] -= KeplerObjectSimonJ2000.Earth(t).Longitude - .047 * ARCSEC_TO_RAD

            val elements = EarthMoonElements(t)

            val LP_equinox = elements.longitude
            val NF_arcsec = elements.ascendingNode

            // phi+psi
            p[2] += LP_equinox + 6.48e5 * ARCSEC_TO_RAD
            if (p[2] < -6.45e5 * ARCSEC_TO_RAD) p[2] += 1.296e6 * ARCSEC_TO_RAD
            if (p[2] > 6.45e5 * ARCSEC_TO_RAD) p[2] -= 1.296e6 * ARCSEC_TO_RAD

            // phi
            p[0] += LP_equinox - NF_arcsec + 6.48e5 * ARCSEC_TO_RAD
            if (p[0] < -6.45e5 * ARCSEC_TO_RAD) p[0] += 1.296e6 * ARCSEC_TO_RAD
            if (p[0] > 6.45e5 * ARCSEC_TO_RAD) p[0] -= 1.296e6 * ARCSEC_TO_RAD
            p[2] -= p[0]

            // From Euler angles to matrix M
            val x = Matrix(Matrix.Axis.Z, p[2])
            val y = Matrix(Matrix.Axis.X, p[1])
            val z = Matrix(Matrix.Axis.Z, p[0])
            val mM = x * y * z

            // Get rotation matrix around x axis with eps
            val mQ = Matrix(Matrix.Axis.X, 84381.406173 * ARCSEC_TO_RAD)

            // Get precession matrix
            var mP = Precession.JPL_DE4xx(t).fromJ2000Matrix

            // Precess Q
            mP = mP * mQ

            // Space to body
            val mM2000 = mM * mP

            // Get back the Euler angles, now equatorial (as JPL ones)
            val M = mM2000.toArray() // mM2000.array
            val phi = Math.atan2(M[2][0], -M[2][1])
            var a = M[0][2]
            val b = M[1][2]
            /* psi = zatan2( b, a ); */
            val psi = Math.atan2(a, b)

            if (Math.abs(a) > Math.abs(b)) a /= Math.sin(psi) else a = b / Math.cos(psi)

            /* theta = zatan2( M[2][2], a ); */
            val theta = Math.atan2(a, M[2][2])

            p[0] = phi
            p[1] = theta
            p[2] = psi
            RectangularVector(p)
        }
}
