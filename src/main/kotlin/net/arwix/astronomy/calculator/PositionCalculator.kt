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

package net.arwix.astronomy.calculator

import net.arwix.astronomy.annotation.Ecliptic
import net.arwix.astronomy.annotation.Equatorial
import net.arwix.astronomy.annotation.Geocentric
import net.arwix.astronomy.annotation.Heliocentric
import net.arwix.astronomy.core.C_Light
import net.arwix.astronomy.core.JULIAN_DAYS_PER_CENTURY
import net.arwix.astronomy.core.coordinates.FunGetGeocentricEclipticCoordinates
import net.arwix.astronomy.core.coordinates.FunGetHeliocentricEclipticCoordinates
import net.arwix.astronomy.core.kepler.KeplerBodyJPL
import net.arwix.astronomy.core.vector.RectangularVector
import net.arwix.astronomy.core.vector.Vector
import net.arwix.astronomy.core.vector.VectorType
import net.arwix.astronomy.ephemeris.Precession

@Geocentric @Equatorial
class PositionCalculator(
        val precession: Precession,
        @Heliocentric @Ecliptic
        private val getEarthHelioEclipticCoordinatesJ2000: FunGetHeliocentricEclipticCoordinates) {

    val obliquity = if (precession.isEcliptic) precession.getNearestObliquityModel() else precession.getNearestObliquityModel(0.0)
    val nutation = precession.getNearsetNutationModel(obliquity)


    @Geocentric @Equatorial
    fun getApparentGeocentricEquatorialPosition(T: Double, getHelioEclipticCoordinatesJ2000: FunGetHeliocentricEclipticCoordinates): Vector {
        @Heliocentric @Ecliptic val earthPosition = getEarthHelioEclipticCoordinatesJ2000(T)
        @Heliocentric @Ecliptic var bodyPosition = getHelioEclipticCoordinatesJ2000(T)
        @Geocentric @Ecliptic var bodyGeoPosition: Vector = bodyPosition - earthPosition
        val lightTime = bodyGeoPosition.normalize() / C_Light

        bodyPosition = getHelioEclipticCoordinatesJ2000(T - lightTime / JULIAN_DAYS_PER_CENTURY)
        bodyGeoPosition = bodyPosition - earthPosition

        if (precession.isEcliptic) bodyGeoPosition = precession.transformFromJ2000(bodyGeoPosition)
        val eOrbit = KeplerBodyJPL.EarthMoonBarycenter(T).getOrbitalPlane()
        bodyGeoPosition = aberration(bodyGeoPosition.getVectorOfType(VectorType.RECTANGULAR) as RectangularVector, eOrbit.velocity as RectangularVector, lightTime)

        @Geocentric @Equatorial var bodyGeoEquPosition = obliquity.rotateFromEclipticToEquatorial(bodyGeoPosition)

        if (!precession.isEcliptic) bodyGeoEquPosition = precession.transformFromJ2000(bodyGeoEquPosition)

        bodyGeoEquPosition = nutation.applyNutationToGeocentricVector(bodyGeoEquPosition)

        return bodyGeoEquPosition
    }

    fun getMoon(T: Double, getGeocentricEclipticCoordinates: FunGetGeocentricEclipticCoordinates): Vector {
        @Geocentric @Ecliptic var bodyGeoPosition = getGeocentricEclipticCoordinates(T)
        val lightTime = bodyGeoPosition.normalize() / C_Light

        bodyGeoPosition = getGeocentricEclipticCoordinates(T - lightTime / JULIAN_DAYS_PER_CENTURY)

        if (precession.isEcliptic) bodyGeoPosition = precession.transformFromJ2000(bodyGeoPosition)
        val eOrbit = KeplerBodyJPL.EarthMoonBarycenter(T).getOrbitalPlane()
        //     bodyGeoPosition = aberration(bodyGeoPosition.getVectorOfType(VectorType.RECTANGULAR) as RectangularVector, eOrbit.velocity as RectangularVector, lightTime)

        @Geocentric @Equatorial var bodyGeoEquPosition = obliquity.rotateFromEclipticToEquatorial(bodyGeoPosition)

        if (!precession.isEcliptic) bodyGeoEquPosition = precession.transformFromJ2000(bodyGeoEquPosition)

        bodyGeoEquPosition = nutation.applyNutationToGeocentricVector(bodyGeoEquPosition)

        return bodyGeoEquPosition
    }

    fun aberration(pObject: RectangularVector, vearth: RectangularVector, light_time: Double): RectangularVector {
        if (light_time <= 0) return pObject

        //    val vearth = doubleArrayOf(earth[3], earth[4], earth[5])
        val p = DoubleArray(3)

        val TL = light_time
        val P1MAG = TL * C_Light
        val VEMAG = vearth.normalize()
        if (VEMAG == 0.0) return pObject
        val BETA = VEMAG / C_Light
        val DOT = pObject[0] * vearth.x + pObject[1] * vearth.y + pObject[2] * vearth.z
        val COSD = DOT / (P1MAG * VEMAG)
        val GAMMAI = Math.sqrt(1.0 - BETA * BETA)
        val P = BETA * COSD
        val Q = (1.0 + P / (1.0 + GAMMAI)) * TL
        val R = 1.0 + P

        for (i in 0..2) {
            p[i] = (GAMMAI * pObject[i] + Q * vearth[i]) / R
        }

        return RectangularVector(p[0], p[1], p[2])
    }


}