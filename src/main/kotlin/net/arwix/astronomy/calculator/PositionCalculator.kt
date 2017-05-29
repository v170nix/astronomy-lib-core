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
import net.arwix.astronomy.core.kepler.KeplerBodySimonJ2000
import net.arwix.astronomy.core.kepler.KeplerElements
import net.arwix.astronomy.core.vector.RectangularVector
import net.arwix.astronomy.core.vector.Vector
import net.arwix.astronomy.ephemeris.Precession


@Geocentric @Equatorial
class PositionCalculator(
        val precession: Precession,
        @Heliocentric @Ecliptic val getEarthCoordinatesJ2000: FunGetHeliocentricEclipticCoordinates) {

    val obliquity = if (precession.isEcliptic) precession.getNearestObliquityModel() else precession.getNearestObliquityModel(0.0)
    val nutation = precession.getNearsetNutationModel(obliquity)

    val transformMatrix = if (precession.isEcliptic)
        nutation.equatorialMatrix * obliquity.eclipticToEquatorialMatrix * precession.fromJ2000Matrix
    else
        nutation.equatorialMatrix * precession.fromJ2000Matrix * obliquity.eclipticToEquatorialMatrix
    /**
     * @param lightTime in days
     */
    @Geocentric @Ecliptic
    private fun getGeoPositionJ2000(T: Double, request: Request, lightTime: Double = 0.0): Vector {
        val dT = lightTime / JULIAN_DAYS_PER_CENTURY
        return when (request) {
            is Request.GeocentricEclipticBody -> request.getCoordinates(T - dT)
            is Request.HeliocentricEclipticBody -> {
                val body = request.getCoordinates(T - dT)
                val earth = getEarthCoordinatesJ2000(T)
                body - earth
            }
        }
    }

    @Heliocentric @Ecliptic
    fun getHeliocentricEclipticPositionJ2000(T: Double, request: Request) = when (request) {
        is Request.HeliocentricEclipticBody -> request.getCoordinates(T)
        is Request.GeocentricEclipticBody -> getEarthCoordinatesJ2000(T) + request.getCoordinates(T)
    }


    @Heliocentric @Ecliptic
    fun getBodyVelocity(T: Double, request: Request, lightTime: Double = 0.01): Vector {
        val body = getHeliocentricEclipticPositionJ2000(T, request)
        val bodyPlus = getHeliocentricEclipticPositionJ2000(T + lightTime / JULIAN_DAYS_PER_CENTURY, request)
        // скорость объекта
        return (bodyPlus - body) / lightTime
    }

    @Heliocentric @Ecliptic
    fun getEarthVelocity(T: Double): Vector {
        val earth = getEarthCoordinatesJ2000(T)
        val earthPlus = getEarthCoordinatesJ2000(T + 0.01 / JULIAN_DAYS_PER_CENTURY)
        //скорость Земли
        return (earthPlus - earth) / .01
    }

    @Geocentric @Equatorial
    fun getGeocentricEquatorialPositionGeometric(T: Double, request: Request): Vector {
        return obliquity.rotateFromEclipticToEquatorial(getGeoPositionJ2000(T, request, 0.0))
    }


    @Geocentric @Equatorial
    fun getGeocentricEquatorialPositionApparent(T: Double, request: Request): Vector {
        val bodyGeo = getGeoPositionJ2000(T, request, 0.0)
        val lightTime = bodyGeo.normalize() / C_Light
        val bodyElements = if (request is Request.HeliocentricEclipticBody && request.keplerElements != null) request.keplerElements else null
        val earthVelocity = KeplerBodySimonJ2000.Earth.getOrbitalPlane(T).velocity

        val oldBodyVelocity = request.velocitySpeed

        val bodyVelocity = if (request.usePreviousVelocitySpeed && oldBodyVelocity != null) oldBodyVelocity else {
            bodyElements?.getOrbitalPlane(T)?.velocity ?: getBodyVelocity(T, request, lightTime)
        }
        request.velocitySpeed = bodyVelocity
        return transformMatrix * (bodyGeo - (bodyVelocity - earthVelocity) * lightTime)
    }

    /**
     * @param light_time in days
     */
    @Deprecated("old function")
    private fun aberration(pObject: RectangularVector, vearth: RectangularVector, light_time: Double): RectangularVector {
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

    @Ecliptic
    sealed class Request(val usePreviousVelocitySpeed: Boolean = false) {

        internal var velocitySpeed: Vector? = null

        @Heliocentric @Ecliptic
        class HeliocentricEclipticBody(
                val getCoordinates: FunGetHeliocentricEclipticCoordinates,
                val keplerElements: KeplerElements? = null,
                usePreviousVelocitySpeed: Boolean = false) : Request(usePreviousVelocitySpeed)

        @Geocentric @Ecliptic
        class GeocentricEclipticBody(
                val getCoordinates: FunGetGeocentricEclipticCoordinates,
                usePreviousVelocitySpeed: Boolean = false) : Request(usePreviousVelocitySpeed)
    }


}