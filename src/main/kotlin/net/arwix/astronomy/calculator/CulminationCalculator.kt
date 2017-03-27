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

import net.arwix.astronomy.core.Epoch
import net.arwix.astronomy.core.MJD_J2000
import net.arwix.astronomy.core.calendar.getGMST
import net.arwix.astronomy.core.calendar.getMJD
import net.arwix.astronomy.core.calendar.resetTime
import net.arwix.astronomy.core.calendar.setHours
import net.arwix.astronomy.core.coordinates.Location
import net.arwix.astronomy.core.vector.SphericalVector
import net.arwix.astronomy.core.vector.Vector
import net.arwix.astronomy.core.vector.VectorType
import net.arwix.astronomy.math.SearchExtremumGoldenSection
import java.lang.Math.cos
import java.lang.Math.sin
import java.util.*


class CulminationCalculator(date: Calendar, location: Location, getGeocentricEquatorialCoordinates: (T: Double, Epoch) -> Vector) :
        BaseEventCalculator(date, location, 0.0, getGeocentricEquatorialCoordinates) {

    sealed class CulminationResult {
        data class Upper(val isAbove: Boolean, val calendar: Calendar) : CulminationResult()
        data class Lower(val isAbove: Boolean, val calendar: Calendar) : CulminationResult()
        data class UpperLover(val upper: Upper, val lower: Lower) : CulminationResult()
    }

    var result: CulminationResult? = null
        get() {
            if (!isValid) result = calls()
            return result
        }


    private fun calls(): CulminationResult {
        isValid = false
        innerDate.resetTime()
        val MJD0 = innerDate.getMJD()
        val cosLatitude = cos(location.latitude)
        val sinLatitude = sin(location.latitude)
        return SearchExtremumGoldenSection(0.0, 24.0, Math.ulp(100.0), 50,
                { x ->
                    val T = (MJD0 + x / 24.0 - MJD_J2000 - deltaT) / 36525.0
                    val p = getGeocentricEquatorialCoordinates(T, Epoch.APPARENT).getVectorOfType(VectorType.SPHERICAL) as SphericalVector
                    // часовой угол
                    val tau = getGMST(MJD0 + x / 24.0) + location.longitude - p.phi
                    return@SearchExtremumGoldenSection sinLatitude * sin(p.theta) + cosLatitude * cos(p.theta) * cos(tau)
                })
                .let {
                    val upperTime = Calendar.getInstance(date.timeZone).apply { time = date.time; setHours(it.max) }
                    val lowerTime = Calendar.getInstance(date.timeZone).apply { time = date.time; setHours(it.min) }
                    isValid = true
                    CulminationResult.UpperLover(
                            CulminationResult.Upper(
                                    getSinAltitude(upperTime.getMJD(), location.longitude, cosLatitude, sinLatitude) - sinRefractionAngle > 0.0,
                                    upperTime
                            ),
                            CulminationResult.Lower(
                                    getSinAltitude(lowerTime.getMJD(), location.longitude, cosLatitude, sinLatitude) - sinRefractionAngle > 0.0,
                                    lowerTime
                            )
                    )
                }
    }

    /**
     * Синус высоты объекта над горизонтом
     * @param MJD         на расчетную дату
     * @param longitude   долгота в радианах
     * @param cosLatitude косинус широты
     * @param sinLatitude синус широты
     * @return cинус высоты Солнца или Луны в момент искомого события
     */
    private fun getSinAltitude(MJD: Double, longitude: Double, cosLatitude: Double, sinLatitude: Double): Double {
        val T = (MJD - MJD_J2000 - deltaT) / 36525.0
        val p = getGeocentricEquatorialCoordinates(T, Epoch.APPARENT).getVectorOfType(VectorType.SPHERICAL) as SphericalVector
        // часовой угол
        val tau = getGMST(MJD) + longitude - p.phi
        return sinLatitude * sin(p.theta) + cosLatitude * cos(p.theta) * cos(tau)
    }
}