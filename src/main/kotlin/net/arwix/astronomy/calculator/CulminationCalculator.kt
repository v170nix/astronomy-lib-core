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
import net.arwix.astronomy.core.calendar.getDeltaT
import net.arwix.astronomy.core.calendar.getMJD
import net.arwix.astronomy.core.calendar.resetTime
import net.arwix.astronomy.core.calendar.setHours
import net.arwix.astronomy.core.coordinates.Location
import net.arwix.astronomy.core.vector.Vector
import net.arwix.astronomy.math.SearchExtremumGoldenMethod
import java.lang.Math.cos
import java.lang.Math.sin
import java.util.*
import java.util.concurrent.TimeUnit

/**
 * @param date дата расчета дня
 * @param location долгота и широта в радианах
 * @param funGeocentricEquatorialCoordinates получение геоцентрических экваториальных координат
 * @param precision точность в долях часа
 */
@Deprecated("use CulminationCalc")
class CulminationCalculator(date: Calendar, location: Location, funGeocentricEquatorialCoordinates: (T: Double, Epoch) -> Vector, private val precision: Double) :
        Calculator<CulminationCalculator.CulminationResult>(date, location, funGeocentricEquatorialCoordinates) {

    sealed class CulminationResult {
        data class Upper(val isAbove: Boolean, val calendar: Calendar) : CulminationResult()
        data class Lower(val isAbove: Boolean, val calendar: Calendar) : CulminationResult()
        data class UpperLower(val upper: Upper, val lower: Lower) : CulminationResult()
    }


    override fun calls(): CulminationResult {
        val innerDate = Calendar.getInstance(date.timeZone).apply { this.timeInMillis = date.timeInMillis }
        deltaT = innerDate.getDeltaT(TimeUnit.DAYS)
        innerDate.resetTime()
        val MJD0 = innerDate.getMJD()
        val cosLatitude = cos(location.latitude)
        val sinLatitude = sin(location.latitude)
        return SearchExtremumGoldenMethod(0.0, 24.0, precision, 50,
                { x -> getSinAltitude(MJD0 + x / 24.0, location.longitude, cosLatitude, sinLatitude) })
                .let {
                    val upperTime = Calendar.getInstance(date.timeZone).apply { time = date.time; setHours(it.max) }
                    val lowerTime = Calendar.getInstance(date.timeZone).apply { time = date.time; setHours(it.min) }
                    CulminationResult.UpperLower(
                            CulminationResult.Upper(
                                    getSinAltitude(upperTime.getMJD(), location.longitude, cosLatitude, sinLatitude) > 0.0,
                                    upperTime
                            ),
                            CulminationResult.Lower(
                                    getSinAltitude(lowerTime.getMJD(), location.longitude, cosLatitude, sinLatitude) > 0.0,
                                    lowerTime
                            )
                    )
                }
    }
}