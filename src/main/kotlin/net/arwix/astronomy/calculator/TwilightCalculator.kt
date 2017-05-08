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
import net.arwix.astronomy.core.coordinates.Location
import net.arwix.astronomy.core.vector.Vector
import java.lang.Math.sin
import java.lang.Math.toRadians
import java.util.*

class TwilightCalculator(date: Calendar, location: Location, funGeocentricEquatorialCoordinates: (T: Double, Epoch) -> Vector, type: TwilightType) :
        Calculator<TwilightCalculator.TwilightResult>(date, location, funGeocentricEquatorialCoordinates) {

    sealed class TwilightResult(open val type: TwilightType) {
        data class Begin(override val type: TwilightType, val calendar: Calendar) : TwilightResult(type)
        data class End(override val type: TwilightType, val calendar: Calendar) : TwilightResult(type)
        data class BeginEnd(override val type: TwilightType, val begin: Calendar, val end: Calendar) : TwilightResult(type)
        data class None(override val type: TwilightType, val isAbove: Boolean) : TwilightResult(type)
    }

    enum class TwilightType constructor(angle: Double) {
        CIVIL(-6.0), NAUTICAL(-12.0), ASTRONOMICAL(-18.0);
        internal val sinRefractionAngle = sin(toRadians(angle))
    }

    private var innerCalculator: RiseSetCalculator = RiseSetCalculator(date, location, getCoordinates, type.sinRefractionAngle)

    var type: TwilightType = type
        set(value) {
            if (value != field) {
                isValidResult = false
                innerCalculator.sinRefractionAngle = value.sinRefractionAngle
                field = value
            }
        }

    override fun calls(): TwilightResult {
        if (innerCalculator.date.timeInMillis != date.timeInMillis || innerCalculator.location != location) {
            innerCalculator = RiseSetCalculator(date, location, getCoordinates, type.sinRefractionAngle)
        }
        val innerResult = innerCalculator.getResult()
        return when (innerResult) {
            is RiseSetCalculator.Result.RiseSet -> TwilightResult.BeginEnd(type, innerResult.rise.calendar, innerResult.set.calendar)
            is RiseSetCalculator.Result.Rise -> TwilightResult.Begin(type, innerResult.calendar)
            is RiseSetCalculator.Result.Set -> TwilightResult.End(type, innerResult.calendar)
            is RiseSetCalculator.Result.None -> TwilightResult.None(type, innerResult.isAbove)
        }
    }


}