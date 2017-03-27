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
import net.arwix.astronomy.core.coordinates.Location
import net.arwix.astronomy.core.vector.Vector
import java.util.*
import java.util.concurrent.TimeUnit


open class BaseEventCalculator {

    protected lateinit var innerDate: Calendar
    protected var deltaT: Double = 0.0
    protected var isValid = false
    protected val getGeocentricEquatorialCoordinates: (T: Double, Epoch) -> Vector
    protected val sinRefractionAngle: Double

    var date: Calendar
        set(value) {
            isValid = false
            innerDate.timeZone = value.timeZone
            innerDate.time = value.time
            deltaT = innerDate.getDeltaT(TimeUnit.DAYS)
            field = value
        }

    var location: Location
        set(value) {
            this.isValid = false
            field = value
        }

    internal constructor(date: Calendar, location: Location, sinRefractionAngle: Double, getGeocentricEquatorialCoordinates: (T: Double, Epoch) -> Vector) {
        this.date = date
        this.location = location
        this.sinRefractionAngle = sinRefractionAngle
        this.getGeocentricEquatorialCoordinates = getGeocentricEquatorialCoordinates
    }
}
