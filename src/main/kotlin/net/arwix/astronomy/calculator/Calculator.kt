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
import net.arwix.astronomy.core.coordinates.Location
import net.arwix.astronomy.core.vector.SphericalVector
import net.arwix.astronomy.core.vector.Vector
import net.arwix.astronomy.core.vector.VectorType
import java.util.*


abstract class Calculator<T> {

    protected var deltaT: Double = 0.0
    protected var lastResult: T? = null

    var date: Calendar
        set(value) {
            isValidResult = false; field = value
        }

    var location: Location
        set(value) {
            isValidResult = false; field = value
        }

    protected var getCoordinates: (T: Double, Epoch) -> Vector
    protected var isValidResult = false


    constructor(date: Calendar, location: Location, funGeocentricEquatorialCoordinates: (T: Double, Epoch) -> Vector) {
        this.date = date
        this.location = location
        this.getCoordinates = funGeocentricEquatorialCoordinates;
    }

    /**
     * Синус высоты объекта над горизонтом
     * @param MJD         на расчетную дату
     * @param longitude   долгота в радианах
     * @param cosLatitude косинус широты
     * @param sinLatitude синус широты
     * @return cинус высоты Солнца или Луны в момент искомого события
     */
    protected fun getSinAltitude(MJD: Double, longitude: Double, cosLatitude: Double, sinLatitude: Double): Double {
        val T = (MJD - MJD_J2000 - deltaT) / 36525.0
        val p = getCoordinates(T, Epoch.APPARENT).getVectorOfType(VectorType.SPHERICAL) as SphericalVector
        // часовой угол
        val tau = getGMST(MJD) + longitude - p.phi
        return sinLatitude * Math.sin(p.theta) + cosLatitude * Math.cos(p.theta) * Math.cos(tau)
    }

    fun getResult(): T {
        if (!isValidResult || lastResult == null) {
            lastResult = calls()
            isValidResult = true
        }
        return lastResult!!
    }

    protected abstract fun calls(): T


//    companion object {
//
//        inline fun <reified K, F> create(date: Calendar, location: Location, noinline getGeocentricEquatorialCoordinates: (T: Double, Epoch) -> Vector): Calculator<F> {
//            when (K::class) {
//                RiseSetCalculator::class -> RiseSetCalculator(date, location, getGeocentricEquatorialCoordinates)
//                else -> throw IllegalAccessException()
//            }
//
//        }
//    }

}



