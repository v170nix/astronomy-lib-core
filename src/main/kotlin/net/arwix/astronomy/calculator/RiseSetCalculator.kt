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
import net.arwix.astronomy.math.QuadraticInterpolator
import java.lang.Math.*
import java.util.*

class RiseSetCalculator : BaseEventCalculator {


    sealed class Result {
        data class Rise(val calendar: Calendar) : Result()
        data class Set(val calendar: Calendar) : Result()
        data class RiseSet(val rise: Rise, val set: Set) : Result()
        data class None(val isAbove: Boolean) : Result()
    }


    enum class ObjectType constructor(angle: Double) {
        SUN(-0.833), MOON(+0.133), DOT(-0.5667);

        internal val sinRefractionAngle = sin(toRadians(angle))
    }


//    private lateinit var innerDate: Calendar
//    private var deltaT: Double = 0.0
//    private var isValid = false
//    private val getGeocentricEquatorialCoordinates: (T: Double, Epoch) -> Vector

    internal constructor (date: Calendar, location: Location, sinRefractionAngle: Double, getGeocentricEquatorialCoordinates: (T: Double, Epoch) -> Vector) :
            super(date, location, sinRefractionAngle, getGeocentricEquatorialCoordinates)

    public constructor (date: Calendar, location: Location, type: ObjectType, getGeocentricEquatorialCoordinates: (T: Double, Epoch) -> Vector) :
            super(date, location, type.sinRefractionAngle, getGeocentricEquatorialCoordinates)

    var result: Result? = null
        get() {
            if (!isValid) result = calls()
            return result
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

    /**
     * Расчет моментов восхода/захода Объекта и наступления сумерек
     * <p/>
     * с учетом рефракции время в UT
     *
     * @return @see RiseSet
     */
    protected fun calls(): Result {
        isValid = false
        // latitude = 65.5;
        // 27.05.2012
//        final double refraction = getSinRefractionAngle(event);
        innerDate.resetTime()
        val MJD0 = innerDate.getMJD()
        val cosLatitude = cos(location.latitude)
        val sinLatitude = sin(location.latitude)

        var hour = 1.0
        var y_0: Double
        var y_plus: Double

        // Инициализация поиска
        var y_minus = getSinAltitude(MJD0 + (hour - 1.0) / 24.0, location.longitude, cosLatitude, sinLatitude) - sinRefractionAngle

        var rise: Result.Rise? = null
        var set: Result.Set? = null

        // перебор интервалов [0h-2h] to [22h-24h]
        do {
            y_0 = getSinAltitude(MJD0 + hour / 24.0, location.longitude, cosLatitude, sinLatitude) - sinRefractionAngle
            y_plus = getSinAltitude(MJD0 + (hour + 1.0) / 24.0, location.longitude, cosLatitude, sinLatitude) - sinRefractionAngle

            // определние параболы по трем значением y_minus,y_0,y_plus
            val quadraticResult = QuadraticInterpolator.getResult(y_minus, y_0, y_plus)

            when (quadraticResult) {
                is QuadraticInterpolator.Result.Root -> {
                    if (y_minus < 0.0)
                        rise = Result.Rise(Calendar.getInstance(date.timeZone).apply { time = date.time; setHours(hour + quadraticResult.root) })
                    else
                        set = Result.Set(Calendar.getInstance(date.timeZone).apply { time = date.time; setHours(hour + quadraticResult.root) })
                }
                is QuadraticInterpolator.Result.Roots -> {
                    val (LT_Rise, LT_Set) = if (quadraticResult.extremum.y < 0.0)
                        Pair(hour + quadraticResult.root2, hour + quadraticResult.root1) else
                        Pair(hour + quadraticResult.root1, hour + quadraticResult.root2)
                    return Result.RiseSet(
                            Result.Rise(Calendar.getInstance(date.timeZone).apply { time = date.time; setHours(LT_Rise) }),
                            Result.Set(Calendar.getInstance(date.timeZone).apply { time = date.time; setHours(LT_Set) }))
                }
            }

            y_minus = y_plus // подготовка к обработке следующего интервала
            hour += 2.0
        } while (!((hour == 25.0) || (rise != null && set != null)))
        isValid = true
        if (rise != null && set != null) return Result.RiseSet(rise, set)
        if (rise != null) return rise
        if (set != null) return set
        return Result.None(y_minus > 0.0)
    }

}