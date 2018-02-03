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

package net.arwix.astronomy.core.calc

import kotlinx.coroutines.experimental.CommonPool
import kotlinx.coroutines.experimental.async
import net.arwix.astronomy.annotation.Apparent
import net.arwix.astronomy.annotation.Equatorial
import net.arwix.astronomy.annotation.Geocentric
import net.arwix.astronomy.core.calendar.getDeltaT
import net.arwix.astronomy.core.calendar.getMJD
import net.arwix.astronomy.core.calendar.resetTime
import net.arwix.astronomy.core.calendar.setHours
import net.arwix.astronomy.core.coordinates.Location
import net.arwix.astronomy.core.vector.Vector
import net.arwix.astronomy.math.QuadraticInterpolator
import java.lang.Math.cos
import java.lang.Math.sin
import java.util.*
import java.util.concurrent.TimeUnit

class RiseSetCalc(val date: Calendar,
                  val location: Location,
                  @Geocentric @Equatorial @Apparent private val funGetCoordinates: (T: Double) -> Vector,
                  private val sinRefractionAngle: Double) {

    constructor (date: Calendar,
                 location: Location,
                 @Geocentric @Equatorial @Apparent funGetCoordinates: (T: Double) -> Vector,
                 type: ObjectType) : this(date, location, funGetCoordinates, type.sinRefractionAngle)

    sealed class Result {
        data class Rise(val calendar: Calendar) : Result()
        data class Set(val calendar: Calendar) : Result()
        data class RiseSet(val rise: Rise, val set: Set) : Result()
        data class None(val isAbove: Boolean) : Result()
    }

    enum class ObjectType constructor(angle: Double) {
        SUN(-0.833), MOON(+0.133), DOT(-0.5667);

        internal val sinRefractionAngle = Math.sin(Math.toRadians(angle))
    }

    /**
     * Расчет моментов восхода/захода Объекта и наступления сумерек
     * <p/>
     * с учетом рефракции время в UT
     *
     * @return @see RiseSet
     */
    suspend fun calculation(): Result {
        val innerDate = Calendar.getInstance(date.timeZone).apply { this.timeInMillis = date.timeInMillis }
        val deltaT = innerDate.getDeltaT(TimeUnit.DAYS)
        innerDate.resetTime()


        //   return runBlocking<Result>(CommonPool) {
        val MJD0 = innerDate.getMJD()
        val cosLatitude = cos(location.latitude)
        val sinLatitude = sin(location.latitude)

        var hour = 1.0
        var y_0: Double
        var y_plus: Double

        // Инициализация поиска
        var y_minus = getSinAltitude(
                MJD0 + (hour - 1.0) / 24.0,
                deltaT,
                location.longitude,
                cosLatitude, sinLatitude,
                funGetCoordinates) - sinRefractionAngle

        var rise: Result.Rise? = null
        var set: Result.Set? = null
        do {
            val defferedY_0 = async(CommonPool) {
                getSinAltitude(MJD0 + hour / 24.0,
                        deltaT,
                        location.longitude,
                        cosLatitude,
                        sinLatitude,
                        funGetCoordinates) - sinRefractionAngle
            }

            val defferedY_plus = async(CommonPool) {
                getSinAltitude(MJD0 + (hour + 1.0) / 24.0, deltaT,
                        location.longitude,
                        cosLatitude,
                        sinLatitude,
                        funGetCoordinates) - sinRefractionAngle
            }

            y_0 = defferedY_0.await()
            y_plus = defferedY_plus.await()

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
        if (rise != null && set != null) return Result.RiseSet(rise, set)
        if (rise != null) return rise
        if (set != null) return set
        return Result.None(y_minus > 0.0)
    }


}