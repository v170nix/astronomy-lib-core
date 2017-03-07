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

package net.arwix.astronomy.core.calendar

import net.arwix.astronomy.core.*
import java.lang.Math.abs
import java.util.Calendar.*
import java.util.Calendar.ZONE_OFFSET
import java.util.Calendar.HOUR_OF_DAY
import java.util.Calendar





object CalendarMath {
    /**
     * Устанавливает время в календаре используя доли часа
     * Новый объект не создается
     * @param calendar календарь в котором будет установлено время
     * @param hours    часы в формате час.десятичные доли часа
     * @return Calendar
     */
    fun setHours(calendar: Calendar, hours: Double): Calendar {
        val hour = hours.toInt()
        val minutes = (hours - hour) * 60.0
        val minute = minutes.toInt()
        val seconds = (minutes - minute) * 60.0
        val second = seconds.toInt()
        val millisecond = ((seconds - second) * 1000.0).toInt()

        calendar.set(HOUR_OF_DAY, hour)
        calendar.set(MINUTE, minute)
        calendar.set(SECOND, second)
        calendar.set(MILLISECOND, millisecond)
        return calendar
    }

    /**
     * Вычисление модифицированного юлианского дня с учетом временной зоны календаря
     * ** Важно: все числа даты и времени должны быть положительные **
     * @param calendar календарь
     * @return Модифицированная юлианская дата
     */
    fun getMJD(calendar: Calendar): Double {
        val b: Long
        var y = calendar[YEAR].toLong()
        var m = calendar[MONTH] + 1
        if (m <= 2) { m += 12; --y }
        b = if (10000L * y + 100L * m + calendar[DATE].toLong() <= 15821004L)
            -2L + (y + 4716L) / 4 - 1179L else (y / 400 - y / 100 + y / 4)

        val MJDN = 365 * y - 679004L + b + (30.6001 * (m + 1)).toInt() + calendar[DATE]
        val MJDF = (abs(calendar[HOUR_OF_DAY]) + Math.abs(calendar[MINUTE]) / 60.0 +
                abs(calendar[SECOND] + calendar[MILLISECOND] / 1000.0) / 3600.0) / 24.0
        return MJDN + MJDF - calendar[ZONE_OFFSET] / 24.0 / 60.0 / 60.0 / 1000.0
    }

    /**
     * Вычисление юлианского столетия на эпоху J2000
     * @param aMJD юлианская дата
     * @return J2000
     */
    fun getJT(aMJD: Double): Double {
        return (aMJD - MJD_J2000) / 36525.0
    }

    /**
     * Вычисление юлианского столетия на эпоху J2000
     * @param calendar календарь может быть в локальном времени
     * @return J2000
     */
    fun getJT(calendar: Calendar): Double {
        return getJT(getMJD(calendar))
    }

    /**
     * [POLYNOMIAL EXPRESSIONS FOR DELTA T (ΔT)](http://eclipse.gsfc.nasa.gov/LEcat5/deltatpoly.html)
     * @return возвращает разницу во времени между TT и UT в секундах
     */
    fun getDeltaTofSecond(calendar: Calendar): Double {
        val y = calendar[YEAR] + (calendar[MONTH] + 1.0 - 0.5) / 12.0
        val u: Double
        val dt: Double
        if (y < 1986 || y > 2050) throw IndexOutOfBoundsException()
        if (y < 2005) {
            u = y - 2000.0
            dt = 63.86 + (0.3345 + (-0.060374 + (0.0017275 + (0.000651814 + 0.00002373599 * u) * u) * u) * u) * u
        } else {
            u = y - 2000.0
            dt = 62.92 + (0.32217 + 0.005589 * u) * u
        }
        return dt
    }

    /**
     * [POLYNOMIAL EXPRESSIONS FOR DELTA T (ΔT)](http://eclipse.gsfc.nasa.gov/LEcat5/deltatpoly.html)
     * @return возвращает разницу во времени между ET и UT в долях суток
     */
    fun getDeltaTofDay(calendar: Calendar): Double {
        return getDeltaTofSecond(calendar) / 60.0 / 60.0 / 24.0
    }

    /**
     * Среднее гринвичское звездное время
     * Greenwich Mean Sidereal Time
     * @param aMJD Время в форме модифицированной юлианской даты
     * @return GMST в радианах
     */
    fun getGMST(aMJD: Double): Double {
        val MJD_0 = Math.floor(aMJD)
        val UT = SECS_DAY * (aMJD - MJD_0) // [сек]
        val T_0 = (MJD_0 - 51544.5) / 36525.0
        val T = (aMJD - 51544.5) / 36525.0
        val gmst = 24110.54841 + 8640184.812866 * T_0 + 1.0027379093 * UT + (0.093104 - 0.0000062 * T) * T * T // [сек]
        return PI2 / SECS_DAY * (gmst % SECS_DAY) // [рад]
    }


}