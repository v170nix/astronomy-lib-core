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
@file:JvmName("CalendarExtensions")

package net.arwix.astronomy.core.calendar

import net.arwix.astronomy.core.MJD_J2000
import net.arwix.astronomy.core.PI2
import net.arwix.astronomy.core.SECS_DAY
import java.lang.Math.*
import java.util.*
import java.util.Calendar.MONTH
import java.util.Calendar.YEAR
import java.util.concurrent.TimeUnit

typealias MJD = Double

/**
 * Устанавливает время в календаре используя доли часа
 * Новый объект не создается
 * @param hours    часы в формате час.десятичные доли часа
 * @return Calendar
 */
fun Calendar.setHours(hours: Double) {
    val hour = hours.toInt()
    val minutes = (hours - hour) * 60.0
    val minute = minutes.toInt()
    val seconds = (minutes - minute) * 60.0
    val second = seconds.toInt()
    val millisecond = ((seconds - second) * 1000.0).toInt()

    set(Calendar.HOUR_OF_DAY, hour)
    set(Calendar.MINUTE, minute)
    set(Calendar.SECOND, second)
    set(Calendar.MILLISECOND, millisecond)
}

inline fun Calendar.getHours(): Double {
    return get(Calendar.HOUR_OF_DAY) + (get(Calendar.MINUTE) + get(Calendar.SECOND) / 60.0) / 60.0
}

/**
 * Вычисление модифицированного юлианского дня
 * @return Модифицированная юлианская дата
 */
fun Calendar.getMJD(): MJD {
    val daylight = if (timeZone.inDaylightTime(time)) timeZone.dstSavings else 0
    var y = get(Calendar.YEAR).toLong()
    var m = get(Calendar.MONTH) + 1
    if (m <= 2) {
        m += 12; --y
    }
    val b = if (10000L * y + 100L * m + get(Calendar.DATE).toLong() <= 15821004L)
        -2L + (y + 4716L) / 4 - 1179L else (y / 400 - y / 100 + y / 4)
    val MJDN = 365 * y - 679004L + b + (30.6001 * (m + 1)).toInt() + get(Calendar.DAY_OF_MONTH)
    val MJDF = (abs(get(Calendar.HOUR_OF_DAY)) + abs(get(Calendar.MINUTE)) / 60.0 +
            abs(get(Calendar.SECOND) + get(Calendar.MILLISECOND) / 1000.0) / 3600.0) / 24.0
    return MJDN + MJDF - (get(Calendar.ZONE_OFFSET) + daylight) / 24.0 / 60.0 / 60.0 / 1000.0
}

/**
 * Вычисление юлианского столетия на эпоху J2000
 * @param aMJD юлианская дата
 * @return J2000
 */
fun getJT(aMJD: MJD) = (aMJD - MJD_J2000) / 36525.0

fun Calendar.getJT(applyDeltaT: Boolean = false) = getJT(getMJD() + if (applyDeltaT) getDeltaT(TimeUnit.DAYS) else 0.0)

/**
 * Среднее гринвичское звездное время
 * Greenwich Mean Sidereal Time
 * @param aMJD Время в форме модифицированной юлианской даты
 * @return GMST в радианах
 */
fun getGMST(aMJD: MJD): Double {
    val MJD_0 = Math.floor(aMJD)
    val UT = SECS_DAY * (aMJD - MJD_0) // [сек]
    val T_0 = (MJD_0 - 51544.5) / 36525.0
    val T = (aMJD - 51544.5) / 36525.0
    val gmst = 24110.54841 + 8640184.812866 * T_0 + 1.0027379093 * UT + (0.093104 - 0.0000062 * T) * T * T // [сек]
    return PI2 / SECS_DAY * (gmst % SECS_DAY) // [рад]
}

fun Calendar.getGMST() = getGMST(getMJD())

/**
 * [POLYNOMIAL EXPRESSIONS FOR DELTA T (ΔT)](http://eclipse.gsfc.nasa.gov/LEcat5/deltatpoly.html)
 * @return возвращает разницу во времени между TT и UT в секундах
 */
fun Calendar.getDeltaT(unit: TimeUnit): Double {
    val y = get(YEAR) + (get(MONTH) + 1.0 - 0.5) / 12.0
    val dt: Double

    dt = when (y) {
        in 2005..2050 -> {
            val t = y - 2000.0
            62.92 + (0.32217 + 0.005589 * t) * t
        }
        in -1999..-500 -> -20.0 + 32.0 * (get(YEAR) - 1820.0) / 100.0
        in -500..500 -> {
            val u = y / 100.0
            10583.6 + (-1014.41 + (33.78311 + (-5.952053 + (-0.1798452 + (0.022174192 + 0.0090316521 * u) * u) * u) * u) * u) * u
        }
        in 500..1600 -> {
            val u = (y - 1000.0) / 100.0
            1574.2 + (-556.01 + (71.23472 + (0.319781 + (-0.8503463 + (-0.005050998 + 0.0083572073 * u) * u) * u) * u) * u) * u
        }
        in 1600..1700 -> {
            val u = y - 1600.0
            120.0 + (-0.9808 + (-0.01532 + u / 7129.0) * u) * u
        }
        in 1700..1800 -> {
            val u = y - 1700.0
            8.83 + (0.1603 + (-0.0059285 + (0.00013336 - u / 1174000.0) * u) * u) * u
        }
        in 1800..1860 -> {
            val u = y - 1800.0
            13.72 + (-0.332447 + (0.0068612 + (0.0041116 + (-0.00037436 + (0.0000121272 + (-0.0000001699 + 0.000000000875 * u) * u) * u) * u) * u) * u) * u
        }
        in 1860..1900 -> {
            val u = y - 1860.0
            7.62 + (0.5737 + (-0.251754 + (0.01680668 + (-0.0004473624 + u / 233174.0) * u) * u) * u) * u
        }
        in 1900..1920 -> {
            val u = y - 1900.0
            -2.79 + (1.494119 + (-0.0598939 + (0.0061966 - 0.000197 * u) * u) * u) * u
        }
        in 1920..1941 -> {
            val u = y - 1920.0
            21.20 + (0.84493 + (-0.076100 + 0.0020936 * u) * u) * u
        }
        in 1941..1961 -> {
            val u = y - 1950.0
            29.07 + (0.407 + (-1 / 233.0 + u / 2547.0) * u) * u
        }
        in 1961..1986 -> {
            val u = y - 1975.0
            45.45 + (1.067 + (-1 / 260.0 + -u / 718.0) * u) * u
        }
        in 1086..2005 -> {
            val u = y - 2000.0
            63.86 + (0.3345 + (-0.060374 + (0.0017275 + (0.000651814 + 0.00002373599 * u) * u) * u) * u) * u
        }
        in 2050..2150 -> -20.0 + 32.0 * pow((y - 1820.0) / 100.0, 2.0) - 0.5628 * (2150.0 - y)
        in 2150..3000 -> -20.0 + 32.0 * pow((y - 1820.0) / 100.0, 2.0)
        else -> -20.0 + 32.0 * pow((y - 1820.0) / 100.0, 2.0)
    }
    return when (unit) {
        TimeUnit.NANOSECONDS -> dt * 1000000000.0
        TimeUnit.MICROSECONDS -> dt * 1000000.0
        TimeUnit.MILLISECONDS -> dt * 1000.0
        TimeUnit.SECONDS -> dt
        TimeUnit.MINUTES -> dt / 60.0
        TimeUnit.HOURS -> dt / 3600.0
        TimeUnit.DAYS -> dt / 86400.0
    }
}

inline fun Calendar.resetTime() {
    set(Calendar.HOUR_OF_DAY, 0)
    set(Calendar.MINUTE, 0)
    set(Calendar.SECOND, 0)
    set(Calendar.MILLISECOND, 0)
}

inline fun Calendar.year() = this.get(Calendar.YEAR)
inline fun Calendar.month() = this.get(Calendar.MONTH)
inline fun Calendar.dayOfMonth() = this.get(Calendar.DAY_OF_MONTH)
inline fun Calendar.dayOfYear() = this.get(Calendar.DAY_OF_YEAR)
inline fun Calendar.hourOfDay() = this.get(Calendar.HOUR_OF_DAY)
inline fun Calendar.minute() = this.get(Calendar.MINUTE)
inline fun Calendar.second() = this.get(Calendar.SECOND)
inline fun Calendar.millisecond() = this.get(Calendar.MILLISECOND)
inline fun Calendar.zoneOffset() = this.get(Calendar.ZONE_OFFSET)

inline fun Calendar.year(year: Int) {
    this.set(Calendar.YEAR, year)
}

inline fun Calendar.month(month: Int) = this.set(Calendar.MONTH, month)
inline fun Calendar.dayOfMonth(day: Int) = this.set(Calendar.DAY_OF_MONTH, day)
inline fun Calendar.dayOfYear(day: Int) = this.set(Calendar.DAY_OF_YEAR, day)
inline fun Calendar.hourOfDay(hour: Int) = this.set(Calendar.HOUR_OF_DAY, hour)
inline fun Calendar.minute(minute: Int) = this.set(Calendar.MINUTE, minute)
inline fun Calendar.second(second: Int) = this.set(Calendar.SECOND, second)
inline fun Calendar.millisecond(millisecond: Int) = this.set(Calendar.MILLISECOND, millisecond)

fun fromMJDToCalendar(mjd: MJD, calendar: Calendar, applyDeltaT: Boolean = false) {

    val a = (mjd + 2400001.0).toLong()
    val c = if (a < 2299161) a + 1524L else {
        val b = ((a - 1867216.25) / 36524.25).toLong()
        a + b - (b / 4) + 1525
    }

    val d = ((c - 122.1) / 365.25).toLong()
    val e = 365 * d + d / 4
    val f = ((c - e) / 30.6001).toLong()

    val day = c - e - (30.6001 * f).toLong()
    val month = f - 1 - 12 * (f / 14)
    val year = d - 4715 - ((7 + month) / 10)
    val hours = mjd - floor(mjd) + calendar.get(Calendar.ZONE_OFFSET) / 24.0 / 60.0 / 60.0 / 1000.0
    calendar.set(year.toInt(), month.toInt() - 1, day.toInt())
    calendar.timeZone.useDaylightTime()
    calendar.setHours(hours * 24.0)
    val daylight = if (calendar.timeZone.inDaylightTime(calendar.time)) calendar.timeZone.dstSavings else 0
    if (daylight != 0) calendar.set(Calendar.MILLISECOND, calendar.get(Calendar.MILLISECOND) + daylight)
    if (applyDeltaT) calendar.second(calendar.second() - calendar.getDeltaT(TimeUnit.SECONDS).toInt())
}

fun MJD.toCalendar(applyDeltaT: Boolean = false) = fromMJDToCalendar(this, Calendar.getInstance(), applyDeltaT)
