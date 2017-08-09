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
import net.arwix.astronomy.annotation.Geocentric
import net.arwix.astronomy.core.JULIAN_DAYS_PER_CENTURY
import net.arwix.astronomy.core.MJD_J2000
import net.arwix.astronomy.core.PI2
import net.arwix.astronomy.core.calendar.fromMJDToCalendar
import net.arwix.astronomy.core.calendar.getJT
import net.arwix.astronomy.core.calendar.resetTime
import net.arwix.astronomy.core.modulo
import net.arwix.astronomy.core.vector.SphericalVector
import net.arwix.astronomy.core.vector.Vector
import net.arwix.astronomy.math.PegasusInterpolator
import java.text.SimpleDateFormat
import java.util.*

@Deprecated("use MoonEventsCalculator")
class MoonPhaseCalculator(
        @Geocentric @Ecliptic val funGetMoonCoordinates: (t: Double) -> Vector,
        @Geocentric @Ecliptic val funGetSunCoordinates: (t: Double) -> Vector) {

    var phase = 0.0

    private fun phaseFunction(t: Double): Double {
        val mjd = t * JULIAN_DAYS_PER_CENTURY + MJD_J2000
        val format1 = SimpleDateFormat("yyyy-MM-dd hh:mm:ss")
        val calendar = Calendar.getInstance(TimeZone.getTimeZone("UTC"))
        fromMJDToCalendar(mjd, calendar)

        val delta = funGetMoonCoordinates(t).toType<SphericalVector>().phi -
                funGetSunCoordinates(t).toType<SphericalVector>().phi

        val result = (funGetMoonCoordinates(t).toType<SphericalVector>().phi -
                funGetSunCoordinates(t).toType<SphericalVector>().phi
                - phase * Math.PI / 2.0 + Math.PI).modulo(PI2) - Math.PI

        val stringTime = format1.format(calendar.time)


//        println("$phase $delta $result $stringTime")

        return result
    }

    companion object {
        val dT = 7.0 / JULIAN_DAYS_PER_CENTURY // Step (1 week)
        val accuracy = (0.5 / 1440.0) / JULIAN_DAYS_PER_CENTURY // Desired Accuracy (0.5 min)
    }

    fun calls(year: Int = 2017) {
        var t0 = Calendar.getInstance()
                .apply {
                    set(year - 1, Calendar.DECEMBER, 1)
                    resetTime()
                }.getJT(false)

        var t1 = t0 + dT
        var tPhase = t0

        (0..13).forEach {

            (0..3).forEach { iPhase ->

                phase = iPhase.toDouble()

                var d0 = phaseFunction(t0)
                var d1 = phaseFunction(t1)

                while (d0 * d1 > 0.0 || d1 < d0) {
                    t0 = t1
                    d0 = d1
                    t1 += dT
                    d1 = phaseFunction(t1)
                }

                PegasusInterpolator.getResult(this::phaseFunction, t0, t1, accuracy).let {
                    if (it is PegasusInterpolator.Result.Root) {
                        tPhase = it.root

                        val mjd = tPhase * JULIAN_DAYS_PER_CENTURY + MJD_J2000

                        val MjdRound = Math.floor(86400.0 * mjd + 0.5) / 86400.0 + 0.000001;

                        val format1 = SimpleDateFormat("yyyy-MM-dd HH:mm:ss")
                        val calendar = Calendar.getInstance(TimeZone.getTimeZone("UTC"))
                        fromMJDToCalendar(MjdRound, calendar)

                        println("$iPhase ${format1.format(calendar.time)}")
                    }
                }

                t0 = tPhase
                t1 = t0 + dT
            }
        }
    }

}