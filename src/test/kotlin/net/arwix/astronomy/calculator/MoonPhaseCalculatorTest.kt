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

import net.arwix.astronomy.core.calendar.*
import net.arwix.astronomy.core.vector.RectangularVector
import net.arwix.astronomy.core.vector.SphericalVector
import net.arwix.astronomy.ephemeris.Obliquity
import net.arwix.astronomy.ephemeris.Precession
import net.arwix.astronomy.math.radians.toDeg
import net.arwix.astronomy.swiss.SwissBody
import net.arwix.astronomy.vsop87.AEarthData
import net.arwix.astronomy.vsop87.VsopData
import org.junit.Test

import java.text.SimpleDateFormat
import java.util.*
import kotlin.collections.ArrayList
import kotlin.system.measureTimeMillis


class MoonPhaseCalculatorTest {
    @Test
    fun calls() {

        val tau_Sun = 8.32 / (1440.0 * 36525.0)

//        val cal = Calendar.getInstance(TimeZone.getTimeZone("PST"))
//        fromMJDToCalendar(57972.77996129096, cal)
//        val format1 = SimpleDateFormat("yyyy-MM-dd HH:mm:ss zzz")
//        println(format1.format(cal.time))

        val t = Calendar.getInstance().apply { setHours(21.0 + 11.0 / 60.0) }.getJT(true)

        //    SwissBody.Earth(Precession.Williams1994(t)).getHeliocentricEclipticCoordinates

        val positionCalculator = PositionCalculator(Precession.Williams1994(t),
                //      SwissBody.Earth(Precession.Williams1994(t)).getHeliocentricEclipticCoordinates
                (AEarthData() as VsopData)::getEclipticCoordinates
        )

        val moonRequest = PositionCalculator.Request.GeocentricEclipticBody(
                SwissBody.Moon(Precession.Williams1994(t)).getGeocentricEclipticCoordinates,
                false
        )

        val sunRequest = PositionCalculator.Request.HeliocentricEclipticBody(
                { _ -> RectangularVector() }
        )

        val oblibity = Obliquity.Williams1994(t)

        val moonEcl = oblibity.rotateFromEquatorialToEcliptic(positionCalculator.getGeocentricEquatorialPositionApparent(t, moonRequest))
        val sunEcl = oblibity.rotateFromEquatorialToEcliptic(positionCalculator.getGeocentricEquatorialPositionApparent(t, sunRequest))

        val delta = moonEcl.toType<SphericalVector>().phi - sunEcl.toType<SphericalVector>().phi

        println("delta position ${delta.toDeg()}")

//Calendar.getInstance().apply { month(12) }
        val cal1 = Calendar.getInstance().apply { year(2000) }
        var result: List<MoonEventsCalculator.MoonEvent> = ArrayList()
        val time = measureTimeMillis {
            result = MoonEventsCalculator(cal1, 1).getNextPhaseEvents(Calendar.getInstance().apply { year(2017); month(12) })
        }

        println("time $time")


        val calendar = Calendar.getInstance()
        val format1 = SimpleDateFormat("yyyy-MM-dd HH:mm:ss zzz")
        result.forEach {
            fromMJDToCalendar(it.mJdET, calendar)
            println("mJdET = ${format1.format(calendar.time)}, ${it.phase}")
            it.eclipse?.let {
                when (it) {
                    is MoonEventsCalculator.SolarEclipse -> {
                        fromMJDToCalendar(it.timeOfMaximumEclipseMJD, calendar, true)
                        println("solar time ${format1.format(calendar.time)}")
                    }
                    is MoonEventsCalculator.MoonEclipse -> {
                        fromMJDToCalendar(it.timeOfMaximumEclipseMJD, calendar, true)
                        println("moon time ${format1.format(calendar.time)}")
                    }
                }
            }
        }


//179.9121105410326
//179.91764596954215


        val d = SwissBody.Moon(
                Precession.Williams1994(t))
                .getGeocentricEclipticCoordinates(t)
                .toType<SphericalVector>().phi -
                (RectangularVector() -
                        //       (CEarthData() as VsopData).getEclipticCoordinates(t )).toType<SphericalVector>().phi
                        SwissBody.Earth(Precession.Williams1994(t)).getHeliocentricEclipticCoordinates(t)).toType<SphericalVector>().phi

        println("d ${d.toDeg()}")


        MoonPhaseCalculator(
                SwissBody.Moon(Precession.Williams1994(t)).getGeocentricEclipticCoordinates,
                { time ->
                    (RectangularVector() -
                            //       (CEarthData() as VsopData).getEclipticCoordinates(t )).toType<SphericalVector>().phi
                            SwissBody.Earth(Precession.Williams1994(t)).getHeliocentricEclipticCoordinates(time - tau_Sun))
//                    return@MoonPhaseCalculator (RectangularVector() -
//                            (AEarthData() as VsopData).getEclipticCoordinates(time))
                }
        ).calls(2017)
//

    }

}