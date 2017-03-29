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
import net.arwix.astronomy.core.Position
import net.arwix.astronomy.core.coordinates.EclipticCoordinates
import net.arwix.astronomy.core.coordinates.Location
import net.arwix.astronomy.core.vector.RectangularVector
import net.arwix.astronomy.core.vector.Vector
import net.arwix.astronomy.vsop87.CEarthData
import org.junit.Test
import java.util.*


class RiseSetCalculatorTest {
    @Test
    fun calls() {
        val date = Calendar.getInstance(TimeZone.getTimeZone("GMT+4"));
        date.set(Calendar.YEAR, 2014);
        date.set(Calendar.MONTH, 8);
        date.set(Calendar.DAY_OF_MONTH, 14);
        date.set(Calendar.HOUR_OF_DAY, 20);
        val location = Location(Math.toRadians(30.3290233), Math.toRadians(59.909328));

        val c = RiseSetCalculator(date, location, { t: Double, e: Epoch ->

            Position(CEarthData()).getGeocentricEquatorialPosition(t, object : EclipticCoordinates<Any> {
                override fun getEpoch(): Epoch {
                    return e
                }

                override fun getIdObject(): Any {
                    return "Sun"
                }

                override fun getEclipticCoordinates(T: Double): Vector {
                    return RectangularVector(0.0, 0.0, 0.0)
                }

            })

        }, RiseSetCalculator.ObjectType.SUN)

        c.getResult().let {
            when (it) {
                is RiseSetCalculator.Result.RiseSet -> {
                    printResult(it.set.calendar.timeInMillis)
                }
            }
        }

    }

    fun printResult(tMill: Long) {
        val calendar = Calendar.getInstance()

        System.out.println("EVENTS - " +
                calendar.apply { this.timeInMillis = tMill }.time.toString() + "; ")

    }
}
