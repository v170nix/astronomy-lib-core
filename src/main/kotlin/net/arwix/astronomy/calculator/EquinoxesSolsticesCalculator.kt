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

import kotlinx.coroutines.experimental.CommonPool
import kotlinx.coroutines.experimental.async
import kotlinx.coroutines.experimental.runBlocking
import net.arwix.astronomy.core.PI_OVER_TWO
import net.arwix.astronomy.core.calendar.*
import net.arwix.astronomy.core.vector.RectangularVector
import net.arwix.astronomy.core.vector.SphericalVector
import net.arwix.astronomy.ephemeris.Precession
import net.arwix.astronomy.swiss.SwissBody
import java.util.*

class EquinoxesSolsticesCalculator {

    companion object {
        private val precision = 0.1 / 24.0 / 3600.0
    }

    data class Result(val request: Request, val calendar: Calendar)

    sealed class Request(val year: Int) : Iterator<Request> {

        companion object {

            fun getAll(year: Int) = listOf(SpringEquinox(year), SummerSolstice(year), AutumnEquinox(year), WinterSolstice(year))

            internal fun getDelta(request: Request, longitude: Double) = when (request) {
                is EquinoxesSolsticesCalculator.Request.SpringEquinox -> 58.13 * Math.sin(-longitude)
                is EquinoxesSolsticesCalculator.Request.SummerSolstice -> 58.13 * Math.sin(PI_OVER_TWO - longitude)
                is EquinoxesSolsticesCalculator.Request.AutumnEquinox -> 58.13 * Math.sin(Math.PI - longitude)
                is EquinoxesSolsticesCalculator.Request.WinterSolstice -> 58.13 * Math.sin(-PI_OVER_TWO - longitude)
            }

            internal fun getInitMonth(request: Request) = when (request) {
                is EquinoxesSolsticesCalculator.Request.SpringEquinox -> 2
                is EquinoxesSolsticesCalculator.Request.SummerSolstice -> 5
                is EquinoxesSolsticesCalculator.Request.AutumnEquinox -> 8
                is EquinoxesSolsticesCalculator.Request.WinterSolstice -> 11
            }
        }

        class SpringEquinox(year: Int) : Request(year) {
            override fun hasNext() = true
            override fun next() = SummerSolstice(year)
        }

        class SummerSolstice(year: Int) : Request(year) {
            override fun hasNext() = true
            override fun next() = AutumnEquinox(year)
        }

        class AutumnEquinox(year: Int) : Request(year) {
            override fun hasNext() = true
            override fun next() = WinterSolstice(year)
        }

        class WinterSolstice(year: Int) : Request(year) {
            override fun hasNext() = false
            override fun next() = SpringEquinox(year)
        }
    }

    private val positionRequest = PositionCalculator.Request.HeliocentricEclipticBody(
            { _ -> RectangularVector() })

    private fun getPositionCalculator(t: Double,
                                      funCreatePositionCalculator: ((t: Double) -> PositionCalculator)? = null): PositionCalculator {
        return funCreatePositionCalculator?.invoke(t) ?:
                Precession.Williams1994(t).let {
                    PositionCalculator(it, SwissBody.Earth(it).getHeliocentricEclipticCoordinates)
                }
    }

    fun getResult(request: Request, funCreatePositionCalculator: ((t: Double) -> PositionCalculator)? = null): Result {

        val calendar = Calendar.getInstance()
        var mjd = calendar
                .apply { year(request.year); resetTime(); month(Request.getInitMonth(request)); dayOfMonth(1) }
                .getMJD()
        var isFirst = true
        var innerPositionCalculator = getPositionCalculator(getJT(mjd), funCreatePositionCalculator)
        var delta: Double = 100.0
        var obliquity = innerPositionCalculator.precession.getNearestObliquityModel(getJT(mjd))
        do {
            if (delta >= 1 && !isFirst) {
                innerPositionCalculator = getPositionCalculator(getJT(mjd), funCreatePositionCalculator)
                obliquity = innerPositionCalculator.precession.getNearestObliquityModel(getJT(mjd))
            }
            val d = innerPositionCalculator.getGeocentricEquatorialPositionApparent(getJT(mjd), positionRequest)
            val p = obliquity.rotateFromEquatorialToEcliptic(d).toType<SphericalVector>().phi

            delta = Request.getDelta(request, p)
            mjd += delta
            isFirst = false
        } while (Math.abs(delta) > precision)
        fromMJDToCalendar(mjd, calendar, true)
        return Result(request, calendar)
    }

    fun getResult(year: Int, funCreatePositionCalculator: ((t: Double) -> PositionCalculator)? = null): List<Result> {
        return runBlocking(CommonPool) {
            Request.getAll(year).let { requests ->
                Array(4) { async(CommonPool) { getResult(requests[it], funCreatePositionCalculator) } }
            }.map { it.await() }
        }
    }

}