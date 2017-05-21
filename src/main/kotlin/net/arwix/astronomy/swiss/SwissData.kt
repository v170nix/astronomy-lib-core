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

package net.arwix.astronomy.swiss

import net.arwix.astronomy.core.coordinates.FunGetHeliocentricEclipticCoordinates
import net.arwix.astronomy.core.vector.RectangularVector


abstract class SwissData {

    abstract protected val tabl: DoubleArray
    abstract protected val tabb: DoubleArray
    abstract protected val tabr: DoubleArray
    abstract protected val args: IntArray
    abstract protected val maxargs: Int
    abstract protected val max_harmonic: IntArray
    abstract protected val max_power_of_t: Int
    abstract protected val distance: Double

    val timescale = 3652500.0
    val trunclvl = 1.0

    open val getCoordinates: FunGetHeliocentricEclipticCoordinates = { t ->
        RectangularVector(
                gplan(t, args, distance, tabb, tabl,
                        tabr, max_harmonic, max_power_of_t,
                        maxargs, timescale, trunclvl)
        )
    }

}



