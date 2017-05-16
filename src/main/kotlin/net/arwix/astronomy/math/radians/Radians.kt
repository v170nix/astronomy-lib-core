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

package net.arwix.astronomy.math.radians

import net.arwix.astronomy.core.PI2
import net.arwix.astronomy.core.PI4
import net.arwix.astronomy.core.RAD_TO_HOUR

typealias Radian = Double

/**
 * Reduce an angle in radians to the range (0 - 2 Pi).
 *
 * @return The reduced radian value.
 */
fun Double.normalize(): Radian {
    if (this >= 0 && this < PI2) return this
    if (this < 0 && this >= -PI2) return this + PI2
    if (this >= PI2 && this < PI4) return this - PI2

    var d = this - PI2 * Math.floor(this / PI2)
    // Can't use Math.IEEE remainder here because remainder differs
    // from modulus for negative numbers.
    if (d < 0.0) d += PI2
    return d
}

/**
 * Radians to hours
 */
fun Radian.toHour(): Double = this * RAD_TO_HOUR