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

package net.arwix.astronomy.math.Degrees

import net.arwix.astronomy.core.DEG_TO_RAD
import net.arwix.astronomy.math.radians.Radian

typealias Degree = Double

fun Degree.normalizeDegree(): Degree {
    if (this < 0.0 && this >= -360.0) return this + 360.0
    if (this >= 360.0 && this < 720) return this - 360.0
    if (this >= 0 && this < 360.0) return this

    var d = this - 360.0 * Math.floor(this / 360.0)
    // Can't use Math.IEEEremainder here because remainder differs
    // from modulus for negative numbers.
    if (d < 0.0) d += 360.0

    return d
}

fun Degree.toRad(): Radian = this * DEG_TO_RAD

