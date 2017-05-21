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

package net.arwix.astronomy.math

fun DoubleArray.polynomialSum(x: Double): Double {
    var t = 1.0
    return this.fold(0.0, { acc, d -> (acc + d * t).let { t *= x; it } })
}

/**
 * Module operation in arcseconds.
 *
 * @param x Value in arcseconds.
 * @return module.
 */
fun Double.mod3600() = this - 1296000.0 * Math.floor(this / 1296000.0)