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


abstract class SwissData {

    abstract internal val tabl: DoubleArray
    abstract internal val tabb: DoubleArray
    abstract internal val tabr: DoubleArray
    abstract internal val args: IntArray
    abstract internal val maxargs: Int
    abstract internal val max_harmonic: IntArray
    abstract internal val max_power_of_t: Int
    abstract internal val distance: Double
    internal val timescale = 3652500.0
    internal val trunclvl = 1.0
}



