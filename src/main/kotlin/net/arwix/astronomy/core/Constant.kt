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
@file:JvmName("Constants")
package net.arwix.astronomy.core

val PI2 = 2.0 * Math.PI
val RAD = Math.PI / 180.0
val DEG = 180.0 / Math.PI
val ARCS = 3600.0 * 180.0 / Math.PI

// радиусы Земли, Солнца, Луны в км
val R_Earth = 6378.137
val R_Sun = 696000.0
val R_Moon = 1738.0

val MJD_J2000 = 51544.5        // MJD на эпоху J2000.0
val T_J2000 = 0.0           // эпоха J2000.0
val T_B1950 = -0.500002108   // эпоха B1950
val JD_SECOND = 0.000011574074074074074074
val JD_MINUTE = 0.00069444444444444444444
val JD_HOUR = 0.041666666666666666666
val JD_DAY = 1.0
val SECS_DAY = 86400.0 // колличество секунд в сутках

val kGauss = 0.01720209895  // гравитационная константа
val GM_Sun = kGauss * kGauss  // [AU^3/d^2]

val AU = 149597870.0    // 1ае

val C_Light = 173.14         // скорость света [AU/d]