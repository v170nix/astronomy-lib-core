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

@JvmField val PI2 = 2.0 * Math.PI
@JvmField val PI4 = 4.0 * Math.PI
@Deprecated("use DEG_TO_RAD")
@JvmField val RAD = Math.PI / 180.0
@Deprecated("use RAD_TO_DEG")
@JvmField val DEG = 180.0 / Math.PI

@Deprecated("use RAD_TO_ARCSEC")
@JvmField val ARCS = 3600.0 * 180.0 / Math.PI

@JvmField val ARCSEC_2RAD = PI2 / (360.0 * 3600.0)

/** Arc minutes in one degree = 60.  */
@JvmField val MINUTES_PER_DEGREE = 60.0
/** Arc seconds in one degree = 3600.  */
@JvmField val SECONDS_PER_DEGREE = 60.0 * MINUTES_PER_DEGREE
/** Arc seconds to radians.  */
@JvmField val ARCSEC_TO_RAD = Math.PI / (180.0 * 3600.0)

/** Radians to arc seconds.  */
@JvmField val RAD_TO_ARCSEC = 1.0 / ARCSEC_TO_RAD
/** Arc seconds to degrees.  */
@JvmField val ARCSEC_TO_DEG = 1.0 / 3600.0
/** Radians to hours.  */
@JvmField val RAD_TO_HOUR = 180.0 / (15.0 * Math.PI)
/** Radians to days.  */
@JvmField val RAD_TO_DAY = RAD_TO_HOUR / 24.0
/** Radians to degrees.  */
@JvmField val RAD_TO_DEG = 180.0 / Math.PI
/** Degrees to radians.  */
@JvmField val DEG_TO_RAD = 1.0 / RAD_TO_DEG

// радиусы Земли, Солнца, Луны в км
// @JvmField val R_Earth = 6378.137
// @JvmField val R_Sun = 696000.0
// @JvmField val R_Moon = 1738.0

@JvmField val MJD_J2000 = 51544.5        // MJD на эпоху J2000.0
@JvmField val T_J2000 = 0.0           // эпоха J2000.0
@JvmField val T_B1950 = -0.500002108   // эпоха B1950
@JvmField val JD_SECOND = 0.000011574074074074074074
@JvmField val JD_MINUTE = 0.00069444444444444444444
@JvmField val JD_HOUR = 0.041666666666666666666
@JvmField val JD_DAY = 1.0

/** Julian century conversion constant = 100 * days per year.  */
val JULIAN_DAYS_PER_CENTURY = 36525.0
/** Our default epoch.
 * The Julian Day which represents noon on 2000-01-01.  */
@JvmField val J2000 = 2451545.0
/** Length of a tropical year in days for B1950.  */
@JvmField val TROPICAL_YEAR = 365.242198781


@JvmField val SECS_DAY = 86400.0 // колличество секунд в сутках

@JvmField val kGauss = 0.01720209895  // гравитационная константа
@JvmField val GM_Sun = kGauss * kGauss  // [AU^3/d^2]

@JvmField val AU = 149597870.7    // 1ае

@JvmField val C_Light = 173.14         // скорость света [AU/d]