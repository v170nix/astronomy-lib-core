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

package net.arwix.astronomy.physical


enum class PhysicalBody(override val equatorialRadius: Double,
                        override val polarRadius: Double,
                        override val massRatio: Double) : PhysicalProperties {

    SUN(696000.0, 696000.0, 1.0),
    MERCURY(2440.0, 2440.0, 6023600.0),
    VENUS(6051.8, 6051.8, 408523.71),
    EARTH(6378.1366, 6356.7519, 332946.050895), //328900.56
    MARS(3396.19, 3376.2, 3098708.0),
    JUPITER(71492.0, 66854.0, 1047.3486),
    SATURN(60268.0, 54364.0, 3497.898),
    URANUS(25559.0, 24973.0, 22902.98),
    NEPTUNE(24764.0, 24341.0, 19412.24),
    Pluto(1195.0, 1195.0, 1.35E8),
    Moon(1737.4, 1737.4, 2.7068700387534E7),
    EarthMoonBarycenter(0.0, 0.0, 328900.5614)
}