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


enum class SwissObject(
        val getHeliocentricEclipticCoordinates: FunGetHeliocentricEclipticCoordinates?) {

    Libration(LibrationObject.getCoordinates),
    MERCURY(MercuryObject.getCoordinates),
    VENUS(VenusObject.getCoordinates),
    Earth_Moon_Barycenter(EarthMoonBarycenterObject.getCoordinates),
    EARTH(null),
    Moon(null),
    MARS(MarsObject.getCoordinates),
    JUPITER(JupiterObject.getCoordinates),
    SATURN(SaturnObject.getCoordinates),
    URANUS(UranusObject.getCoordinates),
    NEPTUNE(NeptuneObject.getCoordinates),
    Pluto(PlutoObject.getCoordinates)

}
