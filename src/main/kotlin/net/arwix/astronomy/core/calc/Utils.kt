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

package net.arwix.astronomy.core.calc

import net.arwix.astronomy.annotation.Apparent
import net.arwix.astronomy.annotation.Equatorial
import net.arwix.astronomy.annotation.Geocentric
import net.arwix.astronomy.core.MJD_J2000
import net.arwix.astronomy.core.calendar.getGMST
import net.arwix.astronomy.core.vector.SphericalVector
import net.arwix.astronomy.core.vector.Vector

/**
 * Синус высоты объекта над горизонтом
 * @param MJD         на расчетную дату
 * @param deltaT      в долях дня
 * @param longitude   долгота в радианах
 * @param cosLatitude косинус широты
 * @param sinLatitude синус широты
 * @return cинус высоты Солнца или Луны в момент искомого события
 */
inline internal fun getSinAltitude(MJD: Double,
                                   deltaT: Double,
                                   longitude: Double,
                                   cosLatitude: Double,
                                   sinLatitude: Double,
                                   @Geocentric @Equatorial @Apparent funGetCoordinates: (T: Double) -> Vector): Double {
    val T = (MJD - MJD_J2000 - deltaT) / 36525.0
    val p: SphericalVector = funGetCoordinates(T).toType()
    // часовой угол
    val tau = getGMST(MJD) + longitude - p.phi
    return sinLatitude * Math.sin(p.theta) + cosLatitude * Math.cos(p.theta) * Math.cos(tau)
}