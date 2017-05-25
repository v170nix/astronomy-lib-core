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

package net.arwix.astronomy.core.coordinates

import net.arwix.astronomy.annotation.Azimuth
import net.arwix.astronomy.annotation.Equatorial
import net.arwix.astronomy.annotation.Geocentric
import net.arwix.astronomy.core.calendar.getGMST
import net.arwix.astronomy.core.vector.Matrix
import net.arwix.astronomy.core.vector.SphericalVector
import net.arwix.astronomy.core.vector.Vector
import java.util.*


data class AzimuthCoordintates(val lambda: Double, val phi: Double) {

    private val hourMatrix = Matrix(Matrix.Axis.Y, Math.PI / 2.0 - phi)

    /**
     * часовой угол
     * @param aMJD модифицированная юлианская дата
     * @return часовой угол в радианах
     */
    fun getHourAngle(aMJD: Double): Double {
        return getGMST(aMJD) + lambda
    }

    fun getHourAngle(calendar: Calendar): Double {
        return calendar.getGMST() + lambda
    }

    /**
     * Получение азимутальных координат
     * @param GMST Среднее гринвичское звездное время
     * @param inVector ГеоЦентрические экваториальные координаты
     * @return азимутальный вектор
     */
    @Azimuth @Geocentric
    fun getCoordinates(GMST: Double, @Equatorial @Geocentric inVector: Vector): Vector {
        val vector: SphericalVector = inVector.toType() // getVectorOfType(VectorType.SPHERICAL) as SphericalVector
        val tau = GMST + lambda - vector.phi
        return hourMatrix * SphericalVector(tau, vector.theta, vector.r)
    }

    fun getCoordinates(calendar: Calendar, geoEquatorialVector: Vector): Vector = getCoordinates(calendar.getGMST(), geoEquatorialVector)


}