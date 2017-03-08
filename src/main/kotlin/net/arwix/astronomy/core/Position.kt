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

package net.arwix.astronomy.core

import net.arwix.astronomy.core.AstronomyMatrix.Coordinates.ECLIPTIC
import net.arwix.astronomy.core.AstronomyMatrix.Coordinates.EQUATORIAL
import net.arwix.astronomy.core.coordinates.EclipticCoordinates
import net.arwix.astronomy.core.vector.Vector
import java.lang.IllegalArgumentException


class Position(val earthEclipticCoordinates: EclipticCoordinates<*>) {

    var cacheListener: CacheListener? = null

    interface CacheListener {
        fun getEclipticCacheData(id: String, T: Double): Vector?
        fun putEclipticCacheData(id: String, T: Double, vector: Vector)
    }

    /**
     * get geometric heliocentric ecliptic coordinates (minus 1-way light-time)
     * @param T Юлианские столетия (ET) Time in Julian centuries since J2000
     * @return Vector
     */
    fun getHeliocentricEclipticPosition(T: Double, objectEclipticCoordinates: EclipticCoordinates<*>): Vector {
        val earthEcliptic = getData(earthEclipticCoordinates, T)
        val objectEcliptic = getData(objectEclipticCoordinates, T)
        val objGeoEcliptic = objectEcliptic - earthEcliptic
        val dT = objGeoEcliptic.normalize() / C_Light / 36525.0
        return objectEclipticCoordinates.getEclipticCoordinates(T - dT)
    }

    /**
     * Вычисляет геоцентрические (в центре Земля) экваториальные координаты
     * @param T     Юлианские столетия (ET) Time in Julian centuries since J2000
     * @return Vector
     */
    fun getGeocentricEquatorialPosition(T: Double, objectEclipticCoordinates: EclipticCoordinates<*>): Vector {
        if (earthEclipticCoordinates.getEpoch() != objectEclipticCoordinates.getEpoch()) throw IllegalArgumentException("Epoch don't equal")

        var earthEcliptic = getData(earthEclipticCoordinates, T)
        var objectEcliptic = getData(objectEclipticCoordinates, T)
        var objGeoEcliptic = objectEcliptic - earthEcliptic
        val dT = objGeoEcliptic.normalize() / C_Light / 36525.0
        val innerT = T - dT
        when (earthEclipticCoordinates.getEpoch()) {
            Epoch.APPARENT -> {
                earthEcliptic = getData(earthEclipticCoordinates, innerT)
                objectEcliptic = getData(objectEclipticCoordinates, innerT)
                objGeoEcliptic = objectEcliptic - earthEcliptic
                return AstronomyMatrix.createNutation(innerT) * AstronomyMatrix.createTransformationCoordinates(innerT, ECLIPTIC, EQUATORIAL) * objGeoEcliptic

            }
            Epoch.J2000 -> {
                objectEcliptic = objectEclipticCoordinates.getEclipticCoordinates(innerT)
                objGeoEcliptic = objectEcliptic - earthEcliptic
                return AstronomyMatrix.createTransformationCoordinates(T_J2000, ECLIPTIC, EQUATORIAL) * objGeoEcliptic
            }
        }

    }

    private fun getData(coordinates: EclipticCoordinates<*>, T: Double): Vector {
        return cacheListener?.let { listener ->
            listener.getEclipticCacheData(coordinates.getIdObject().toString(), T)
                    ?.let { return it }
                    ?: coordinates.getEclipticCoordinates(T)
                    .let { listener.putEclipticCacheData(coordinates.getIdObject().toString(), T, it); return it }
        } ?: coordinates.getEclipticCoordinates(T)
    }

}
