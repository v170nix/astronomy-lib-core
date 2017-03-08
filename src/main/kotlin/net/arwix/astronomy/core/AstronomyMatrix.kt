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

import net.arwix.astronomy.core.vector.Matrix

object AstronomyMatrix {

    enum class Coordinates {ECLIPTIC, EQUATORIAL }

    /**
     * Преобразование экваториальных координат в эклиптические и обратно
     * Transformation of equatorial to ecliptical coordinates
     *
     * @param T Время в юлианских столетиях с эпохи J2000
     * @return Матрица преобразования
     */
    fun createTransformationCoordinates(T: Double, fromType: Coordinates, toType: Coordinates): Matrix {
        if (fromType == toType) throw IllegalArgumentException()
        return when (fromType) {
            Coordinates.EQUATORIAL -> {
                // Преобразование экваториальных координат в эклиптические
                // Transformation of equatorial to ecliptical coordinates
                // https://en.wikipedia.org/wiki/Axial_tilt
                val eps = Math.toRadians(23.43929111 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T / 3600.0)
                Matrix(Matrix.Axis.X, eps)
            }
            Coordinates.ECLIPTIC -> {
                Matrix.transpose(createTransformationCoordinates(T, Coordinates.EQUATORIAL, Coordinates.ECLIPTIC))
            }
        }
    }

    /**
     * Прецессия в эклиптических координатах Precession of ecliptic coordinates
     * Прецессия в экваториальных координатах Precession of equatorial coordinates
     * @param T1 Исходная эпоха в J2000
     * @param T2 Требуемая эпоха в J2000
     * @return Матрица преобразования координат
     */
    fun createPrecession(T1: Double, T2: Double, type: Coordinates) = when (type) {
        Coordinates.ECLIPTIC -> {
            val dT = T2 - T1
            val Pi: Double
            val pi: Double
            val p_a: Double
            Pi = 174.876383889 * RAD + ((3289.4789 + 0.60622 * T1) * T1 + (-869.8089 - 0.50491 * T1 + 0.03536 * dT) * dT) / ARCS
            pi = (47.0029 - (0.06603 - 0.000598 * T1) * T1 + (-0.03302 + 0.000598 * T1 + 0.000060 * dT) * dT) * dT / ARCS
            p_a = (5029.0966 + (2.22226 - 0.000042 * T1) * T1 + (1.11113 - 0.000042 * T1 - 0.000006 * dT) * dT) * dT / ARCS
            Matrix(Matrix.Axis.Z, -(Pi + p_a)) * Matrix(Matrix.Axis.X, pi) * Matrix(Matrix.Axis.Z, Pi)
        }
        Coordinates.EQUATORIAL -> {
            val dT = T2 - T1
            val zeta: Double
            val z: Double
            val theta: Double
            zeta = (2306.2181 + (1.39656 - 0.000139 * T1) * T1 + (0.30188 - 0.000344 * T1 + 0.017998 * dT) * dT) * dT / ARCS
            z = zeta + (0.79280 + 0.000411 * T1 + 0.000205 * dT) * dT * dT / ARCS
            theta = (2004.3109 - (0.85330 + 0.000217 * T1) * T1 - (0.42665 + 0.000217 * T1 + 0.041833 * dT) * dT) * dT / ARCS
            Matrix(Matrix.Axis.Z, -z) * Matrix(Matrix.Axis.Y, theta) * Matrix(Matrix.Axis.Z, -zeta)
        }
    }

    /**
     * Преобразование от средних к истинным экватору и равноденствию
     * Transformation from mean to true equator and equinox
     * @param T Время в юлианских столетиях от эпохи J2000
     * @return матрица нутации Nutation matrix
     */
    fun createNutation(T: Double): Matrix {
        val ls: Double = PI2 * frac(0.993133 + 99.997306 * T)
        val D: Double = PI2 * frac(0.827362 + 1236.853087 * T)
        val F: Double = PI2 * frac(0.259089 + 1342.227826 * T)
        val N: Double = PI2 * frac(0.347346 - 5.372447 * T)
        val eps: Double = 0.4090928 - 2.2696E-4 * T
        val dpsi = (-17.200 * Math.sin(N) - 1.319 * Math.sin(2 * (F - D + N)) - 0.227 * Math.sin(2 * (F + N)) + 0.206 * Math.sin(2 * N) + 0.143 * Math.sin(ls)) / ARCS
        val deps = (+9.203 * Math.cos(N) + 0.574 * Math.cos(2 * (F - D + N)) + 0.098 * Math.cos(2 * (F + N)) - 0.090 * Math.cos(2 * N)) / ARCS
        return Matrix(Matrix.Axis.X, -eps - deps) * Matrix(Matrix.Axis.Z, -dpsi) * Matrix(Matrix.Axis.X, eps)
    }


}