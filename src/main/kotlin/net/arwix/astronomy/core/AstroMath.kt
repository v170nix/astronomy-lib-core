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
import net.arwix.astronomy.core.vector.Matrix.Axis.*
import java.lang.Math.cos
import java.lang.Math.sin


object AstroMath {

    fun frac(x: Double): Double {
        return x - Math.floor(x)
    }

    infix fun Double.modulo(y: Double) = y * frac (this / y)

    /**
     * Преобразование экваториальных координат в эклиптические
     * Transformation of equatorial to ecliptical coordinates
     * @param T Время в юлианских столетиях с эпохи J2000
     * @return Матрица преобразования
     */
    fun getMatrixEquatorialToEclipticCoordinates(T: Double): Matrix {
        val eps = Math.toRadians(23.43929111 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T / 3600.0)
        return Matrix(X, eps)
    }

    /**
     * Преобразование эклиптических координат в экваториальные
     * @param T Время в юлианских столетиях с эпохи J2000
     * @return Матрица преобразования
     */
    fun getMatrixEclipticToEquatorialCoordinates(T: Double): Matrix {
        return Matrix.transpose(getMatrixEquatorialToEclipticCoordinates(T))
    }

    /**
     * Прецессия в эклиптических координатах Precession of ecliptic coordinates
     * @param T1 Исходная эпоха в J2000
     * @param T2 Требуемая эпоха в J2000
     * @return Матрица преобразования координат
     */
    fun getMatrixEclipticPrecession(T1: Double, T2: Double): Matrix {
        val dT = T2 - T1
        val Pi: Double
        val pi: Double
        val p_a: Double
        Pi = 174.876383889 * RAD + ((3289.4789 + 0.60622 * T1) * T1 + (-869.8089 - 0.50491 * T1 + 0.03536 * dT) * dT) / ARCS
        pi = (47.0029 - (0.06603 - 0.000598 * T1) * T1 + (-0.03302 + 0.000598 * T1 + 0.000060 * dT) * dT) * dT / ARCS
        p_a = (5029.0966 + (2.22226 - 0.000042 * T1) * T1 + (1.11113 - 0.000042 * T1 - 0.000006 * dT) * dT) * dT / ARCS
        return Matrix(Z, -(Pi + p_a)) * Matrix(X, pi) * Matrix(Z, Pi)
    }

    /**
     * Прецессия в экваториальных координатах Precession of equatorial coordinates
     * @param T1 Исходная эпоха в J2000
     * @param T2 Требуемая эпоха в J2000
     * @return Матрица преобразования координат
     */
    fun getMatrixEquatorialPrecession(T1: Double, T2: Double): Matrix {
        val dT = T2 - T1
        val zeta: Double
        val z: Double
        val theta: Double
        zeta = (2306.2181 + (1.39656 - 0.000139 * T1) * T1 + (0.30188 - 0.000344 * T1 + 0.017998 * dT) * dT) * dT / ARCS
        z = zeta + (0.79280 + 0.000411 * T1 + 0.000205 * dT) * dT * dT / ARCS
        theta = (2004.3109 - (0.85330 + 0.000217 * T1) * T1 - (0.42665 + 0.000217 * T1 + 0.041833 * dT) * dT) * dT / ARCS
        return Matrix(Z, -z) * Matrix(Y, theta) * Matrix(Z, -zeta)
    }

    /**
     * Преобразование от средних к истинным экватору и равноденствию
     * Transformation from mean to true equator and equinox
     * @param T Время в юлианских столетиях от эпохи J2000
     * @return матрица нутации Nutation matrix
     */
    fun getMatrixNutation(T: Double): Matrix {
        val ls: Double = PI2 * frac(0.993133 + 99.997306 * T)
        val D: Double = PI2 * frac(0.827362 + 1236.853087 * T)
        val F: Double = PI2 * frac(0.259089 + 1342.227826 * T)
        val N: Double = PI2 * frac(0.347346 - 5.372447 * T)
        val eps: Double = 0.4090928 - 2.2696E-4 * T
        val dpsi = (-17.200 * sin(N) - 1.319 * sin(2 * (F - D + N)) - 0.227 * sin(2 * (F + N)) + 0.206 * sin(2 * N) + 0.143 * sin(ls)) / ARCS
        val deps = (+9.203 * cos(N) + 0.574 * cos(2 * (F - D + N)) + 0.098 * cos(2 * (F + N)) - 0.090 * cos(2 * N)) / ARCS
        return Matrix(X, -eps - deps) * Matrix(Z, -dpsi) * Matrix(X, eps)
    }
}