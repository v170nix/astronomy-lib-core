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

package net.arwix.astronomy.math

import net.arwix.astronomy.core.vector.PointD


/**
 * Квадратичная интерполяция
 * Находит корни и экстремум параболы, интерполирующей три эквидистантных значения функции f при x = {-1, 0, 1}
 * @param yMinus Значение функции в точке x = -1
 * @param y0     Значение функции в точке x = 0
 * @param yPlus  Значение функции в точке x = +1
 */
abstract class QuadraticInterpolator {
    /**
     * Результат может быть и без корней тогда Roots = null или roots.count = 0
     *
     * @return {@link Result} корни и экстремум
     * элементы
     * 0 - абсцисса экстремума
     * 1 - ордината экстремума
     * 2,3 - корни
     * 4 - колличество найденых корней
     */
    sealed class Result(val extremum: PointD) {
        class Roots(extremum: PointD, val root1: Double, val root2: Double) : Result(extremum)
        class Root(extremum: PointD, val root: Double) : Result(extremum)
        class None(extremum: PointD) : Result(extremum)
    }

    companion object {
        /**
         * Квадратичная интерполяция
         * Находит корни и экстремум параболы, интерполирующей три эквидистантных значения функции f при x = {-1, 0, 1}
         * @param yMinus Значение функции в точке x = -1
         * @param y0     Значение функции в точке x = 0
         * @param yPlus  Значение функции в точке x = +1
         */
        fun getResult(y_minus: Double, c: Double, y_plus: Double): Result {
            // Коэффициенты итерполирующей параболы y=a*x^2+b*x+c
            val a = 0.5 * (y_plus + y_minus) - c
            val b = 0.5 * (y_plus - y_minus)

            // Находим экстремум
            val xe = -b / (2.0 * a)
            val extremum = PointD(xe, (a * xe + b) * xe + c)
            val dis = b * b - 4.0 * a * c // дискриминант уравнения y=a*x^2+b*x+c
            if (dis >= 0)
            // Парабола имеет корни roots
            {
                val dx = 0.5 * Math.sqrt(dis) / Math.abs(a)
                var root1 = extremum.x - dx
                val root2 = extremum.x + dx
                var count = 0

                if (Math.abs(root1) <= 1.0) count++
                if (Math.abs(root2) <= 1.0) count++
                if (root1 < -1.0) {
                    root1 = root2
                }
                if (count > 1) return Result.Roots(extremum, root1, root2) else return Result.Root(extremum, root1)
            } else
                return Result.None(extremum)
        }
    }

}