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

/**
 * @param f        функция реализующа интерфейс [Function]
 * @param a        начальная точка отразка
 * @param b        конечная точка отрезка
 * @param e        желаемая ночность
 * @param maxSteps максимальное число шагов
 */
class SearchExtremumGoldenMethod(private val a: Double, private val b: Double, private val e: Double, private val maxSteps: Int, private val function: (x: Double) -> Double)
//  this.function = f
{

    companion object {
        val GOLDEN_RATIO = 0.5 + Math.sqrt(5.0) / 2.0
    }

    val min: Double by lazy { doMin(a, b) }
    val max: Double by lazy { doMax(a, b) }

    private fun doMax(a: Double, b: Double): Double {
        var a = a
        var b = b
        var step = 0
        do {
            step++
            val d = (b - a) / GOLDEN_RATIO
            val x1 = b - d
            val x2 = a + d
            val y1 = function(x1)
            val y2 = function(x2)
            if (y1 <= y2) {
                a = x1
            } else {
                b = x2
            }
        } while (Math.abs(a - b) > this.e && step < maxSteps)
        return (a + b) / 2.0
    }

    private fun doMin(a: Double, b: Double): Double {
        var a = a
        var b = b
        var step = 0
        do {
            step++
            val d = (b - a) / GOLDEN_RATIO
            val x1 = b - d
            val x2 = a + d
            val y1 = function(x1)
            val y2 = function(x2)
            if (y1 >= y2) {
                a = x1
            } else {
                b = x2
            }
        } while (Math.abs(a - b) > this.e && step < maxSteps)
        return (a + b) / 2.0
    }


}