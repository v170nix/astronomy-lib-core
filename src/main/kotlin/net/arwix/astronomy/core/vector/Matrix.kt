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

package net.arwix.astronomy.core.vector

import java.lang.Math.cos
import java.lang.Math.sin


class Matrix() {

    enum class Axis {X, Y, Z}

    companion object {
        /**
         * матрицы поворта вокруг осей базиса
         * elementary rotations
         */
        private fun getRotateX(angle: Double): Matrix {
            val s = sin(angle)
            val c = cos(angle)
            return Matrix(
                    doubleArrayOf(1.0, 0.0, 0.0),
                    doubleArrayOf(0.0, c, s),
                    doubleArrayOf(0.0, -s, c)
            )
        }

        private fun getRotateY(angle: Double): Matrix {
            val s = sin(angle)
            val c = cos(angle)
            return Matrix(
                    doubleArrayOf(c, 0.0, -s),
                    doubleArrayOf(0.0, 1.0, 0.0),
                    doubleArrayOf(s, 0.0, c)
            )
        }

        private fun getRotateZ(angle: Double): Matrix {
            val s = sin(angle)
            val c = cos(angle)
            return Matrix(
                    doubleArrayOf(c, s, 0.0),
                    doubleArrayOf(-s, c, 0.0),
                    doubleArrayOf(0.0, 0.0, 1.0)
            )
        }

        operator fun invoke(axis: Axis, angle: Double): Matrix =
            when (axis) {
                Matrix.Axis.X -> getRotateX(angle)
                Matrix.Axis.Y -> getRotateY(angle)
                Matrix.Axis.Z -> getRotateZ(angle)
            }

        /**
         * Transpose of matrix
         */
        fun transpose(matrix: Matrix): Matrix {
            val out = Matrix()
            (0..2).forEach { i -> (0..2).forEach { j -> out[i, j] = matrix[j, i] } }
            for (i in 0..2)
                for (j in 0..2)
                    out[i, j] = matrix[j, i]
            return out
        }

        /**
         * матричное векторное умножение, матрицы на вектор
         * @param matrix матрица
         * @param vector вектор
         * @return новый вектор
         */
        internal fun timesMV(matrix: Matrix, vector: Vector): Vector {
            val v = (vector.getVectorOfType(VectorType.RECTANGULAR) as RectangularVector)
            return RectangularVector().apply {
                (0..2).forEach { i -> this[i] = (0..2).sumByDouble { j -> matrix[i, j] * v[j] } }
            }
        }

        /**
         * матричное векторное умножение, вектора на матрицу
         * @param vector вектор
         * @param matrix матрица
         * @return новый вектор
         */
        internal fun timesVM(vector: Vector, matrix: Matrix): Vector {
            val v = (vector.getVectorOfType(VectorType.RECTANGULAR) as RectangularVector)
            return RectangularVector().apply {
                (0..2).forEach { j -> this[j] = (0..2).sumByDouble { i -> v[i] * matrix[i, j] } }
            }
        }

        /**
         * матричное умножение
         * @param left  матрица
         * @param right матрица
         * @return новая матрица
         */
        private fun times(left: Matrix, right: Matrix): Matrix {
            return Matrix().apply {
                (0..2).forEach { i ->
                    (0..2).forEach { j -> this[i, j] = (0..2).sumByDouble { k -> left[i, k] * right[k, j] } }
                }
            }
        }

    }

    private var elements = arrayOf(doubleArrayOf(0.0, 0.0, 0.0), doubleArrayOf(0.0, 0.0, 0.0), doubleArrayOf(0.0, 0.0, 0.0))

    constructor(array1: DoubleArray, array2: DoubleArray, array3: DoubleArray) : this() {
        this[0] = array1
        this[1] = array2
        this[2] = array3
    }

    constructor(vector1: Vector, vector2: Vector, vector3: Vector) : this() {
        this[0] = vector1
        this[1] = vector2
        this[2] = vector3
    }

    fun toArray() = elements.copyOf()

    fun transpose() = Companion.transpose(this)

    operator fun times(right: Matrix) = Companion.times(this, right)

    operator fun times(right: Vector) = Companion.timesMV(this, right)

    operator fun timesAssign(matrix: Matrix) {
        elements = (this * matrix).elements
    }

    operator fun get(i: Int, j: Int) = elements[i][j]

    operator fun get(i: Int) = RectangularVector(elements[i])

    operator fun set(i: Int, array: DoubleArray) {
        elements[i] = array
    }

    operator fun set(i: Int, vector: Vector) {
        elements[i] = vector.getVectorOfType(VectorType.RECTANGULAR).toArray()
    }

    operator fun set(i: Int, j: Int, element: Double) {
        elements[i][j] = element
    }

}