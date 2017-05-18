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

class RectangularVector(@JvmField var x: Double, @JvmField var y: Double, @JvmField var z: Double): Vector() {

    constructor(array: DoubleArray) : this(array[0], array[1], array[2])
    constructor(): this(0.0, 0.0, 0.0)
    constructor(vector: Vector): this() {
        set(vector)
    }

    override fun getType() = VectorType.RECTANGULAR

    override fun toArray() = doubleArrayOf(x ,y ,z)

    override fun set(vector: Vector) {
        val rectangularVector = (if (vector.getType() === VectorType.RECTANGULAR) vector else convert(vector, VectorType.RECTANGULAR)) as RectangularVector
        set(rectangularVector.x, rectangularVector.y, rectangularVector.z)
    }

    fun set(x: Double, y: Double, z: Double) {
        this.x = x
        this.y = y
        this.z = z
    }

    override fun set(i: Int, element: Double) {
        when(i) {
            0 -> x = element
            1 -> y = element
            2 -> z = element
            else -> throw IndexOutOfBoundsException()
        }
    }

    override fun get(index: Int): Double = when(index) {
        0 -> x
        1 -> y
        2 -> z
        else -> throw IndexOutOfBoundsException()
    }

    override fun component1() = x

    override fun component2() = y

    override fun component3() = z

}