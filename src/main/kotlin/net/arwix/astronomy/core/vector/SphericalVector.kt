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


class SphericalVector(@JvmField var phi: Double, @JvmField var theta: Double, @JvmField var r: Double) : Vector() {

    constructor() : this(0.0, 0.0, 0.0)
    constructor(phi: Double, theta: Double) : this(phi, theta, 1.0)
    constructor(vector: Vector) : this() {
        set(vector)
    }


    override fun getType() = VectorType.SPHERICAL

    override fun toArray() = doubleArrayOf(phi, theta, r)

    fun set(phi: Double, theta: Double, r: Double) {
        this.phi = phi
        this.theta = theta
        this.r = r
    }

    override fun set(vector: Vector) {
        val sphericalVector = (if (vector.getType() === VectorType.SPHERICAL) vector else convert(vector, VectorType.SPHERICAL)) as SphericalVector
        set(sphericalVector.phi, sphericalVector.theta, sphericalVector.r)
    }

    override fun set(i: Int, element: Double) {
        when(i) {
            0 -> phi = element
            1 -> theta = element
            2 -> r = element
            else -> throw IndexOutOfBoundsException()
        }
    }

    override fun get(index: Int): Double = when(index) {
            0 -> phi
            1 -> theta
            2 -> r
            else -> throw IndexOutOfBoundsException()
        }



}