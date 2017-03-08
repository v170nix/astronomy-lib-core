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

package net.arwix.astronomy.vsop87

import net.arwix.astronomy.core.Epoch
import net.arwix.astronomy.core.coordinates.EclipticCoordinates
import net.arwix.astronomy.core.vector.RectangularVector
import net.arwix.astronomy.core.vector.Vector


abstract class VsopData(private val epoch: Epoch, private val obj: VsopObject) : EclipticCoordinates<VsopObject> {

    override fun getIdObject() = obj

    override fun getEpoch(): Epoch = epoch

    override fun getEclipticCoordinates(T: Double): Vector {
        val innerT = T / 10.0
        return RectangularVector(X0(innerT) + X1(innerT) + X2(innerT) + X3(innerT) + X4(innerT) + X5(innerT),
                Y0(innerT) + Y1(innerT) + Y2(innerT) + Y3(innerT) + Y4(innerT) + Y5(innerT),
                Z0(innerT) + Z1(innerT) + Z2(innerT) + Z3(innerT) + Z4(innerT) + Z5(innerT))
    }

    abstract fun X0(T: Double): Double
    abstract fun X1(T: Double): Double
    abstract fun X2(T: Double): Double
    abstract fun X3(T: Double): Double
    abstract fun X4(T: Double): Double
    abstract fun X5(T: Double): Double

    abstract fun Y0(T: Double): Double
    abstract fun Y1(T: Double): Double
    abstract fun Y2(T: Double): Double
    abstract fun Y3(T: Double): Double
    abstract fun Y4(T: Double): Double
    abstract fun Y5(T: Double): Double

    abstract fun Z0(T: Double): Double
    abstract fun Z1(T: Double): Double
    abstract fun Z2(T: Double): Double
    abstract fun Z3(T: Double): Double
    abstract fun Z4(T: Double): Double
    abstract fun Z5(T: Double): Double
}