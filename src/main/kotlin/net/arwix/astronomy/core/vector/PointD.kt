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


/**
 * PointD holds two double coordinates
 */
class PointD {
    @JvmField
    var x: Double = 0.0
    @JvmField
    var y: Double = 0.0

    constructor() {}

    constructor(x: Double, y: Double) {
        this.x = x
        this.y = y
    }

    /**
     * Set the point's x and y coordinates
     */
    fun set(x: Double, y: Double) {
        this.x = x
        this.y = y
    }

    /**
     * Set the point's x and y coordinates to the coordinates of p
     */
    fun set(p: PointD) {
        this.x = p.x
        this.y = p.y
    }

    fun negate() {
        x = -x
        y = -y
    }

    fun offset(dx: Double, dy: Double) {
        x += dx
        y += dy
    }

    /**
     * Returns true if the point's coordinates equal (x,y)
     */
    fun equals(x: Double, y: Double): Boolean {
        return this.x == x && this.y == y
    }

    /**
     * Return the euclidean distance from (0,0) to the point
     */
    fun length(): Double {
        return length(x, y)
    }

    companion object {

        /**
         * Returns the euclidean distance from (0,0) to (x,y)
         */
        fun length(x: Double, y: Double): Double {
            return Math.sqrt(x * x + y * y)
        }
    }
}