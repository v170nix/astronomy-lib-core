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

import net.arwix.astronomy.core.calendar.getJT
import net.arwix.astronomy.core.vector.Matrix
import net.arwix.astronomy.testing.NutationMatrix
import org.junit.Assert.assertArrayEquals
import org.junit.Before
import org.junit.Test
import java.util.*

class AstronomyMatrixTest {
    var t: Double = 0.0

    @Before
    fun setUp() {
        val calendar = Calendar.getInstance(TimeZone.getTimeZone("GMT")).apply {
            this[Calendar.YEAR] = 2014
            this[Calendar.MONTH] = 8
            this[Calendar.DAY_OF_MONTH] = 17
            this[Calendar.HOUR_OF_DAY] = 0
            this[Calendar.MINUTE] = 0
            this[Calendar.SECOND] = 0
            this[Calendar.MILLISECOND] = 0
        }
        t = calendar.getJT(true)
    }

    @Test
    fun createTransformationCoordinates() {
        var def = arrayOf(
                doubleArrayOf(1.0, 0.0, 0.0),
                doubleArrayOf(0.0, 0.9174953412097739, 0.39774652589100074),
                doubleArrayOf(0.0, -0.39774652589100074, 0.9174953412097739)
        )
        var m = AstronomyMatrix.createTransformationCoordinates(t, AstronomyMatrix.Coordinates.EQUATORIAL, AstronomyMatrix.Coordinates.ECLIPTIC).toArray()
        assertArrayEquals("not equal", def[0], m[0], 1E-14)
        assertArrayEquals("not equal", def[1], m[1], 1E-14)
        assertArrayEquals("not equal", def[2], m[2], 1E-14)

        def = arrayOf(
                doubleArrayOf(1.0, 0.0, 0.0),
                doubleArrayOf(0.0, 0.9174953412097739, -0.39774652589100074),
                doubleArrayOf(0.0, 0.39774652589100074, 0.9174953412097739)
        )
        m = AstronomyMatrix.createTransformationCoordinates(t, AstronomyMatrix.Coordinates.ECLIPTIC, AstronomyMatrix.Coordinates.EQUATORIAL).toArray()
        assertArrayEquals("not equal", def[0], m[0], 1E-14)
        assertArrayEquals("not equal", def[1], m[1], 1E-14)
        assertArrayEquals("not equal", def[2], m[2], 1E-14)
    }

    @Test
    fun createPrecession() {
        var def = arrayOf(
                doubleArrayOf(0.9999935686776528, -0.0035864460064648773, -2.8940355041349935E-6),
                doubleArrayOf(0.0035864459078496086, 0.9999935681247464, -3.338994676584286E-5),
                doubleArrayOf(3.013768131294084E-6, 3.3379352722541076E-5, 0.999999999438368)
        )
        var m = AstronomyMatrix.createPrecession(0.0, t, AstronomyMatrix.Coordinates.ECLIPTIC).toArray()
//        var m1 = Precession.get(t).toArray()
//        m1.forEach { it.forEach { System.out.println(it.toString()) } }
        assertArrayEquals("createPrecession not equal", def[0], m[0], 2E-5)
        assertArrayEquals("e not equal", def[1], m[1], 2E-5)
        assertArrayEquals("e not equal", def[2], m[2], 2E-5)

        def = arrayOf(
                doubleArrayOf(0.999993568677741, -0.0032893486856280085, -0.0014292614807592046),
                doubleArrayOf(0.0032893486855430664, 0.9999945900752159, -2.350736676170233E-6),
                doubleArrayOf(0.0014292614809546928, -2.3506178152075035E-6, 0.9999989786025252)
        )
        m = AstronomyMatrix.createPrecession(0.0, t, AstronomyMatrix.Coordinates.EQUATORIAL).toArray()
        assertArrayEquals("not equal", def[0], m[0], 1E-14)
        assertArrayEquals("not equal", def[1], m[1], 1E-14)
        assertArrayEquals("not equal", def[2], m[2], 1E-14)
    }

    @Test
    fun createNutation() {
//        val t = 0.4779055751655893
        val def = arrayOf(
                doubleArrayOf(0.9999999995457635, -2.765412347493946E-5, -1.1988433012695878E-5),
                doubleArrayOf(2.7654600894686653E-5, 0.9999999988245951, 3.982502882166061E-5),
                doubleArrayOf(1.1987331672340185E-5, -3.98253603389187E-5, 0.9999999991351223)
        )
        val m = AstronomyMatrix.createNutation(t).toArray()
        assertArrayEquals("not equal", def[0], m[0], 1E-14)
        assertArrayEquals("not equal", def[1], m[1], 1E-14)
        assertArrayEquals("not equal", def[2], m[2], 1E-14)
        NutationMatrix.getNutationAngles(t)
    }

    private fun printMatrix(m: Matrix) {
        (0..2).forEach { i ->
            (0..2).forEach { j ->
                System.out.print(m[i, j].toString() + "  ");
            }
            System.out.println();
        }
    }

}