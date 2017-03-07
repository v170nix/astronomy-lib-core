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

import net.arwix.astronomy.core.calendar.CalendarMath
import org.junit.After
import org.junit.Before
import org.junit.Test

import org.junit.Assert.*
import java.util.*
import java.util.Calendar.*
import java.util.Locale
import net.arwix.astronomy.core.vector.VectorType
import net.arwix.astronomy.core.vector.SphericalVector




class AstroMathTest {

    var t: Double =0.0

    @Before
    fun setUp() {
        val calendar =  Calendar.getInstance(TimeZone.getTimeZone("GMT")).apply {
            this[YEAR] = 2014
            this[MONTH] = 8
            this[DAY_OF_MONTH] = 17
            this[HOUR_OF_DAY] = 0
            this[MINUTE] = 0
            this[SECOND] = 0
            this[MILLISECOND] = 0
        }
        val mjd = CalendarMath.getMJD(calendar)
        val deltaT = CalendarMath.getDeltaTofDay(calendar)
        t =CalendarMath.getJT(mjd + deltaT)
    }

    @After
    fun tearDown() {

    }

    @Test
    fun getMatrixEquatorialToEclipticCoordinates() {
        val def = arrayOf(
                doubleArrayOf(1.0, 0.0, 0.0),
                doubleArrayOf(0.0, 0.9174953412097739, 0.39774652589100074),
                doubleArrayOf(0.0, -0.39774652589100074, 0.9174953412097739)
        )
        val m = AstroMath.getMatrixEquatorialToEclipticCoordinates(t).toArray()
        assertArrayEquals("not equal", def[0], m[0], 1E-14)
        assertArrayEquals("not equal", def[1], m[1], 1E-14)
        assertArrayEquals("not equal", def[2], m[2], 1E-14)
    }

    @Test
    fun getMatrixEclipticToEquatorialCoordinates() {
        val def = arrayOf(
                doubleArrayOf(1.0, 0.0, 0.0),
                doubleArrayOf(0.0, 0.9174953412097739, -0.39774652589100074),
                doubleArrayOf(0.0, 0.39774652589100074, 0.9174953412097739)
        )
        val m = AstroMath.getMatrixEclipticToEquatorialCoordinates(t).toArray()
        assertArrayEquals("not equal", def[0], m[0], 1E-14)
        assertArrayEquals("not equal", def[1], m[1], 1E-14)
        assertArrayEquals("not equal", def[2], m[2], 1E-14)
    }

    @Test
    fun getMatrixEclipticPrecession() {
        val def = arrayOf(
                doubleArrayOf(0.9999935686776528, -0.0035864460064648773, -2.8940355041349935E-6),
                doubleArrayOf(0.0035864459078496086, 0.9999935681247464, -3.338994676584286E-5),
                doubleArrayOf(3.013768131294084E-6, 3.3379352722541076E-5, 0.999999999438368)
        )
        val m = AstroMath.getMatrixEclipticPrecession(0.0, t).toArray()
        assertArrayEquals("not equal", def[0], m[0], 1E-14)
        assertArrayEquals("not equal", def[1], m[1], 1E-14)
        assertArrayEquals("not equal", def[2], m[2], 1E-14)
    }

    @Test
    fun getMatrixEquatorialPrecession() {
        val def = arrayOf(
                doubleArrayOf(0.999993568677741, -0.0032893486856280085, -0.0014292614807592046),
                doubleArrayOf(0.0032893486855430664, 0.9999945900752159, -2.350736676170233E-6),
                doubleArrayOf(0.0014292614809546928, -2.3506178152075035E-6, 0.9999989786025252)
        )
        val m = AstroMath.getMatrixEquatorialPrecession(0.0, t).toArray()
        assertArrayEquals("not equal", def[0], m[0], 1E-14)
        assertArrayEquals("not equal", def[1], m[1], 1E-14)
        assertArrayEquals("not equal", def[2], m[2], 1E-14)
    }

    @Test
    fun getMatrixNutation() {
        val def = arrayOf(
                doubleArrayOf(0.9999999995457635, -2.765412347493946E-5, -1.1988433012695878E-5),
                doubleArrayOf(2.7654600894686653E-5,  0.9999999988245951, 3.982502882166061E-5),
                doubleArrayOf(1.1987331672340185E-5, -3.98253603389187E-5, 0.9999999991351223)
        )
        val m = AstroMath.getMatrixNutation(t).toArray()
        assertArrayEquals("not equal", def[0], m[0], 1E-14)
        assertArrayEquals("not equal", def[1], m[1], 1E-14)
        assertArrayEquals("not equal", def[2], m[2], 1E-14)
    }

}