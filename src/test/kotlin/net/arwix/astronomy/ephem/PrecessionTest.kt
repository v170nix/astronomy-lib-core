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

package net.arwix.astronomy.ephem

import net.arwix.astronomy.core.MJD_J2000
import net.arwix.astronomy.core.calendar.getJT
import net.arwix.astronomy.core.vector.Matrix
import net.arwix.astronomy.core.vector.RectangularVector
import org.junit.Before
import org.junit.Test

/**
 * Created by dreamirus on 14.05.17.
 */
class PrecessionTest {
    @Before
    fun setUp() {

    }

    @Test
    fun getVodak() {
//        val vector = Precession.precessionVondrak2011(getJT(MJD_J2000+2000), RectangularVector(0.1, 0.7, 0.6) ).getVectorOfType(VectorType.RECTANGULAR) as RectangularVector
//
//
//        System.out.println("new Precession ${vector[0]} ${vector[1]} ${vector[2]}")


    }

    @Test
    fun precessFromJ2000() {
//        val vector = Precession.precessFromJ2000(getJT(MJD_J2000 + 200), RectangularVector(0.1, 0.7, 0.6));
//
//        System.out.println("new Precession ${vector[0]} ${vector[1]} ${vector[2]}")
    }

    @Test
    fun precessToJ2000() {
//        val vector = Precession.precessToJ2000(getJT(MJD_J2000 + 200), RectangularVector(0.1, 0.7, 0.6));
//
//        System.out.println("new Precession ${vector[0]} ${vector[1]} ${vector[2]}")
    }

    @Test
    fun precess() {
        val vector = Precession.precess(getJT(MJD_J2000 - 2000), getJT(MJD_J2000 + 2000), RectangularVector(0.1, 0.7, 0.6));

        System.out.println("new Precession ${vector[0]} ${vector[1]} ${vector[2]}")
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