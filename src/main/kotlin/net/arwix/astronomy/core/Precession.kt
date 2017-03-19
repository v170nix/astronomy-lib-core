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

import net.arwix.astronomy.core.vector.Matrix
import java.lang.Math.cos
import java.lang.Math.sin


object Precession {

    val precVals = arrayListOf(
            // 1/Pn         psi_A:Cn       om_A:Cn       chi_A:Cn         psi_A:Sn       om_A:Sn       chi_A:Sn
            doubleArrayOf(1.0 / 402.90, -22206.325946, 1267.727824, -13765.924050, -3243.236469, -8571.476251, -2206.967126),
            doubleArrayOf(1.0 / 256.75, 12236.649447, 1702.324248, 13511.858383, -3969.723769, 5309.796459, -4186.752711),
            doubleArrayOf(1.0 / 292.00, -1589.008343, -2970.553839, -1455.229106, 7099.207893, -610.393953, 6737.949677),
            doubleArrayOf(1.0 / 537.22, 2482.103195, 693.790312, 1054.394467, -1903.696711, 923.201931, -856.922846),
            doubleArrayOf(1.0 / 241.45, 150.322920, -14.724451, 0.0, 146.435014, 3.759055, 0.0),
            doubleArrayOf(1.0 / 375.22, -13.632066, -516.649401, -112.300144, 1300.630106, -40.691114, 957.149088),
            doubleArrayOf(1.0 / 157.87, 389.437420, -356.794454, 202.769908, 1727.498039, 80.437484, 1709.440735),
            doubleArrayOf(1.0 / 274.20, 2031.433792, -129.552058, 1936.050095, 299.854055, 807.300668, 154.425505),
            doubleArrayOf(1.0 / 203.00, 363.748303, 256.129314, 0.0, -1217.125982, 83.712326, 0.0),
            doubleArrayOf(1.0 / 440.00, -896.747562, 190.266114, -655.484214, -471.367487, -368.654854, -243.520976),
            doubleArrayOf(1.0 / 170.72, -926.995700, 95.103991, -891.898637, -441.682145, -191.881064, -406.539008),
            doubleArrayOf(1.0 / 713.37, 37.070667, -332.907067, 0.0, -86.169171, -4.263770, 0.0),
            doubleArrayOf(1.0 / 313.00, -597.682468, 131.337633, 0.0, -308.320429, -270.353691, 0.0),
            doubleArrayOf(1.0 / 128.38, 66.282812, 82.731919, -333.322021, -422.815629, 11.602861, -446.656435),
            doubleArrayOf(1.0 / 202.00, 0.0, 0.0, 327.517465, 0.0, 0.0, -1049.071786),
            doubleArrayOf(1.0 / 315.00, 0.0, 0.0, -494.780332, 0.0, 0.0, -301.504189),
            doubleArrayOf(1.0 / 136.32, 0.0, 0.0, 585.492621, 0.0, 0.0, 41.348740),
            doubleArrayOf(1.0 / 490.00, 0.0, 0.0, 110.512834, 0.0, 0.0, 142.525186))

    val p_epsVals = arrayListOf(
            //  1/Pn         p_A:Cn     eps_A:Cn        p_A:Sn      eps_A:Sn
            doubleArrayOf(1.0 / 409.90, -6908.287473, 753.872780, -2845.175469, -1704.720302),
            doubleArrayOf(1.0 / 396.15, -3198.706291, -247.805823, 449.844989, -862.308358),
            doubleArrayOf(1.0 / 537.22, 1453.674527, 379.471484, -1255.915323, 447.832178),
            doubleArrayOf(1.0 / 402.90, -857.748557, -53.880558, 886.736783, -889.571909),
            doubleArrayOf(1.0 / 417.15, 1173.231614, -90.109153, 418.887514, 190.402846),
            doubleArrayOf(1.0 / 288.92, -156.981465, -353.600190, 997.912441, -56.564991),
            doubleArrayOf(1.0 / 4043.00, 371.836550, -63.115353, -240.979710, -296.222622),
            doubleArrayOf(1.0 / 306.00, -216.619040, -28.248187, 76.541307, -75.859952),
            doubleArrayOf(1.0 / 277.00, 193.691479, 17.703387, -36.788069, 67.473503),
            doubleArrayOf(1.0 / 203.00, 11.891524, 38.911307, -170.964086, 3.014055))

    fun get(T: Double): Matrix {
        val T2pi = T * (2.0 * Math.PI) // Julian centuries from J2000.0, premultiplied by 2Pi
        // these are actually small greek letters in the papers.
        var Psi_A = 0.0
        var Omega_A = 0.0
        var Chi_A = 0.0
        var Epsilon_A = 0.0
        var p_A = 0.0; // currently unused. The data don't disturb.

        (0..17).forEach { i ->
            val phase = T2pi * precVals[i][0]
            val sin2piT_P = sin(phase)
            val cos2piT_P = cos(phase)
            Psi_A += precVals[i][1] * cos2piT_P + precVals[i][4] * sin2piT_P
            Omega_A += precVals[i][2] * cos2piT_P + precVals[i][5] * sin2piT_P
            Chi_A += precVals[i][3] * cos2piT_P + precVals[i][6] * sin2piT_P
        }

        (0..9).forEach { i ->
            val invP = p_epsVals[i][0]
            val phase = T2pi * invP;
            val sin2piT_P = sin(phase);
            val cos2piT_P = cos(phase);
            p_A += p_epsVals[i][1] * cos2piT_P + p_epsVals[i][3] * sin2piT_P;
            Epsilon_A += p_epsVals[i][2] * cos2piT_P + p_epsVals[i][4] * sin2piT_P
        }

        Psi_A += ((289.0E-9 * T - 0.00740913) * T + 5042.7980307) * T + 8473.343527
        Omega_A += ((151.0E-9 * T + 0.00000146) * T - 0.4436568) * T + 84283.175915
        Chi_A += ((-61.0E-9 * T + 0.00001472) * T + 0.0790159) * T - 19.657270
        Epsilon_A += ((-110.0E-9 * T - 0.00004039) * T + 0.3624445) * T + 84028.206305
        p_A += ((271.0e-9 * T - 0.00710733) * T + 5043.0520035) * T + 8134.017132;
        System.out.println("Psi " + Psi_A * ARC_SEC_2RAD)
        System.out.println("Omega_A " + Omega_A * ARC_SEC_2RAD)
        System.out.println("Chi_A " + Chi_A * ARC_SEC_2RAD)
        System.out.println("Epsilon_A " + Epsilon_A * ARC_SEC_2RAD)
        System.out.println("p_A " + p_A * ARC_SEC_2RAD)

        return Matrix(Matrix.Axis.Z, -Psi_A * ARC_SEC_2RAD) * Matrix(Matrix.Axis.X, -Omega_A * ARC_SEC_2RAD) * Matrix(Matrix.Axis.Z, Chi_A * ARC_SEC_2RAD);
//        c_psi_A     = ARC_SEC_2RAD * Psi_A
//        c_omega_A   = ARC_SEC_2RAD * Omega_A
//        c_chi_A     = ARC_SEC_2RAD * Chi_A
//        c_epsilon_A = ARC_SEC_2RAD * Epsilon_A;
    }

    fun getPrecessionAngleVondrakEpsilon(T: Double): Double {
        val T2pi = T * (2.0 * Math.PI) // Julian centuries from J2000.0, premultiplied by 2Pi
        return (((-110.0E-9 * T - 0.00004039) * T + 0.3624445) * T + 84028.206305 +
                (0..9).sumByDouble { i ->
                    val invP = p_epsVals[i][0]
                    val phase = T2pi * invP
                    val sin2piT_P = sin(phase)
                    val cos2piT_P = cos(phase)
                    p_epsVals[i][2] * cos2piT_P + p_epsVals[i][4] * sin2piT_P
                }) * ARC_SEC_2RAD
    }


}