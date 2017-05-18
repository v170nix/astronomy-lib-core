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

package net.arwix.astronomy.ephemeris

import net.arwix.astronomy.core.ARCSEC_TO_RAD
import net.arwix.astronomy.core.PI2
import net.arwix.astronomy.core.vector.Matrix
import net.arwix.astronomy.core.vector.RectangularVector
import net.arwix.astronomy.core.vector.Vector
import net.arwix.astronomy.core.vector.VectorType


object Precession {

    private val xypol = arrayOf(
            doubleArrayOf(8473.343527, 5042.7980307, -0.00740913, 289E-9),
            doubleArrayOf(84283.175915, -0.4436568, 0.00000146, 151E-9),
            doubleArrayOf(-19.657270, 0.0790159, 0.00001472, -61E-9))

    private val xyper = arrayOf(
            doubleArrayOf(402.90, -22206.325946, 1267.727824, -3243.236469, -8571.476251),
            doubleArrayOf(256.75, 12236.649447, 1702.324248, -3969.723769, 5309.796459),
            doubleArrayOf(292.00, -1589.008343, -2970.553839, 7099.207893, -610.393953),
            doubleArrayOf(537.22, 2482.103195, 693.790312, -1903.696711, 923.201931),
            doubleArrayOf(241.45, 150.322920, -14.724451, 146.435014, 3.759055),
            doubleArrayOf(375.22, -13.632066, -516.649401, 1300.630106, -40.691114),
            doubleArrayOf(157.87, 389.437420, -356.794454, 1727.498039, 80.437484),
            doubleArrayOf(274.20, 2031.433792, -129.552058, 299.854055, 807.300668),
            doubleArrayOf(203.00, 363.748303, 256.129314, -1217.125982, 83.712326),
            doubleArrayOf(440.00, -896.747562, 190.266114, -471.367487, -368.654854),
            doubleArrayOf(170.72, -926.995700, 95.103991, -441.682145, -191.881064),
            doubleArrayOf(713.37, 37.070667, -332.907067, -86.169171, -4.263770),
            doubleArrayOf(313.00, -597.682468, 131.337633, -308.320429, -270.353691),
            doubleArrayOf(128.38, 66.282812, 82.731919, -422.815629, 11.602861))

    private val zper = arrayOf(
            doubleArrayOf(402.90, -13765.924050, -2206.967126),
            doubleArrayOf(256.75, 13511.858383, -4186.752711),
            doubleArrayOf(292.00, -1455.229106, 6737.949677),
            doubleArrayOf(537.22, 1054.394467, -856.922846),
            doubleArrayOf(375.22, -112.300144, 957.149088),
            doubleArrayOf(157.87, 202.769908, 1709.440735),
            doubleArrayOf(274.20, 1936.050095, 154.425505),
            doubleArrayOf(202.00, 327.517465, -1049.071786),
            doubleArrayOf(440.00, -655.484214, -243.520976),
            doubleArrayOf(170.72, -891.898637, -406.539008),
            doubleArrayOf(315.00, -494.780332, -301.504189),
            doubleArrayOf(136.32, 585.492621, 41.348740),
            doubleArrayOf(128.38, -333.322021, -446.656435),
            doubleArrayOf(490.00, 110.512834, 142.525186))


    fun getMatrix(T: Double): Matrix {
        var x = 0.0
        var y = 0.0
        var z = 0.0
        var w = PI2 * T

        for (i in 0..13) {
            var a = w / xyper[i][0]
            var s = Math.sin(a)
            var c = Math.cos(a)
            x += c * xyper[i][1] + s * xyper[i][3]
            y += c * xyper[i][2] + s * xyper[i][4]

            a = w / zper[i][0]
            s = Math.sin(a)
            c = Math.cos(a)
            z += c * zper[i][1] + s * zper[i][2]
        }

        w = 1.0
        for (j in 0..3) {
            x += xypol[0][j] * w
            y += xypol[1][j] * w
            z += xypol[2][j] * w
            w *= T
        }
        x *= ARCSEC_TO_RAD
        y *= ARCSEC_TO_RAD
        z *= ARCSEC_TO_RAD
        //w = x * x + y * y;
        //double z = 0.0;
        //if (w < 1.0) z = -Math.sqrt(1.0 - w);

        val PSIA = x
        val OMEGAA = y
        val CHIA = z
        val EPS0 = 84381.406 * ARCSEC_TO_RAD
        val SA = Math.sin(EPS0)
        val CA = Math.cos(EPS0)
        val SB = Math.sin(-PSIA)
        val CB = Math.cos(-PSIA)
        val SC = Math.sin(-OMEGAA)
        val CC = Math.cos(-OMEGAA)
        val SD = Math.sin(CHIA)
        val CD = Math.cos(CHIA)

        // COMPUTE ELEMENTS OF PRECESSION ROTATION MATRIX
        // EQUIVALENT TO R3(CHI_A)R1(-OMEGA_A)R3(-PSI_A)R1(EPSILON_0)
        val XX = CD * CB - SB * SD * CC
        val YX = CD * SB * CA + SD * CC * CB * CA - SA * SD * SC
        val ZX = CD * SB * SA + SD * CC * CB * SA + CA * SD * SC
        val XY = -SD * CB - SB * CD * CC
        val YY = -SD * SB * CA + CD * CC * CB * CA - SA * CD * SC
        val ZY = -SD * SB * SA + CD * CC * CB * SA + CA * CD * SC
        val XZ = SB * SC
        val YZ = -SC * CB * CA - SA * CC
        val ZZ = -SC * CB * SA + CC * CA

        return Matrix(doubleArrayOf(XX, YX, ZX), doubleArrayOf(XY, YY, ZY), doubleArrayOf(XZ, YZ, ZZ))
    }

    /**
     * Precession following Vondrak et al. 2011. See A&amp;A 534, A22.
     *
     * @param JD0 Julian day of input vector (equatorial rectangular).
     * @param JD Julian day of output. Either JD or JD0 must be equal to Constant.J2000.
     * @param R Input vector.
     * @param isFromJ2000ToApparent если true из эпохи J2000 в текущую, иначе из текущей в J2000
     * @return Vector referred to mean equinox and equator of JD.
     */
    fun precessionVondrak2011(T: Double, R: Vector, isFromJ2000ToApparent: Boolean = true): Vector {

        var x = 0.0
        var y = 0.0
        var z = 0.0
        var w = PI2 * T

        for (i in 0..13) {
            var a = w / xyper[i][0]
            var s = Math.sin(a)
            var c = Math.cos(a)
            x += c * xyper[i][1] + s * xyper[i][3]
            y += c * xyper[i][2] + s * xyper[i][4]

            a = w / zper[i][0]
            s = Math.sin(a)
            c = Math.cos(a)
            z += c * zper[i][1] + s * zper[i][2]
        }

        w = 1.0
        for (j in 0..3) {
            x += xypol[0][j] * w
            y += xypol[1][j] * w
            z += xypol[2][j] * w
            w *= T
        }
        x *= ARCSEC_TO_RAD
        y *= ARCSEC_TO_RAD
        z *= ARCSEC_TO_RAD
        //w = x * x + y * y;
        //double z = 0.0;
        //if (w < 1.0) z = -Math.sqrt(1.0 - w);

        val PSIA = x
        val OMEGAA = y
        val CHIA = z
        val EPS0 = 84381.406 * ARCSEC_TO_RAD
        val SA = Math.sin(EPS0)
        val CA = Math.cos(EPS0)
        val SB = Math.sin(-PSIA)
        val CB = Math.cos(-PSIA)
        val SC = Math.sin(-OMEGAA)
        val CC = Math.cos(-OMEGAA)
        val SD = Math.sin(CHIA)
        val CD = Math.cos(CHIA)

        // COMPUTE ELEMENTS OF PRECESSION ROTATION MATRIX
        // EQUIVALENT TO R3(CHI_A)R1(-OMEGA_A)R3(-PSI_A)R1(EPSILON_0)
        val XX = CD * CB - SB * SD * CC
        val YX = CD * SB * CA + SD * CC * CB * CA - SA * SD * SC
        val ZX = CD * SB * SA + SD * CC * CB * SA + CA * SD * SC
        val XY = -SD * CB - SB * CD * CC
        val YY = -SD * SB * CA + CD * CC * CB * CA - SA * CD * SC
        val ZY = -SD * SB * SA + CD * CC * CB * SA + CA * CD * SC
        val XZ = SB * SC
        val YZ = -SC * CB * CA - SA * CC
        val ZZ = -SC * CB * SA + CC * CA

        val precMatrix = Matrix(doubleArrayOf(XX, YX, ZX), doubleArrayOf(XY, YY, ZY), doubleArrayOf(XZ, YZ, ZZ))

        val result = R * precMatrix

        System.out.println("new Precession ${result[0]} ${result[1]} ${result[2]}")

        var px = 0.0
        var py = 0.0
        var pz = 0.0

        val r = R.getVectorOfType(VectorType.RECTANGULAR) as RectangularVector

        if (isFromJ2000ToApparent) {
            // PERFORM ROTATION FROM J2000.0 TO EPOCH
            px = XX * r.x + YX * r.y + ZX * r.z
            py = XY * r.x + YY * r.y + ZY * r.z
            pz = XZ * r.x + YZ * r.y + ZZ * r.z
        } else {
            // PERFORM ROTATION FROM EPOCH TO J2000.0
            px = XX * r.x + XY * r.y + XZ * r.z
            py = YX * r.x + YY * r.y + YZ * r.z
            pz = ZX * r.x + ZY * r.y + ZZ * r.z
        }

        return RectangularVector(px, py, pz)
    }

    /**
     * Precess rectangular equatorial coordinates from J2000 epoch.
     *
     * @param R Array with x, y, z.
     * @return Array with corrected x, y, z.
     */
    fun precessFromJ2000(T: Double, R: Vector): Vector {
        //      return precessionVondrak2011(T, R, true)
        return getMatrix(T) * R
    }

    /**
     * Precess rectangular equatorial coordinates to J2000 epoch.
     *
     * @param JD Equinox of the input in Julian day (TT).
     * @param R Array with x, y, z.
     * @param eph The ephemeris properties.
     * @return Array with corrected x, y, z.
     * @throws JPARSECException If an error occurs.
     */
    fun precessToJ2000(T: Double, R: Vector): Vector {
        return R * getMatrix(T)
//        return precessionVondrak2011(T, R, false)
    }

    /**
     * Precess rectangular equatorial coordinates between two epochs.
     *
     * @param fromDatejT Equinox of the input in Julian day (TT).
     * @param toDatejT Equinox of the output in Julian day (TT).
     * @param R Array with x, y, z.
     * @return Array with corrected x, y, z.
     */
    fun precess(fromDatejT: Double, toDatejT: Double, R: Vector): Vector {
        return getMatrix(toDatejT) * (R * getMatrix(fromDatejT))
//        val v = precessToJ2000(fromDatejT, R)
//        return precessFromJ2000(toDatejT, v)
    }

}