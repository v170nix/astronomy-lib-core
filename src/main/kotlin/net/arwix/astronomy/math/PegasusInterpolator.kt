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

package net.arwix.astronomy.math


//------------------------------------------------------------------------------
//
// Pegasus: Root finder using the Pegasus method
//
// Input:
//
//   PegasusFunct  Pointer to the function to be examined
//
//   LowerBound    Lower bound of search interval
//   UpperBound    Upper bound of search interval
//   Accuracy      Desired accuracy for the root
//
// Output:
//
//   Root          Root found (valid only if Success is true)
//   Success       Flag indicating success of the routine
//
// References:
//
//   Dowell M., Jarratt P., 'A modified Regula Falsi Method for Computing
//     the root of an equation', BIT 11, p.168-174 (1971).
//   Dowell M., Jarratt P., 'The "PEGASUS Method for Computing the root
//     of an equation', BIT 12, p.503-508 (1972).
//   G.Engeln-Muellges, F.Reutter, 'Formelsammlung zur Numerischen
//     Mathematik mit FORTRAN77-Programmen', Bibliogr. Institut,
//     Zuerich (1986).
//
// Notes:
//
//   Pegasus assumes that the root to be found is bracketed in the interval
//   [LowerBound, UpperBound]. Ordinates for these abscissae must therefore
//   have different signs.
//
//------------------------------------------------------------------------------
abstract class PegasusInterpolator {

    sealed class Result {
        data class Root(val root: Double) : Result()
        object None : Result()
    }

    companion object {
        fun getResult(function: (x: Double) -> Double,
                      lowerBound: Double,
                      upperBound: Double,
                      accuracy: Double,
                      maxSteps: Int = 30): Result {

            var x1 = lowerBound
            var x2 = upperBound
            var y1 = function(x1)
            var y2 = function(x2)
            var x3: Double
            var y3: Double

            var success = false
            var root = x1
            var i = 0

            if (y1 * y2 < 0.0) {
                do {
                    x3 = x2 - y2 / ((y2 - y1) / (x2 - x1))
                    y3 = function(x3)

                    if (y3 * y2 <= 0.0) {
                        x1 = x2; y1 = y2
                        x2 = x3; y2 = y3
                    } else {
                        y1 = y1 * y2 / (y2 + y3)
                        x2 = x3; y2 = y3
                    }

                    root = if (Math.abs(y1) < Math.abs(y2)) x1 else x2

                    success = Math.abs(x2 - x1) <= accuracy
                    i++
                } while (!success && (i < maxSteps))
            }

            return if (success) Result.Root(root) else Result.None


        }
    }


}
