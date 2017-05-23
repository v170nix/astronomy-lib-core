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

package net.arwix.astronomy.physical


interface PhysicalProperties {
    /**
     * The equatorial radius in km.
     */
    val equatorialRadius: Double
    /**
     * The polar radius in km.
     */
    val polarRadius: Double

    val massRatio: Double
    /**
     * Returns flattening factor = (equatorial radius - polar radius ) /
     * equatorial radius.
     *
     * @return Flattening factor. Set to 0 if the object size is unknown.
     */
    fun getFlatteningFactor() = if (equatorialRadius != 0.0) (equatorialRadius - polarRadius) / equatorialRadius else 0.0
}