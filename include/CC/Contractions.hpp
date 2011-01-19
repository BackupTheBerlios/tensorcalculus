/*
 * Copyright (C) 2011 Henry Auer
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CONTRACTIONS_HPP_
#define CONTRACTIONS_HPP_

#include "Representation/MPSRepresentation.hpp"

namespace TensorCalculus{
   //only for d=4
   MPSRepresentation<double> contract_21_43 (MPSRepresentation<double> &T,
                                             MPSRepresentation<double> &V);
   double contract_11_24_33_42 (MPSRepresentation<double> &T, MPSRepresentation<double> &V);
   double contract_14_22_33_41 (MPSRepresentation<double> &T, MPSRepresentation<double> &V);
}

#endif /* CONTRACTIONS_HPP_ */
