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

#ifndef MP2_HPP_
#define MP2_HPP_

#include "Representation/MPSRepresentation.hpp"

namespace TensorCalculus {


   double mp2_mps (const double eps, const char *out, const bool chemical = false); // using TWOINT, EPSILON, ...
   double mp2_mps_canonical (const double eps, const char *out, const bool chemical = false); // using v_abij.ften, NOCC, EPSILON

   double mp2_mps (const char *v_abij_mps, const char *t_abij_mps, const bool chemical = false); // using v and t in MPS format
   double mp2_mps (MPSRepresentation<double> &v_abij, MPSRepresentation<double> &t_abij, const bool chemical = false); //computes all contractions
   double mp2_mps_phys (MPSRepresentation<double> &v_abij, MPSRepresentation<double> &t_abij); //computes all contractions
   double mp2_mps_chem (MPSRepresentation<double> &v_abij, MPSRepresentation<double> &t_abij);

   double mp2_ft (const char *v_abij, const char *t_abij); // using v and t form C4

}






#endif /* MP2_HPP_ */
