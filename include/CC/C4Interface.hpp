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

#ifndef C4INTERFACE_HPP_
#define C4INTERFACE_HPP_


#include<vector>
#include "Tensor/FullTensor.hpp"
//#include"MatrixOperators.hpp"

namespace TensorCalculus {

  void readOccN(int &occA, int &occB, int &virtA, int &virtB);
  void readAO(const int aoCount, std::vector<double> &fullao); //C4: AOs symmetrieeffizient in TWOINT gespeichert --> 8-fach permutieren
  void readAO(std::vector<double> &fullao);
  void readAO(std::vector<double> &fullao, const bool chemical);
  void readAOchemical(const int aoCount, std::vector<double> &fullao); // chemical notation [i1 i2 | i3 i4], int(phi_1(x) phi_2(x) phi_3(y) phi_4(y) 1/||x-y||)dxdy
  void readAOchemical(std::vector<double> &fullao);
  void readTransformationMatrix(const int virt, const int occ, std::vector<double> &virtAOMO, std::vector<double> &occAOMO); //C4: NEWMOS
  void readTransformationMatrix(const int aoCount, std::vector<double> &fullTransposeAOMO); //C4: NEWMOS
  FullTensor<double> makeFullTensorAO (std::vector<double> &fullao);

  FullTensor<double> Vabij2Viabj ();
  FullTensor<double> Vabij2Viabj (const char * in);

  void readOrbitalEnergies (std::vector<double> &f); //Diagonal element of Fock matrix in MO basis C4: EPSILON

} /*end namespace TensorCalculus */

#endif /* C4INTERFACE_HPP_ */
