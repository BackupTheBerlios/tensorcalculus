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

#ifndef AMPLITUDES_HPP_
#define AMPLITUDES_HPP_

#include<vector>
#include "Representation/CPTensorRepresentation.hpp"
#include "DKTS.hpp"
#include "Representation/MPSRepresentation.hpp"

namespace TensorCalculus{

CPTensorRepresentation<double> readEpsilonInvDKTS ();
CPTensorRepresentation<double> makeEpsilonInv(const std::vector<char> &occState);
void V2T (MPSRepresentation<double> &V_mnop);
void V2T (MPSRepresentation<double> &V_mnop,const std::vector<char> &occState);


class Amplitudes{
private:
   int occ, virt, aoCount;
   std::vector<double> f;

   std::vector<double> fTensor;
   std::vector<double> f_vv_oo;
   std::vector<double> v_vv_oo; //spinadapted : normal computation: abij; chemical computation aibj;
   std::vector<double> v_vv_oo_simple;
   std::vector<double> t_vv_oo; // normal computation: abij; chemical computation aibj;

   //CPTensorRepresentation<double> fCP;

   void init (const std::vector<double> &spinAdaptV, const std::vector<double> &v);
   void makeOrbitalEnergieTensor ();
   void makeOrbitalEnergieCPTensor ();
   void makeF_vv_oo ();
   //void makeOrbitalEnergieTensorChemical ();
   void makeT_vv_oo ();
   //void makeTchemical ();

   //chemical
   //todo
public:
   Amplitudes (const std::vector<double> &spinAdaptV, const std::vector<double> &v){
      init(spinAdaptV, v);
   }
   //void computeT ();
   //void computeTchemical ();
   double MP2Energie ();
}; // end class Amplitudes
}// end namespace TensorCalculus

#endif /* AMPLITUDES_HPP_ */
