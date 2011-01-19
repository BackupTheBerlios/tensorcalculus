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

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <Utilities/Utilities.hpp>
#include <Tensor/FullTensor.hpp>
#include <Representation/MPSRepresentation.hpp>
#include <CC/C4Interface.hpp>
#include <CC/AOMOtransformation.hpp>
#include <Vector/VectorOperators.hpp>

using namespace TensorCalculus;

void MPSTransformation (std::vector<char> occState, const char* AO,
                        const char* MO_vv_oo,
                        const char* MO_vv_vv,
                        const char* MO_oo_oo, const double eps);


int main(){

   std::vector<char> occState(4);
   occState[0] = 'v';
   occState[1] = 'v';
   occState[2] = 'o';
   occState[3] = 'o';

   MPSTransformation (occState, "compute/fullAOphys.mps.1e-12",
            "compute/AOphys2Vabij.mps.1e-12_new", "compute/AOphys2Vabcd.mps.1e-12",
            "compute/AOphys2Vijkl.mps.1e-12",1e-12);
   MPSTransformation (occState, "compute/fullAOphys.mps.1e-6",
            "compute/AOphys2Vabij.mps.1e-6_new", "compute/AOphys2Vabcd.mps.1e-6",
            "compute/AOphys2Vijkl.mps.1e-6",1e-6);
   MPSTransformation (occState, "compute/fullAOphys.mps.1e-4",
            "compute/AOphys2Vabij.mps.1e-4_new", "compute/AOphys2Vabcd.mps.1e-4",
            "compute/AOphys2Vijkl.mps.1e-4",1e-4);
   MPSTransformation (occState, "compute/fullAOphys.mps.1e-2",
            "compute/AOphys2Vabij.mps.1e-2_new", "compute/AOphys2Vabcd.mps.1e-2",
            "compute/AOphys2Vijkl.mps.1e-2",1e-2);

   occState[0] = 'v';
   occState[1] = 'o';
   occState[2] = 'v';
   occState[3] = 'o';

   MPSTransformation (occState, "compute/fullAOchem.mps.1e-12",
            "compute/AOchem2Vaibj.mps.1e-12_new", "compute/AOchem2Vacbd.mps.1e-12",
            "compute/AOchem2Vikjl.mps.1e-12",1e-12);
   MPSTransformation (occState, "compute/fullAOchem.mps.1e-6",
               "compute/AOchem2Vaibj.mps.1e-6_new", "compute/AOchem2Vacbd.mps.1e-6",
               "compute/AOchem2Vikjl.mps.1e-6",1e-6);
   MPSTransformation (occState, "compute/fullAOchem.mps.1e-4",
               "compute/AOchem2Vaibj.mps.1e-4_new", "compute/AOchem2Vacbd.mps.1e-4",
               "compute/AOchem2Vikjl.mps.1e-4",1e-4);
   MPSTransformation (occState, "compute/fullAOchem.mps.1e-2",
               "compute/AOchem2Vaibj.mps.1e-2_new", "compute/AOchem2Vacbd.mps.1e-2",
               "compute/AOchem2Vikjl.mps.1e-2",1e-2);


   return 0;
}

void MPSTransformation (std::vector<char> occState, const char* AO,
                        const char* MO_vv_oo,
                        const char* MO_vv_vv,
                        const char* MO_oo_oo, const double eps){
   int occ, virt;
   readOccN(occ, occ, virt, virt);

   std::ifstream fin(&AO[0]);
   if(fin.fail()){
      throw std::invalid_argument("no input file");
   }

   std::ifstream fin2 (&MO_vv_oo[0]);
   if(!fin2.fail()){
         std::cout<<"output already exist"<<std::endl;
         return;
   }

   MPSRepresentation<double> aoMPS(&AO[0]);

   std::vector<int> compDimVirt(4);
   std::vector<int> compDimMO(4);
   std::vector<int> compDimOcc(4);
   compDimVirt[0] = virt;
   compDimVirt[1] = virt;
   compDimVirt[2] = virt;
   compDimVirt[3] = virt;

   compDimOcc[0] = occ;
   compDimOcc[1] = occ;
   compDimOcc[2] = occ;
   compDimOcc[3] = occ;

   for(int mu=0; mu < 4; mu++){
      if (occState[mu] == 'v'){
         compDimMO[mu] = virt;
      }else if (occState[mu] == 'o'){
         compDimMO[mu] = occ;
      }else{
         std::cout<<"undefined occumpation State at mu = "<<mu<<std::endl;
         return;
      }
   }

   MPSRepresentation<double> occMPS (aoMPS.getSummations(), compDimOcc);
   MPSRepresentation<double> virtMPS (aoMPS.getSummations(), compDimVirt);
   MPSRepresentation<double> moMPS (aoMPS.getSummations(), compDimMO);

   AO2MO (aoMPS, virtMPS, occMPS);
   makeV(occState, moMPS, virtMPS, occMPS);

   moMPS.truncateMPS(eps, 'm', true);
   moMPS.write2disk(&MO_vv_oo[0]);

   virtMPS.truncateMPS(eps, 'm', true);
   virtMPS.write2disk(&MO_vv_vv[0]);

   occMPS.truncateMPS(eps, 'm', true);
   occMPS.write2disk(&MO_oo_oo[0]);
}




