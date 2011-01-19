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
#include <iostream>
#include <Utilities/Utilities.hpp>
#include <Tensor/FullTensor.hpp>
#include <Representation/MPSRepresentation.hpp>
#include <CC/C4Interface.hpp>
#include <Vector/VectorOperators.hpp>

#include <CC/Amplitudes.hpp>

using namespace TensorCalculus;
using namespace VectorOperators;

void checkEpsInv ();
void checkT (const char *inFile, const std::vector<char> occState);
void checkT (const char *inFile, const char occState1, const char occState2, const char occState3, const char occState4){
   std::vector<char> occState(4);
   occState[0] = occState1;
   occState[1] = occState2;
   occState[2] = occState3;
   occState[3] = occState4;
   checkT (inFile, occState);
}


int main () {
   checkEpsInv();
   checkT("compute/phys_Vabij.mps.1e-12_1e-12", 'v', 'v', 'o', 'o');
   return 0;
}

void checkEpsInv (){
   std::cout<<"checkEpsInv ... ";
   CPTensorRepresentation<double> EpsInv(readEpsilonInvDKTS ());
   std::vector<double> f;
   readOrbitalEnergies(f);
   int occ, virt;
   readOccN(occ, occ, virt, virt);
   int maxOccN = std::max(occ, virt);

   std::vector<double> fOcc(maxOccN);
   std::vector<double> fVirt(maxOccN);

   for(int i=0; i < occ; i++){
      fOcc[i] = f[i];
   }
   for(int i=0; i < virt; i++){
      fVirt[i] = f[occ+i];
   }
   //std::vector<double> fullEpsInv (maxOccN*maxOccN*maxOccN*maxOccN);
   std::vector<double> fullEpsInv (virt*virt*occ*occ);
   for(int i=0; i< occ; i++){
      for(int j=0; j< occ; j++){
         for(int k=0; k< virt; k++){
            for(int l=0; l< virt; l++){
               if(i<occ && j<occ && k<virt && l<virt){
                  fullEpsInv[l + virt * (k + virt * (j + occ * i))]
                          = 1/(-1.0 * fVirt[l] - fVirt[k] + fOcc[j] + fOcc[i]);
               }
            }
         }
      }
   }
   /*
   std::cout<<EpsInv.evaluate()<<std::endl;
   std::cout<<"#############################################"<<std::endl;
   std::cout<<fullEpsInv<<std::endl;
   std::cout<<"size1 : "<<fullEpsInv.size()<<"size2 : "<<EpsInv.evaluate().size()<<std::endl;
   */
   std::cout<<"rel dist : "<<l2_norm(fullEpsInv - EpsInv.evaluate())/l2_norm(fullEpsInv)<<std::endl;
   std::cout<<"finish"<<std::endl;
}

void checkT (const char *inFile, const std::vector<char> occState){
   std::cout<<"checkT ... ";
   std::ifstream f (&inFile[0]);
   if(f.fail()){
      std::cout<<"no input file"<<std::endl;
      return;
   }
   MPSRepresentation<double> in(&inFile[0]);
   V2T(in);
   std::ifstream f2 (&inFile[0]);
   if(f2.fail()){
      std::cout<<"no t_abij.ten file"<<std::endl;
      return;
   }
   FullTensor<double> fullT("t_abij.ten");
   std::cout<< "rel dist : "<<l2_norm(fullT.getV()-in.evaluate()) / l2_norm(fullT.getV())<<std::endl;
   std::cout<<"finish"<<std::endl;
}
