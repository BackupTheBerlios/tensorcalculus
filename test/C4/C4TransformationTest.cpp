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
#include <C4Interface.hpp>
#include <AOMOtransformation.hpp>
#include <Vector/VectorOperators.hpp>

using namespace TensorCalculus;
using namespace VectorOperators;

void MPSao2mo (const MPSRepresentation<double> input,const MPSRepresentation<double> mpsVirt, const MPSRepresentation<double> mpsOcc , std::vector<char> occState, const double eps, const char *out){
   std::cout<<"AO2MO for MPS ... ";
   std::ifstream f(&out[0]);
   if(!(f.fail())){
      std::cout<<"output file "<< &out[0]<<" already exist"<<std::endl;
      return;
   }

   bool fullTransformation = false; // only in MPSphys possible

   int occ, virt;
   readOccN(occ, occ, virt, virt);
   std::vector<int> compDimV_abij (4);

   for(int mu=0; mu < 4; mu++){
      if (occState[mu] == 'v'){
         compDimV_abij[mu] = virt;
      }else if (occState[mu] == 'o'){
         compDimV_abij[mu] = occ;
      }else{
         std::cout<<"undefined occumpation State at mu = "<<mu<<std::endl;
         return;
      }
   }
   MPSRepresentation<double> mpsV_abij (input.getSummations(), compDimV_abij);

   makeV(occState, mpsV_abij, mpsVirt, mpsOcc);
   std::cout<<"componentDimensions : "<<mpsV_abij.getComponentDimensions()<<std::endl;
   std::cout<<"summations old : "<<mpsV_abij.getSummations()<<std::endl;
   std::vector<double> mps_old (mpsV_abij.evaluate());
   mpsV_abij.truncateMPS(eps, 'm', true);

   mpsV_abij.write2disk(&out[0]);

   std::cout<<"summations new : "<<mpsV_abij.getSummations()<<std::endl;
   std::vector<double> mps_new (mpsV_abij.evaluate());
   double dist = l2_norm(mps_new - mps_old);
   std::cout<< "dist old-new :"<<dist<<"   rel dist : "<<dist / l2_norm(mps_old)<<std::endl;


   if(fullTransformation){
      std::string filename = "v_abij.ften";
      std::ifstream fin2(&filename[0]);
      if(fin2.fail()){
         throw std::invalid_argument("no v_abij.ften file");
      }
      FullTensor<double> ftV_abij (&filename[0]);

      double normFT = l2_norm(ftV_abij.getV());
      std::cout<<" norm fullTensor : "<<normFT<<std::endl;
      std::cout<<" rel dist        : "<<l2_norm(mpsV_abij.evaluate()-ftV_abij.getV())<<std::endl;
   }
   /*
   //testing...
   std::vector<std::vector<double> > v(4);
   std::vector<int> summations(3);
   summations[0]=1;
   summations[1]=1;
   summations[2]=1;
   std::vector<int> size(4);
   size[0] = summations[0];
   size[1] = summations[0] * summations[1];
   size[2] = summations[1] * summations[2];
   size[3] = summations[2];
   for(int mu=0; mu<4; mu++){
      size[mu] *= mpsV_abij.getComponentDimension(mu);
      v[mu].resize(size[mu]);
      for(int i=0; i<size[mu]; i++)
         v[mu][i] = Utilities<double>::rand();
   }
   //rank 1 starttensor
   MPSRepresentation<double> mpsV_abij2 (summations, v, mpsV_abij.getComponentDimensions());
   //mpsV_abij.computeOrthonormalBasis(); //unterraumdarstellung;

   mpsV_abij2.performDMRG(mpsV_abij, 1e-12);
   mpsV_abij2.performDMRG(mpsV_abij, 1e-12);

   std::cout<<"summations old : "<<mpsV_abij.getSummations()<<std::endl;
   std::cout<<"summations new : "<<mpsV_abij2.getSummations()<<std::endl;

   //mpsV_abij2.write2disk(&out[0]);
    */
   std::cout<<"AO2MO for MPS ... finish"<<std::endl;
}

MPSRepresentation<double> input (const char *in){
   std::ifstream fin(&in[0]);
   if(fin.fail()){
      throw std::invalid_argument("no input file");
   }
    return(MPSRepresentation<double>(&in[0]));
}

int MPStransform (std::vector<char> occState, const char * inFile, std::vector<std::string> out, std::vector<double> eps){
   MPSRepresentation<double> in (input(&inFile[0]));

   int occ, virt;
   readOccN(occ, occ, virt, virt);
   std::vector<int> compDimVirt (4);
   std::vector<int> compDimOcc (4);

   for(int mu=0; mu < 4; mu++){
      compDimVirt[mu] = virt;
      compDimOcc[mu]  = occ;
   }

   MPSRepresentation<double> mpsVirt (in.getSummations(), compDimVirt);
   MPSRepresentation<double> mpsOcc (in.getSummations(), compDimOcc);
   AO2MO(in, mpsVirt, mpsOcc);

   for(int i=0; i<out.size(); i++){
      MPSao2mo(in, mpsVirt, mpsOcc, occState, eps[i], &out[i][0]);
   }

   return 0;
}



int MPStransform (std::vector<char> occState, const char * inFile, const char* out, double eps){
   std::vector<double> eps1(1);
   eps1[0] = eps;
   std::vector<std::string> out1(1);
   out1[0] = &out[0];
   MPStransform(occState, inFile, out1, eps1);
}

int main (){

   std::vector<char> occState(4);
   occState[0] = 'v';
   occState[1] = 'v';
   occState[2] = 'o';
   occState[3] = 'o';
   std::cout<<"tranform: phys, eps = 1e-12 ..."<<std::endl;
   MPStransform(occState,"compute/fullAOphys.mps.1e-12", "compute/AOphys2Vabij.mps.1e-12", 1e-12);
   std::cout<<"tranform: phys, eps = 1e-6 ..."<<std::endl;
   MPStransform(occState,"compute/fullAOphys.mps.1e-6", "compute/AOphys2Vabij.mps.1e-6", 1e-6);
   std::cout<<"tranform: phys, eps = 1e-4 ..."<<std::endl;
   MPStransform(occState,"compute/fullAOphys.mps.1e-4", "compute/AOphys2Vabij.mps.1e-4", 1e-4);
   std::cout<<"tranform: phys, eps = 1e-2 ..."<<std::endl;
   MPStransform(occState,"compute/fullAOphys.mps.1e-2", "compute/AOphys2Vabij.mps.1e-2", 1e-2);

   occState[0] = 'v';
   occState[1] = 'o';
   occState[2] = 'v';
   occState[3] = 'o';

   std::cout<<"tranform: chem, eps = 1e-12 ..."<<std::endl;
   MPStransform(occState,"compute/fullAOchem.mps.1e-12", "compute/AOchem2Vaibj.mps.1e-12", 1e-12);
   std::cout<<"tranform: chem, eps = 1e-6 ..."<<std::endl;
   MPStransform(occState,"compute/fullAOchem.mps.1e-6", "compute/AOchem2Vaibj.mps.1e-6", 1e-6);
   std::cout<<"tranform: chem, eps = 1e-4 ..."<<std::endl;
   MPStransform(occState,"compute/fullAOchem.mps.1e-4", "compute/AOchem2Vaibj.mps.1e-4", 1e-4);
   std::cout<<"tranform: chem, eps = 1e-2 ..."<<std::endl;
   MPStransform(occState,"compute/fullAOchem.mps.1e-2", "compute/AOchem2Vaibj.mps.1e-2", 1e-2);

}

