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

using namespace TensorCalculus;
using namespace VectorOperators;

void prepareFullAO (std::vector<int> &componentDimensions, std::vector <double> &v);
void writeFullAO ( const std::vector<int> &componentDimensions, const std::vector <double> &v, const char *ft);
void makeFullAOphys (const char *ft);
void makeFullAOchem (const char *ft);

void DecomposeFullAO2MPS (const char *in, const char *out, const double eps);
void DecomposeFullAO2MPS (const FullTensor<double> in, const char *out, const double eps);
void DecomposeFullAO2CC (const char *in, const char * out, const int rank1, const int rank2, double alsEps);
void DecomposeFullAO2CC (const char *in, const char * out, const int rank, double alsEps = 1e-6){
   DecomposeFullAO2CC(in,out,rank,rank, alsEps);
}

TensorRepresentation<double> CCstartValue (int rank1, int rank2, MPSRepresentation<double> &in);

void checkNorm (const char *ft, const char *cc);
void checkNorm (const char *ft, const char *mps, const char *cc);


int main (){
    int aoCount, temp;
    readOccN(aoCount,aoCount, temp,temp);
    aoCount += temp;
    makeFullAOphys ("compute/fullAOphys.ft");
    makeFullAOchem ("compute/fullAOchem.ft");

    FullTensor<double> fullAOphys("compute/fullAOphys.ft");
    DecomposeFullAO2MPS (fullAOphys, "compute/fullAOphys.mps.1e-12", 1e-12);
    DecomposeFullAO2MPS (fullAOphys, "compute/fullAOphys.mps.1e-6", 1e-6);
    DecomposeFullAO2MPS (fullAOphys, "compute/fullAOphys.mps.1e-4", 1e-4);
    DecomposeFullAO2MPS (fullAOphys, "compute/fullAOphys.mps.1e-2", 1e-2);


    FullTensor<double> fullAOchem("compute/fullAOchem.ft");
    DecomposeFullAO2MPS (fullAOchem, "compute/fullAOchem.mps.1e-12", 1e-12);
    DecomposeFullAO2MPS (fullAOchem, "compute/fullAOchem.mps.1e-6", 1e-6);
    DecomposeFullAO2MPS (fullAOchem, "compute/fullAOchem.mps.1e-4", 1e-4);
    DecomposeFullAO2MPS (fullAOchem, "compute/fullAOchem.mps.1e-2", 1e-2);

   // DecomposeFullAO2CC ("compute/fullAOchem.mps.1e-12", "compute/fullAOchem.cc", 4*aoCount);

   // checkNorm("compute/fullAOchem.ft","compute/fullAOchem.mps.1e-12","compute/fullAOchem.cc");
}

void makeFullAOphys (const char *ft){
   std::cout<<"makeFullAOphys ..."<<std::endl;
   std::ifstream fin(ft);
   if( !(fin.fail()) ){
      std::cout<<ft<<" already exists"<<std::endl;
      return;
   }
   std::vector<int> compDim;
   std::vector<double> v;
   prepareFullAO(compDim, v);
   readAO(compDim[0], v);
   writeFullAO(compDim, v, ft);
   std::cout<<"makeFullAOphys ... finish"<<std::endl;
}

void makeFullAOchem (const char *ft){
   std::cout<<"makeFullAOchem ..."<<std::endl;
   std::ifstream fin(ft);
   if( !(fin.fail()) ){
      std::cout<<ft<<" already exists"<<std::endl;
      return;
   }
   std::vector<int> compDim;
   std::vector<double> v;
   prepareFullAO(compDim, v);
   readAOchemical(compDim[0], v);
   writeFullAO(compDim, v, ft);
   std::cout<<"makeFullAOchem ... finish"<<std::endl;
}

void prepareFullAO (std::vector<int> &componentDimensions, std::vector <double> &v){
   int virt, occ;
   readOccN(occ, occ, virt, virt);
   int aoCount = virt + occ;
   componentDimensions.resize(4);
   componentDimensions[0] = aoCount;
   componentDimensions[1] = aoCount;
   componentDimensions[2] = aoCount;
   componentDimensions[3] = aoCount;
   v.resize(aoCount * aoCount * aoCount * aoCount);
}

void writeFullAO ( const std::vector<int> &componentDimensions, const std::vector <double> &v, const char *ft){
   FullTensor<double> fullTensor(componentDimensions, v);
   fullTensor.write2disk(ft);
}

void DecomposeFullAO2MPS (const char *in, const char *out, const double eps){
   std::cout<<"DecomposeFullAO2MPS "<< in <<" ..."<<std::endl;
   std::ifstream fin(in);
   if(fin.fail()){
      throw std::invalid_argument("no inputfile");
   }
   FullTensor<double> ft(in);
   DecomposeFullAO2MPS (ft, out, eps);
}
void DecomposeFullAO2MPS (const FullTensor<double> in, const char *out, const double eps){
   std::cout<<"DecomposeFullAO2MPS "<<out <<" ..."<<std::endl;
   std::ifstream fin2(out);
   if(!(fin2.fail())){
      std::cout<<"output file already exists"<<std::endl;
      return;
   }
   MPSRepresentation<double> mps(in, eps);
   mps.write2disk(out);
   std::cout<<"DecomposeFullAO2MPS "<<out <<" ... finish"<<std::endl;
}

void DecomposeFullAO2CC (const char *in, const char * out, const int rank1, const int rank2, double alsEps){
   std::cout<<"DecomposeFullAO2CC "<<out <<" ..."<<std::endl;
   std::ifstream fin(in);
   if(fin.fail()){
      throw std::invalid_argument("no inputfile");
   }
   std::ifstream fin2(out);
   if(!(fin2.fail())){
      std::cout<<"output file already exists"<<std::endl;
      return;
   }
   MPSRepresentation<double> mps;
   mps.read(in);
   TensorRepresentation<double> cc(CCstartValue(rank1, rank2, mps));

   cc.performALS(mps, alsEps);
   cc.writeToFile(out);
   std::cout<<"DecomposeFullAO2CC "<<out <<" ... finish"<<std::endl;
}

TensorRepresentation<double> CCstartValue (int rank1, int rank2, MPSRepresentation<double> &in){

   std::vector<std::vector<double> > v(4);
   v[0].resize(in.getComponentDimension(0)*rank1);
   v[1].resize(in.getComponentDimension(0)*rank1);
   v[2].resize(in.getComponentDimension(0)*rank2);
   v[3].resize(in.getComponentDimension(0)*rank2);
   for(int mu=0; mu<2; mu++){
      for(int i=0; i<in.getComponentDimension(0)*rank1; i++){
         v[mu][i] = Utilities<double>::rand();
      }
      for(int i=0; i<in.getComponentDimension(0)*rank2; i++){
         v[mu+2][i] = Utilities<double>::rand();
      }
   }

   std::vector< std::vector<double> > w(1);
   w[0].resize(rank1*rank2);
   for(int i=0; i<rank1 * rank2; i++){
      w[0][i] = Utilities<double>::rand();
   }

   std::vector<int> summations (2);
   summations[0] = rank1;
   summations[1] = rank2;

   std::vector<std::vector<int> > incidenceMatrix(5);
   for(int mu=0; mu<4; mu++){
      incidenceMatrix[mu].resize(1);
   }
   incidenceMatrix[2][0] = 1;
   incidenceMatrix[3][0] = 1;
   incidenceMatrix[4].resize(2);
   incidenceMatrix[4][1] = 1;

   return(TensorRepresentation<double> (summations, v, in.getComponentDimensions(), incidenceMatrix, w));
}

void checkNorm (const char *ft, const char *cc){
   std::cout<<"check Norm "<<cc <<" ..."<<std::endl;
   std::ifstream fin1(cc);
   std::ifstream fin2(ft);
   if(fin1.fail()|| fin2.fail()){
      throw std::invalid_argument("no inputfile");
   }
   FullTensor<double> fTensor(ft);
   TensorRepresentation<double> ccTensor(cc);
   double normFt = l2_norm(fTensor.getV());
   double dist   = l2_norm(fTensor.getV()-ccTensor.evaluate());
   std::cout<<"norm fullAO = "<<normFt<<std::endl;
   std::cout<<"norm(fullAO-CC)/norm(fullAO) = "<<dist/normFt<<std::endl;
}

void checkNorm (const char *ft, const char *mps, const char *cc){
   MPSRepresentation<double> mpsT;
   FullTensor<double> fullT(ft);
   TensorRepresentation<double> ccT(cc);
   mpsT.read(mps);
   std::vector<double> ccV(ccT.evaluate());
   std::vector<double> mpsV(mpsT.evaluate());
   std::vector<double> fullV(fullT.getV());
   std::cout<<"norm fullAO = "<<l2_norm(fullV)<<std::endl;
   std::cout<<"relNorm (full-cc)"<<l2_norm(fullV-ccV)<<std::endl;
   std::cout<<"relNorm (full-mps)"<<l2_norm(fullV-mpsV)<<std::endl;
}
