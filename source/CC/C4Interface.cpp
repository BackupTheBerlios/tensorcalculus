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

#include "CC/C4Interface.hpp"
#include <iostream>
#include <fstream>
//#include "VectorOperators.hpp"

namespace TensorCalculus{

void readOccN(int &occA, int &occB, int &virtA, int &virtB){
   std::ifstream fin("NOCC");
   //fin.open("TWOINT");
   if(fin.fail()){
      throw std::invalid_argument("No NOCC file");
   }
   fin >> occA;
   fin >> occB;
   fin >> virtA;
   fin >> virtB;
}

void readAO(std::vector<double> &fullao, const bool chemical){
   if(chemical){
      readAOchemical(fullao);
   } else {
      readAO(fullao);
   }
}

 void readAO(const int aoCount, std::vector<double> &fullao){
   //tested

  std::ifstream fin("TWOINT");
  //fin.open("TWOINT");
  if(fin.fail()){
     throw std::invalid_argument("No TWOINT file");
  }
  fin.setf(std::ios::scientific, std::ios::floatfield);
  int i1 = 0;
  int i2 = 0;
  int i3 = 0;
  int i4 = 0;
  double value;
  unsigned int aoCount2 = aoCount *aoCount;
  unsigned int aoCount3 = aoCount2*aoCount;
  unsigned int aoCount4 = aoCount3*aoCount;

  if(fullao.size()<aoCount4){
     fullao.resize(aoCount4);
  }
  std::vector<int> index(8);


  // der index zählt von rechts nach links hoch: (1,1,1,1) (2,1,1,1) ... (1,2,1,1) (2,2,1,1) ...
  //todo: switchanweisungen für paarweise gleiche

  //todo forschleife für alle, wie viele???
  int i = 0;
  int sum = 0;
  while(sum < 4*aoCount){
     i++;
     fin>> i1 >> i2 >> i3 >> i4 >> value;
     //std::cout<<i1<<"  "<<i2<<"  "<<i3<<"  "<<i4<<" : "<<i<<" = "<<value<<"   ||   ";
     sum = i1+i2+i3+i4;
     i1--;
     i2--;
     i3--;
     i4--;

     /*
     // nach dem Code von Udo --> chemische notation with permutation [i1 i3 | i2 i4]
     index[0] = i1+aoCount*i3+aoCount2*i2+aoCount3*i4;
     index[1] = i3+aoCount*i1+aoCount2*i2+aoCount3*i4;
     index[2] = i1+aoCount*i3+aoCount2*i4+aoCount3*i2;
     index[3] = i3+aoCount*i1+aoCount2*i4+aoCount3*i2;

     index[4] = i2+aoCount*i4+aoCount2*i1+aoCount3*i3;
     index[5] = i2+aoCount*i4+aoCount2*i3+aoCount3*i1;
     index[6] = i4+aoCount*i2+aoCount2*i1+aoCount3*i3;
     index[7] = i4+aoCount*i2+aoCount2*i3+aoCount3*i1;
      */
     // physical notation with permutation of 2nd and 3rd direction
     //<i1 i3 | i2 i4>
       index[0] = i1+aoCount*i3+aoCount2*i2+aoCount3*i4;
       index[1] = i2+aoCount*i3+aoCount2*i1+aoCount3*i4;
       index[2] = i1+aoCount*i4+aoCount2*i2+aoCount3*i3;
       index[3] = i2+aoCount*i4+aoCount2*i1+aoCount3*i3;

       index[4] = i3+aoCount*i2+aoCount2*i4+aoCount3*i1;
       index[5] = i3+aoCount*i1+aoCount2*i4+aoCount3*i2;
       index[6] = i4+aoCount*i2+aoCount2*i3+aoCount3*i1;
       index[7] = i4+aoCount*i1+aoCount2*i3+aoCount3*i2;


     /*
     // <i1 i2 | i3 i4 > physical notation
     index[0] = i1+aoCount*i2+aoCount2*i3+aoCount3*i4;
     index[1] = i1+aoCount*i4+aoCount2*i3+aoCount3*i2;
     index[2] = i3+aoCount*i2+aoCount2*i1+aoCount3*i4;
     index[3] = i3+aoCount*i4+aoCount2*i1+aoCount3*i2;

     index[4] = i2+aoCount*i1+aoCount2*i4+aoCount3*i3;
     index[5] = i2+aoCount*i3+aoCount2*i4+aoCount3*i1;
     index[6] = i4+aoCount*i1+aoCount2*i2+aoCount3*i3;
     index[7] = i4+aoCount*i3+aoCount2*i2+aoCount3*i1;
     */
     /*
      // [i1 i2 | i3 i4] chemical notation
     index[0] = i1+aoCount*i2+aoCount2*i3+aoCount3*i4;
     index[1] = i2+aoCount*i1+aoCount2*i3+aoCount3*i4;
     index[2] = i1+aoCount*i2+aoCount2*i4+aoCount3*i3;
     index[3] = i2+aoCount*i1+aoCount2*i4+aoCount3*i3;

     index[4] = i3+aoCount*i4+aoCount2*i1+aoCount3*i2;
     index[5] = i4+aoCount*i3+aoCount2*i1+aoCount3*i2;
     index[6] = i3+aoCount*i4+aoCount2*i2+aoCount3*i1;
     index[7] = i4+aoCount*i3+aoCount2*i2+aoCount3*i1;
      */

     for(int i=0; i<8; i++){
          fullao[index[i]]   =  value;
          //fullao[index[i+4]] = -1.0*value;
     }
  }
 }; /* readAO */

 void readAO(std::vector<double> &fullao){
    int occA, occB, virtA, virtB;
    readOccN(occA, occB, virtA, virtB);
    int aoCount = occA + virtA;
    readAO(aoCount, fullao);
 }; // end readAO(std::vector<double> &fullao)

 void readAOchemical(const int aoCount, std::vector<double> &fullao) {
    std::ifstream fin("TWOINT");
   //fin.open("TWOINT");
   if(fin.fail()){
      throw std::invalid_argument("No TWOINT file");
   }
   fin.setf(std::ios::scientific, std::ios::floatfield);
   int i1 = 0;
   int i2 = 0;
   int i3 = 0;
   int i4 = 0;
   double value;
   unsigned int aoCount2 = aoCount *aoCount;
   unsigned int aoCount3 = aoCount2*aoCount;
   unsigned int aoCount4 = aoCount3*aoCount;

   if(fullao.size()<aoCount4){
      fullao.resize(aoCount4);
   }
   std::vector<int> index(8);


   // der index zählt von rechts nach links hoch: (1,1,1,1) (2,1,1,1) ... (1,2,1,1) (2,2,1,1) ...
   //todo: switchanweisungen für paarweise gleiche

   //todo forschleife für alle, wie viele???
   int i = 0;
   int sum = 0;
   while(sum < 4*aoCount){
      i++;
      fin>> i1 >> i2 >> i3 >> i4 >> value;
      //std::cout<<i1<<"  "<<i2<<"  "<<i3<<"  "<<i4<<" : "<<i<<" = "<<value<<"   ||   ";
      sum = i1+i2+i3+i4;
      i1--;
      i2--;
      i3--;
      i4--;


       // [i1 i2 | i3 i4] chemical notation
      index[0] = i1+aoCount*i2+aoCount2*i3+aoCount3*i4;
      index[1] = i2+aoCount*i1+aoCount2*i3+aoCount3*i4;
      index[2] = i1+aoCount*i2+aoCount2*i4+aoCount3*i3;
      index[3] = i2+aoCount*i1+aoCount2*i4+aoCount3*i3;

      index[4] = i3+aoCount*i4+aoCount2*i1+aoCount3*i2;
      index[5] = i4+aoCount*i3+aoCount2*i1+aoCount3*i2;
      index[6] = i3+aoCount*i4+aoCount2*i2+aoCount3*i1;
      index[7] = i4+aoCount*i3+aoCount2*i2+aoCount3*i1;

      for(int i=0; i<8; i++){
           fullao[index[i]]   =  value;
      }
   }
 }; /* readAOchemical */

 void readAOchemical(std::vector<double> &fullao){
     int occA, occB, virtA, virtB;
     readOccN(occA, occB, virtA, virtB);
     int aoCount = occA + virtA;
     readAOchemical(aoCount, fullao);
  }; // end readAOchemical(std::vector<double> &fullao)


/*
 void readAO(const int aoCount, TensorCalculus::FullTensor<double> &fullao){
    std::vector<int> readIndex(4);
    std::vector<int> index(4);
    double value;

    std::ifstream fin("TWOINT");
    //fin.open("TWOINT");
    if(fin.fail()){
       std::cout<<"No TWOINT file"<<std::endl;
    }
    fin.setf(std::ios::scientific, std::ios::floatfield);

    while(fin >> readIndex[0]){
       fin >> readIndex[1];
       fin >> readIndex[2];
       fin >> readIndex[3];
       fin >> value;

    }

 }; // readAO
*/

 void readTransformationMatrix(const int virt,const int occ, std::vector<double> &virtAOMO, std::vector<double> &occAOMO){
    // tested
   const unsigned int ao = virt+occ;

   std::ifstream fin("NEWMOS");
   if(fin.fail()){
      std::cout<<"No NEWMOS file"<<std::endl;
      return;
   }
   fin.setf(std::ios::scientific, std::ios::floatfield);

   if(virtAOMO.size() < ao*virt){
      virtAOMO.resize(ao*virt);
   }
   if(occAOMO.size() < ao*occ){
      occAOMO.resize(ao*occ);
   }

   //NEWMOS = M^T Matrix(ao;mo) : M * AO = MO
   // occAOMO   in Matrix(occ;ao)
   // virtAOMO  in Matrix(virt;ao);
   for(int i=0; i<occ; i++){
      for(unsigned int j=0; j<ao; j++){
         fin >> occAOMO[i+j*occ];
      }
   }

   for(int i=0; i<virt; i++){
      for(unsigned int j=0; j<ao; j++){
         fin >> virtAOMO[i+j*virt];
      }
   }
 }; /*end readTransformationMatrix */

 void readTransformationMatrix(const int aoCount, std::vector<double> &fullTransposeAOMO){
    //tested
    std::ifstream fin("NEWMOS");
    if(fin.fail()){
       std::cout<<"No NEWMOS file"<<std::endl;
       return;
    }
    fin.setf(std::ios::scientific, std::ios::floatfield);

    if(fullTransposeAOMO.size() < static_cast<unsigned int>(aoCount*aoCount)){
       fullTransposeAOMO.resize(aoCount* aoCount);
    }

    //NEWMOS = M^T  Matrix(ao;mo) : M * AO = MO
    // fullTransposeAOMO = M^T
    for(int i=0; i<aoCount; i++){
       for(int j=0; j<aoCount; j++){
          fin >> fullTransposeAOMO[j+i*aoCount];
       }
    }
  }; /*end readTransformationMatrix */

 void readOrbitalEnergies (std::vector<double> &f){
    int virt,occ;
    readOccN(occ,occ,virt,virt);
    unsigned int aoCount = virt+occ;
    if(f.size()<aoCount){
       f.resize(aoCount);
    }
    std::ifstream fin("EPSILON");
    if(fin.fail()){
       std::cout<<"No EPSILON file"<<std::endl;
       return;
    }
    double temp;
    for(unsigned int i=0; i<aoCount; i++){
       fin >> temp;
       f[i] = std::fabs(temp);
    }
 }//end readOrbitalEnergies (std::vector<double> &f)

 FullTensor<double> makeFullTensorAO (std::vector<double> &fullao){
    int occ, virt;
    readOccN(occ, occ, virt, virt);
    int aoCount = occ + virt;
    std::vector<int> componentDimensions(4);
    componentDimensions[0] = aoCount;
    componentDimensions[1] = aoCount;
    componentDimensions[2] = aoCount;
    componentDimensions[3] = aoCount;

    return(FullTensor<double> (componentDimensions, fullao));
 }
 FullTensor<double> Vabij2Viabj (){
    return Vabij2Viabj("v_abij.ften");
 }
 FullTensor<double> Vabij2Viabj (const char * in){
    FullTensor<double> Vabij (&in[0]);
    std::vector<int> perm (4);
    perm[0] = 2;
    perm[1] = 0;
    perm[2] = 1;
    perm[3] = 3;
    std::vector<int> componentDimensionsNew(4);
    int size = 1;
    for(int mu=0; mu<4; mu++){
       componentDimensionsNew[mu] = Vabij.getComponentDimension(perm[mu]);
       size *= componentDimensionsNew[mu];
    }
    int index;
    int indexNew;
    std::vector<double> vNew(size);
    for(int i = 0; i< Vabij.getComponentDimension(3); i++){
       for(int j = 0; j< Vabij.getComponentDimension(2); j++){
          for(int k = 0; k< Vabij.getComponentDimension(1); k++){
             for(int l = 0; l< Vabij.getComponentDimension(0); l++){

                index = l + Vabij.getComponentDimension(0) * (k + Vabij.getComponentDimension(1)
                           * (j + Vabij.getComponentDimension(2) * i));
                indexNew = j + Vabij.getComponentDimension(2) * (l + Vabij.getComponentDimension(0)
                      * (k + Vabij.getComponentDimension(1) * i));
                vNew[indexNew] = Vabij.getVofIndex(index);
                //std::cout<<"index : "<<index<<"   indexNew : "<<indexNew;
                //std::cout<<"index :"<<index<<" value : "<< Vabij.getVofIndex(index);
             }
          }
       }
    }
    return(FullTensor<double> (componentDimensionsNew, vNew));
 }

} /* end namespace */

