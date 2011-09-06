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

#include "CC/SpinAdapt.hpp"
#include "CC/C4Interface.hpp"

namespace TensorCalculus{

void readAOSpinAdapted (std::vector<double> &fullAdaptAO){
   std::vector<double> fullAO;
   readAO(fullAO);
   int occA, occB, virtA, virtB;
   readOccN(occA, occB, virtA, virtB);
   int aoCount = occA + virtA;
   spinAdaptedAO(aoCount, fullAO, fullAdaptAO);
} //end readAOSpinAdapted

void spinAdaptedAO (const int aoCount, const std::vector<double> &fullAO, std::vector<double> &fullAdaptAO){
   unsigned int size = aoCount * aoCount * aoCount * aoCount;
   if(fullAdaptAO.size() < size){
      fullAdaptAO.resize(size);
   }
   Blas<double>::copy(size, &fullAO[0], 1, &fullAdaptAO[0], 1);
   Blas<double>::scal(size, 2.0, &fullAdaptAO[0], 1);
   int index1;
   int index2;
   for (int i=0; i<aoCount; i++){
      for (int j=0; j<aoCount; j++){
         for (int k=0; k<aoCount; k++){
            for (int l=0; l<aoCount; l++){
               index1 = l+aoCount*(k + aoCount * (j + aoCount*i));
               index2 = l+aoCount*(k + aoCount * (i + aoCount*j));
               fullAdaptAO[index1] -= fullAO[index2];
            }
         }
      }
   }
} //end spinAdaptedAO


void readAOSpinAdaptedchemical (std::vector<double> &fullAdaptAO){
   std::vector<double> fullAO;
   readAOchemical(fullAO);
   int occA, occB, virtA, virtB;
   readOccN(occA, occB, virtA, virtB);
   int aoCount = occA + virtA;
   spinAdaptedAOchemical(aoCount, fullAO, fullAdaptAO);
} // end readAOSpinAdaptedchemical

void spinAdaptedAOchemical (const int aoCount, const std::vector<double> &fullAO, std::vector<double> &fullAdaptAO){
   unsigned int size = aoCount * aoCount * aoCount * aoCount;
   if(fullAdaptAO.size() < size){
      fullAdaptAO.resize(size);
   }
   Blas<double>::copy(size, &fullAO[0], 1, &fullAdaptAO[0], 1);
   //Blas<double>::scal(size, 2.0, &fullAdaptAO[0], 1);
   int index1;
   int index2;
   for (int i=0; i<aoCount; i++){
      for (int j=0; i<aoCount; i++){
         for (int k=0; i<aoCount; i++){
            for (int l=0; i<aoCount; i++){
               index1 = l+aoCount*(k + aoCount * (j + aoCount*i));
               index2 = l+aoCount*(i + aoCount * (j + aoCount*k));
               fullAdaptAO[index1] -= fullAO[index2];
            }
         }
      }
   }
} // end spinAdaptedAOchemical


} // end namespace TensorCalculus
