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

#include "CC/Contractions.hpp"

namespace TensorCalculus{

MPSRepresentation<double> contract_21_43 (MPSRepresentation<double> &T,
                                          MPSRepresentation<double> &V){
#ifdef ARGUMENT_CHECKS_ON
   if(T.getComponentDimension(1) != V.getComponentDimension(0)){
      throw std::invalid_argument("componentDimensions for contraction 21 do not fit");
   }
   if(T.getComponentDimension(3) != V.getComponentDimension(2)){
      throw std::invalid_argument("componentDimensions for contraction 43 do not fit");
   }
#endif



   int max1 = std::max(T.getSummation(0), V.getSummation(0));;
   int max2 = std::max(T.getSummation(1), V.getSummation(1));;
   int max3 = std::max(T.getSummation(2), V.getSummation(2));;
   int wSize = std::max(max1,max3);
   wSize *= wSize * max2;
   std::vector<double> w(wSize);

   std::vector< std::vector <double> > newV(4);
   newV[0] = T[0];
   newV[1].resize(V.getComponentDimension(1) * V.getSummation(1) * T.getSummation(1) * T.getSummation(0));
   newV[2].resize(T.getComponentDimension(2) * T.getSummation(1) * V.getSummation(1) * V.getSummation(2));
   newV[3] = V[3];

   std::vector<int> componentDimensions(4);
   componentDimensions[0] = T.getComponentDimension(0);
   componentDimensions[1] = V.getComponentDimension(1);
   componentDimensions[2] = T.getComponentDimension(2);
   componentDimensions[3] = V.getComponentDimension(3);

   std::vector<int> summations(3);
   summations[0] = T.getSummation(0);
   summations[1] = T.getSummation(1) * V.getSummation(1);
   summations[2] = V.getSummation(2);

   // w = w(j'1,j1,j2) = v_1(n,j'1)^T * t_2(n,j1,j2)
   Blas<double>::gemm('t', 'n', V.getSummation(0), T.getSummation(0)*T.getSummation(1), V.getComponentDimension(0),
         1.0, &V(0,0), V.getComponentDimension(0), &T(1,0), T.getComponentDimension(1), 0.0, &w[0], V.getSummation(0));

   //newV[1] = v'2 (i2, j1, j2, j'2) = w (j'1,j1,j2) * v2(i2,j'1,j'2)
   for(int j=0; j<V.getSummation(1); j++){
      Blas<double>::gemm('n', 'n', V.getComponentDimension(1), T.getSummation(0)*T.getSummation(1), V.getSummation(0),
            1.0, &V(1,V.getComponentDimension(1)*V.getSummation(0)*j), V.getComponentDimension(1), &w[0], V.getSummation(0),
            0.0, &newV[1][V.getComponentDimension(1)*T.getSummation(0)*T.getSummation(1)*j], V.getComponentDimension(1));
   }

   // w = w(j3, j'2, j'3)
   Blas<double>::gemm('t', 'n', T.getSummation(2), V.getSummation(1)*V.getSummation(2), V.getComponentDimension(2),
         1.0, &T(3,0), T.getComponentDimension(3), &V(2,0), V.getComponentDimension(2), 0.0, &w[0], T.getSummation(2));

   //newV[2] = t'3(i3, j2, j'2, j'3) = t3(i3,j2,j3) * w(j3,j'2,j'3)
      Blas<double>::gemm('n', 'n', T.getComponentDimension(2)*T.getSummation(1), V.getSummation(1)*V.getSummation(2), T.getSummation(2),
         1.0, &T(2,0), T.getComponentDimension(2)*T.getSummation(1), &w[0], T.getSummation(2),
         0.0, &newV[2][0], T.getComponentDimension(2)*T.getSummation(1));

   return(MPSRepresentation<double>(summations, newV, componentDimensions));
}

double contract_11_24_33_42 (MPSRepresentation<double> &T, MPSRepresentation<double> &V){
   double value = 0;
#ifdef ARGUMENT_CHECKS_ON
   if(T.getComponentDimension(0) != V.getComponentDimension(0)){
      throw std::invalid_argument("componentDimensions for contraction 11 do not fit");
   }
   if(T.getComponentDimension(1) != V.getComponentDimension(3)){
      throw std::invalid_argument("componentDimensions for contraction 24 do not fit");
   }
   if(T.getComponentDimension(2) != V.getComponentDimension(2)){
      throw std::invalid_argument("componentDimensions for contraction 33 do not fit");
   }
   if(T.getComponentDimension(3) != V.getComponentDimension(1)){
      throw std::invalid_argument("componentDimensions for contraction 42 do not fit");
   }
#endif

   throw std::invalid_argument("contract 11_24_33_42 not implemented yet");

   std::vector<int> perm(4); // apply on V
   perm[0] = 0;
   perm[1] = 3;
   perm[2] = 2;
   perm[3] = 1;

   std::vector<int> sumT(4);
   sumT[0] = T.getSummation(0);
   sumT[1] = T.getSummation(0) * T.getSummation(1);
   sumT[2] = T.getSummation(1) * T.getSummation(2);
   sumT[3] = T.getSummation(2);

   std::vector<int> sumV(4);
   sumV[0] = V.getSummation(0);
   sumV[3] = V.getSummation(0) * V.getSummation(1);
   sumV[2] = V.getSummation(1) * V.getSummation(2);
   sumV[1] = V.getSummation(2);

   std::vector< std::vector<double> > w(4);
   w[0].resize(sumT[0]*sumV[0]);
   w[1].resize(sumT[1]*sumV[1]);
   w[2].resize(sumT[2]*sumV[2]);
   w[3].resize(sumT[3]*sumV[3]);

   for(int mu=0; mu<4; mu++){
      Blas<double>::gemm('t', 'n', sumT[mu], sumV[mu], T.getComponentDimension(mu),
            1.0, &T(mu,0), T.getComponentDimension(mu), &V(perm[mu], 0), V.getComponentDimension(perm[mu]),
            0.0, &w[mu][0], sumT[mu]);
   }


   /*
   std::vector<double> t4v2(V.getSummation(0)*V.getSummation(1)*T.getSummation(2));
   std::vector<double> t2v4(T.getSummation(0)*T.getSummation(1)*V.getSummation(2));

   //v2t4 = v2t4(j3, j'1,j'2) = t4(i4,j3)^t * v2(i2, j'1,j'2)
   Blas<double>::gemm('t', 'n', T.getSummation(2), V.getSummation(0)*V.getSummation(1), T.getComponentDimension(3),
         1.0, &T(3,0), T.getComponentDimension(3), &V(1,0), V.getComponentDimension(1), 0.0, &t4v2[0], T.getSummation(2));
   //v4t2 = v4t2(j'3, j1,j2) = v4(i4,j'3)^t * t2(i2, j1,j2)
   Blas<double>::gemm('t', 'n', V.getSummation(2), T.getSummation(0)*T.getSummation(1), V.getComponentDimension(3),
            1.0, &V(3,0), V.getComponentDimension(3), &T(1,0), T.getComponentDimension(1), 0.0, &t2v4[0], V.getSummation(2));
   */
   return(value);
}

double contract_14_22_33_41 (MPSRepresentation<double> &T, MPSRepresentation<double> &V){
#ifdef ARGUMENT_CHECKS_ON
   if(T.getComponentDimension(0) != V.getComponentDimension(3)){
      throw std::invalid_argument("componentDimensions for contraction 14 do not fit");
   }
   if(T.getComponentDimension(1) != V.getComponentDimension(1)){
      throw std::invalid_argument("componentDimensions for contraction 22 do not fit");
   }
   if(T.getComponentDimension(2) != V.getComponentDimension(2)){
      throw std::invalid_argument("componentDimensions for contraction 33 do not fit");
   }
   if(T.getComponentDimension(3) != V.getComponentDimension(0)){
      throw std::invalid_argument("componentDimensions for contraction 41 do not fit");
   }
#endif

   //maybe need to check something here...
   int wSize = V.getSummation(0)*V.getSummation(1)*T.getSummation(0)*T.getSummation(1);
   std::vector<double> w1(wSize);
   std::vector<double> w2(wSize);

   std::vector<double> newV3 (V.getComponentDimension(2)*V.getSummation(1)*T.getSummation(0));
   std::vector<double> newT3 (T.getComponentDimension(2)*T.getSummation(1)*V.getSummation(0));

   //w1(j'3,j1) = v4(n,j'3)^t * t1(n, j1)
   Blas<double>::gemm('t','n', V.getSummation(2), T.getSummation(0), T.getComponentDimension(0),
         1.0, &V(3,0), V.getComponentDimension(3), &T(0,0), T.getComponentDimension(0),
         0.0, &w1[0], V.getSummation(2));
   //w2(j3,j'1) = t4(n,j3)^t * v1(n, j'1)
   Blas<double>::gemm('t', 'n', T.getSummation(2), V.getSummation(0), V.getComponentDimension(0),
         1.0, &T(3,0), T.getComponentDimension(0), &V(0,0), V.getComponentDimension(0),
         0.0, &w2[0], T.getSummation(2));

   //newV3 (n3, j'2, j1) = v3(n3,j'2 |j'3) * w1(j'3, j1)
   Blas<double>::gemm('n', 'n', V.getComponentDimension(2)*V.getSummation(1), T.getSummation(0), V.getSummation(2),
         1.0, &V(2,0), V.getComponentDimension(2)*V.getSummation(1), &w1[0], V.getSummation(2),
         0.0, &newV3[0], V.getComponentDimension(2)*V.getSummation(1));

   //newT3 (n3, j2, j'1) = t3(n3,j2 |j3) * w2(j3, j'1)
   Blas<double>::gemm('n', 'n', T.getComponentDimension(2)*T.getSummation(1), V.getSummation(0), T.getSummation(2),
         1.0, &T(2,0), T.getComponentDimension(2)*T.getSummation(1), &w2[0], T.getSummation(2),
         0.0, &newT3[0], T.getComponentDimension(2)*T.getSummation(1));

   //w1 (j1,j2,j'1,j'2) = ...
   for(int j=0; j<V.getSummation(1); j++){
      Blas<double>::gemm('t', 'n', T.getSummation(0), V.getSummation(0)*T.getSummation(1), V.getComponentDimension(2),
            1.0, &newV3[V.getComponentDimension(2)*j], V.getComponentDimension(2)* V.getSummation(1),
            &newT3[0], T.getComponentDimension(2),
            0.0, &w1[T.getSummation(0)*V.getSummation(0)*T.getSummation(1)*j], T.getSummation(0));
   }

   //w2 (j1,j2,j'1,j'2) = ...
   Blas<double>::gemm('t', 'n', T.getSummation(0)*T.getSummation(1), V.getSummation(0)*V.getSummation(1), V.getComponentDimension(1),
         1.0, &T(1,0), T.getComponentDimension(1), &V(1,0), V.getComponentDimension(1),
         0.0, &w2[0], T.getSummation(0)*T.getSummation(1));

   // < w1 | w2 >
   return innerProduct(w1,w2);
}

}//end namespace
