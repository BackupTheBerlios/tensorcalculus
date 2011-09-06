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

#ifndef MP2_CPP_
#define MP2_CPP_

#include "CC/MP2.hpp"
#include "Tensor/FullTensor.hpp"
#include "CC/C4Interface.hpp"
#include "CC/AOMOtransformation.hpp"
#include "CC/Amplitudes.hpp"
#include <time.h>
#include "Vector/VectorOperators.hpp"
#include "CC/Contractions.hpp"

namespace TensorCalculus{

   double mp2_ft (const char *v_abij, const char *t_abij){
      std::ifstream fin1 (&v_abij[0]);
      if(fin1.fail()){
         throw std::invalid_argument("no v file");
      }
      std::ifstream fin2 (&t_abij[0]);
      if(fin2.fail()){
         throw std::invalid_argument("no t file");
      }

      FullTensor<double> v(&v_abij[0]);
      FullTensor<double> t(&t_abij[0]);

      double E2 = innerProduct(v.getV(), t.getV());
      std::cout.precision(13);
      std::cout<<"E2(AB)  = "<<E2<<std::endl;


      double Etemp = 0.0;
      int indexV;
      int indexT;
      int n1 = v.getComponentDimensions()[0];
      int n2 = v.getComponentDimensions()[1];
      int n3 = v.getComponentDimensions()[2];
      int n4 = v.getComponentDimensions()[3];
      for (int i=0; i < n4; i++){
         for (int j=0; j < n3; j++){
            for (int k=0; k < n2; k++ ){
               for (int l=0; l < n1; l++){
                  indexV = l + n1 * (k + n2 * (j + n3 * i));
                  indexT = l + n1 * (k + n2 * (i + n4 * j));
                  Etemp -= v.getV()[indexV] * t.getV()[indexT];
               }
            }
         }
      }
      Etemp += E2;
      E2 += Etemp;
      Etemp /= 2.0;
      std::cout<<"E2(AA)  = "<<Etemp<<std::endl;
      std::cout<<"E2(tot) = "<<E2<<std::endl;

     return (E2);
   }

   double mp2_mps (MPSRepresentation<double> &v_abij, MPSRepresentation<double> &t_abij, const bool chemical){
      double E2;
      if (!chemical){
         E2 = mp2_mps_phys (v_abij, t_abij);
      } else {
         E2 = mp2_mps_chem (v_abij, t_abij);
      }
      return E2;
   }
   // testing
   double mp2_mps_chem (MPSRepresentation<double> &v_abij, MPSRepresentation<double> &t_abij){
      std:: cout<<"mp2_mps_chem ..."<<std::endl;
      double E2;
      const std::vector<int> &vs = v_abij.getSummations();
      const std::vector<int> &ts = t_abij.getSummations();
      //const std::vector<int> &n  = v_abij.getComponentDimensions();

      std::vector< std::vector <double> > w(4);
      w[0].resize(vs[0]*ts[0]);


      ///////////////////////////////////////////////////////////
      double E2AB = v_abij.tensorScalarProduct(t_abij);

      std::cout<<"E2(AB)  = "<<E2AB<<std::endl;

      //E2 = v_abij.TensorRepresentation<double>::tensorScalarProduct(t_abij);

      //std::cout<<"E2 = "<<E2<<std::endl;

      E2 = E2AB;
      if (v_abij.getComponentDimension(0) == v_abij.getComponentDimension(2)){
         std::cout<<"using (ai|bj)"<<std::endl;
         std::vector<int> perm(4);
         perm[0] = 0;
         perm[1] = 3;
         perm[2] = 2;
         perm[3] = 1;

         std::vector< std::vector <double> > v(4);
         std::vector<std::vector<int> > incidenceMatrix(4);
         //std::vector<int> summations(3);
         std::vector<int> compDim(4);
         for(int mu=0; mu < 4; mu++){
            v[mu] = t_abij.getV()[perm[mu]];
            incidenceMatrix[mu] = t_abij.getIncidenceMatrix()[perm[mu]];
            compDim[mu] = t_abij.getComponentDimension(perm[mu]);
         }
         TensorRepresentation<double> t_perm(t_abij.getSummations(), v, compDim, incidenceMatrix);



         E2 -= v_abij.TensorRepresentation<double>::tensorScalarProduct(t_perm);
      }


      if (v_abij.getComponentDimension(0) == v_abij.getComponentDimension(3)){
         std::cout<<"using (ia|bj)"<<std::endl;
         double testE2 = contract_14_22_33_41(t_abij, v_abij);
         E2 -= testE2;
         std::cout<<"testE2(partial AA)  = "<<testE2<<std::endl;
      }

      std::cout<<"E2(AA)  = "<<E2<<std::endl;
      E2 += E2AB;
      std::cout<<"E2(tot) = "<<E2<<std::endl;
      return E2;
   }

   double mp2_mps_phys (MPSRepresentation<double> &v_abij, MPSRepresentation<double> &t_abij){
      // there is a bug, (maybe lapack reads wrong memory?)
      double E2;

      int vs1 = v_abij.getSummation(0);
      int vs2 = v_abij.getSummation(1);
      int vs3 = v_abij.getSummation(2);
      int ts1 = t_abij.getSummation(0);
      int ts2 = t_abij.getSummation(1);
      int ts3 = t_abij.getSummation(2);

      int n4 = v_abij.getComponentDimension(3);
      int n3 = v_abij.getComponentDimension(2);
      int n2 = v_abij.getComponentDimension(1);
      int n1 = v_abij.getComponentDimension(0);

      if (n1 != n2 || n3 != n4){
         throw std::invalid_argument("wrong format");
      }

      std::vector<double> w12(vs2*ts2);
      std::vector<double> w34(vs2*ts2);
      std::vector<double> w1 ( std::max(vs1, vs3) * std::max(ts1, ts3) * vs2);
      std::vector<double> &w4 = w1;
      std::vector<double> &w3 = w1;

      double tempSize = std::max(n2 * vs2 * ts1, n3 * std::max(vs3 * ts2, ts3 * vs2) );
      std::vector<double> temp(tempSize);

      // w12 //
      //compute w1 = <v1(j1),t1(j'1)>_i=1..n
      Blas<double>::gemm('t', 'n', vs1, ts1, n1, 1.0, &v_abij(0,0), n1, &t_abij(0,0), n1, 0.0, &w1[0], vs1);
      //compute v'2 = sum (j1) w1 * v2;
#pragma omp parallel for
      for(int j=0; j < vs2; j++){
         Blas<double>::gemm('n', 'n', n2, ts1, vs1, 1.0, &v_abij(1,n2*vs1*j), n2, &w1[0], vs1, 0.0, &temp[n2*ts1*j], n2);
      }
      // compute w12 = <temp(j2), t2(j2)>_i,j'=1..n*r'1
      Blas<double>::gemm('t', 'n', vs2, ts2, n1*ts1, 1.0, &temp[0], n1*ts1, &t_abij(1,0), n1*ts1, 0.0, &w12[0], vs2);

      // w34 //
      //compute w4 = <v4(j3), t4(j3)>_i=1..n
      Blas<double>::gemm('t', 'n', vs3, ts3, n4, 1.0, &v_abij(3,0), n4, &t_abij(3,0), n4, 0.0, &w4[0], vs3);
      //compute v'3 = sum (j3) w4 * v3;
      Blas<double>::gemm('n', 'n', n3*vs2, ts3, vs3, 1.0, &v_abij(2,0), n3*vs2, &w4[0], vs3, 0.0, &temp[0], n3*vs2);
      //compute w34 = sum(j'3) <temp(j2), t(j'2)>_i=1..n
#pragma omp parallel for
      for (int j=0; j < ts3; j++){
         Blas<double>::gemm('t', 'n', vs2, ts2, n3, 1.0, &temp[j*n3*vs2], n3, &t_abij(2,j*n3*ts2), n3, 1.0, &w34[0], vs2);
      }

      /*
      E2 = innerProduct (w12, w34);
      std::cout<<"E2(AB) = "<< E2 <<std::endl;
      */

      // w'34 //
      //compute w3 = <v3(j2,j3), t4(j'3)>_i=1..n
      Blas<double>::gemm('t', 'n', ts3, vs2*vs3, n3, 1.0, &t_abij(3,0), n3, &v_abij(2,0), n3, 0.0, &w3[0], ts3);
      //compute v'4 = sum (j3) v4 * t3
      Blas<double>::gemm('n', 't', n4, ts3*vs2, vs3, 1.0, &v_abij(3,0), n4, &w3[0], ts3*vs2, 0.0, &temp[0], n4);
      //compute w34 = 2.0 * w34 - sum(j3) <v'4 , t3>_i=1..n
      int factor = 2.0;
      for (int j=0; j < ts3; j++){
         if (j == 1){
            factor = 1.0;
         }
         Blas<double>::gemm('t', 'n', vs2, ts2, n4, -1.0, &temp[n4*j], n4*ts3, &t_abij(2, n3*ts2*j), n3,
               factor, &w34[0], vs2);
      }

      E2 = innerProduct (w12, w34);
      // if factor == 0 for j == 0 you get E2(AA) * E2(BB)
      std::cout<<"E2(tot) = "<< E2 <<std::endl;

     return E2;
   }

   double mp2_mps (const char *v_abij_mps, const char *t_abij_mps, const bool chemical){
      std::ifstream fin1(&v_abij_mps[0]);
      std::ifstream fin2(&t_abij_mps[0]);

      if(fin1.fail()){
         throw std::invalid_argument("no v_abij file");
      }

      if(fin2.fail()){
         throw std::invalid_argument("no t_abij file");
      }
      MPSRepresentation<double> v_abij(&v_abij_mps[0]);
      MPSRepresentation<double> t_abij(&t_abij_mps[0]);

      return mp2_mps(v_abij, t_abij, chemical);
   }
/*
   double mp2_mps (const double eps, const char *out, const bool chemical){
      if(chemical){
         std::cout<<"not implenemted"<<std::endl;
      }
      else{
         mp2_mps_phys(eps, out);
      }
   }
*/
   double mp2_mps (const double eps, const char *out, const bool chemical){
      std::ofstream fout(&out[0]);
      fout.precision(10);
      time_t start, end;
      int occCount, virtCount, aoCount;
      readOccN(occCount, occCount, virtCount, virtCount);
      aoCount = virtCount + occCount;
      if(chemical){
         fout<<"chemical format used: index monp = (mn|op) = <mo|np> = "
               "int{m(1) o(2) ||1-2||^(-1) n(1) p(2) d12}"<<std::endl;
         fout<<"mps is connected: m-n-o-p"<<std::endl;
      } else {
         fout<<"physical format used: index mnop = <mn|op> = "
                        "int{m(1) n(2) ||1-2||^(-1) o(1) p(2) d12}"<<std::endl<<std::endl;
      }
      fout<<"epsilon : "<<eps<<std::endl;
      fout<<"virt = "<<virtCount<<"   occ = "<<occCount<<std::endl<<std::endl;

      std::vector<double> fullAO;
      time(&start);
      readAO(fullAO, chemical);
      time(&end);
      fout<<"### AO ###"<<std::endl;
      fout<<"readAO :"<<difftime(end, start)<<" s"<<std::endl;

      std::vector<int> compDim(4);
      compDim[0] = aoCount;
      compDim[1] = aoCount;
      compDim[2] = aoCount;
      compDim[3] = aoCount;
      time(&start);
      std::cout<<"MPSdecomposition :"<<std::endl;
      MPSRepresentation<double> mpsAO(fullAO, compDim, eps);
      time(&end);
      using namespace VectorOperators;
      fout<<"MPSdecomposition :"<<difftime(end, start)<<" s"<<std::endl;
      fout<<"componentDimensions :"<<mpsAO.getComponentDimensions()<<std::endl;
      fout<<"summations :"<<mpsAO.getSummations()<<std::endl<<std::endl;

      //todo: free fullAO

      std::vector<char> occState(4);
      if(chemical){
         //test
         //occState[0] = 'v';
         //occState[1] = 'o';
         occState[0] = 'o';
         occState[1] = 'v';
         occState[2] = 'v';
         occState[3] = 'o';
      } else {
         occState[0] = 'v';
         occState[1] = 'v';
         occState[2] = 'o';
         occState[3] = 'o';
      }
      fout<<"### occState ###"<<std::endl;
      fout<<occState<<std::endl<<std::endl;

      std::vector<int> compDimVirt(4);
      std::vector<int> compDimOcc(4);
      std::vector<int> compDimV(4);

      for(int mu=0; mu<4; mu++){
         compDimVirt[mu] = virtCount;
         compDimOcc[mu]  = occCount;
         if(occState[mu] == 'v'){
            compDimV[mu] = virtCount;
         } else {
            compDimV[mu] = occCount;
         }
      }

      MPSRepresentation<double> mpsOcc(mpsAO.getSummations(), compDimOcc);
      MPSRepresentation<double> mpsVirt(mpsAO.getSummations(), compDimVirt);
      MPSRepresentation<double> mpsV(mpsAO.getSummations(), compDimV);

      time(&start);
      AO2MO(mpsAO, mpsVirt, mpsOcc);
      makeV(occState, mpsV, mpsVirt, mpsOcc);
      time(&end);
      fout<<"### Transformation ###"<<std::endl;
      fout<<"AO2MO + makeV :"<<difftime(end, start)<<" s"<<std::endl;

      std::cout<<"### V_ab_ij ###"<<std::endl;
      time(&start);
      mpsV.truncateMPS(eps, 'e', true);
      time(&end);
      fout<<"### V_ab_ij ###"<<std::endl;
      fout<<"truncation : "<<difftime(end, start)<<" s"<<std::endl;
      fout<<"componentDimensions :"<< mpsV.getComponentDimensions()<<std::endl;
      fout<<"summations :"<< mpsV.getSummations()<<std::endl<<std::endl;

      std::cout<<"### V_ij_kl ###"<<std::endl;
      time(&start);
      mpsOcc.truncateMPS(eps, 'e', true);
      time(&end);
      fout<<"### V_ij_kl ###"<<std::endl;
      fout<<"truncation : "<<difftime(end, start)<<" s"<<std::endl;
      fout<<"componentDimensions :"<< mpsOcc.getComponentDimensions()<<std::endl;
      fout<<"summations :"<< mpsOcc.getSummations()<<std::endl<<std::endl;

      std::cout<<"### V_ab_cd ###"<<std::endl;
      time(&start);
      mpsVirt.truncateMPS(eps, 'e', true);
      time(&end);
      fout<<"### V_ab_cd ###"<<std::endl;
      fout<<"truncation : "<<difftime(end, start)<<" s"<<std::endl;
      fout<<"componentDimensions :"<< mpsVirt.getComponentDimensions()<<std::endl;
      fout<<"summations :"<< mpsVirt.getSummations()<<std::endl<<std::endl;

      MPSRepresentation<double> mpsT (mpsV);
      std::cout<<"### T_ab_ij ###"<<std::endl;
      time(&start);
      V2T(mpsT, occState);
      time(&end);
      fout<<"### T_ab_ij ###"<<std::endl;
      fout<<"V2T : "<<difftime(end, start)<<" s"<<std::endl;
      fout<<"componentDimensions :"<< mpsT.getComponentDimensions()<<std::endl;
      fout<<"summations :"<< mpsT.getSummations()<<std::endl;

      time(&start);
      mpsT.truncateMPS(eps, 'e', true);
      time(&end);
      fout<<"truncation : "<<difftime(end, start)<<" s"<<std::endl;
      fout<<"componentDimensions :"<< mpsT.getComponentDimensions()<<std::endl;
      fout<<"summations :"<< mpsT.getSummations()<<std::endl<<std::endl;

      time(&start);
      double E2 = mp2_mps(mpsV, mpsT, chemical);
      time(&end);
      fout<<"### MP2 Energy ###"<<std::endl;
      fout<<"contraction : "<<difftime(end, start)<<" s"<<std::endl;
      fout<<"E2 = "<<E2<<std::endl<<std::endl;

      //test one term of amplitude update
      /*
      if(chemical){
         time(&start);
         MPSRepresentation<double> T_new (contract_21_43(mpsT, mpsOcc));
         time(&end);
         fout<<"### T_ab_ij partial update(t_ab^ef * v_ef^ij ###"<<std::endl;
         fout<<"contraction 21_43: "<<difftime(end, start)<< " s"<<std::endl;
         fout<<"componentDimensions : "<<T_new.getComponentDimensions()<<std::endl;
         fout<<"summations : "<<T_new.getSummations()<<std::endl;
         time(&start);
         T_new.truncateMPS(eps);
         time(&end);
         fout<<"truncation : "<<difftime(end, start)<<" s"<<std::endl;
         fout<<"summations : "<<T_new.getSummations()<<std::endl;
      }
      */
      return E2;
   }

   double mp2_mps_canonical (const double eps, const char *out, const bool chemical){
      std::ofstream fout (&out[0]);
      fout.precision(12);
      std::vector<char> occState (4);
      if(chemical){
         occState[0] = 'o';
         occState[1] = 'v';
         occState[2] = 'v';
         occState[3] = 'o';
      } else {
         occState[0] = 'v';
         occState[1] = 'v';
         occState[2] = 'o';
         occState[3] = 'o';
      }
      if(chemical){fout<<"### chemical ###"<<std::endl;}
      else {fout<<"### physical ###"<<std::endl;}

      using namespace VectorOperators;
      fout << "occState : " << occState << std::endl;

      fout<<"eps : "<<eps<<std::endl<<std::endl;
      FullTensor<double> ftV;
      FullTensor<double> ftT;
      if(chemical){
         ftV = Vabij2Viabj ();
         ftT = Vabij2Viabj ("t_abij.ten");
      } else {
         ftV = FullTensor<double>("v_abij.ften");
         ftT = FullTensor<double>("t_abij.ten");
      }
      MPSRepresentation<double> mpsV (ftV, eps);
      fout<<"V : componentDimensions : "<<mpsV.getComponentDimensions()<<std::endl;
      MPSRepresentation<double> mpsT (ftT, eps);
      fout<<"T : componentDimensions : "<<mpsT.getComponentDimensions()<<std::endl;
      /*
      MPSRepresentation<double> mpsT (mpsV);
      V2T (mpsT, occState);
      mpsT.truncateMPS(eps, 'e', true);
      */
      /*
      //to check
      if(!chemical){
         std::vector<double> fullT (mpsT.evaluate());
         FullTensor<double> ftT ("t_abij.ten");
         using namespace VectorOperators;
         double normT = l2_norm(ftT.getV());
         double dist = l2_norm(ftT.getV()-fullT);
         std::cout<<"dist : "<< dist<<"    rel dist : "<<dist / normT<<std::endl;
         double E = mp2_ft("v_abij.ften", "t_abij.ten");
         std::cout<<"E2 using v_abij.ften and t_abij.ten : "<<E<<std::endl;
      }
      */
      double E2 = mp2_mps(mpsV, mpsT, chemical);

      std::cout<<"E(tot) = "<<E2<<std::endl;
      fout<<"E2 (tot) :"<<E2<<std::endl;
      return E2;
   }

}//end namespace

#endif
