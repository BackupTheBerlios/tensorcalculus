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

#include "CC/AOMOtransformation.hpp"

#include "Representation/MPSRepresentation.hpp"
#include "Representation/CPTensorRepresentation.hpp"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>


#include <iostream>
#include "Tensor/FullTensor.hpp"
#include "Matrix/MatrixOperators.hpp"
#include "Vector/VectorOperators.hpp"

#include "BlasInterface.hpp"
#include "CC/C4Interface.hpp"
#include "Representation/MPSDecomposition.hpp"
#include "CC/Amplitudes.hpp"
#include "CC/SpinAdapt.hpp"

#include "Utilities/Random.hpp"

using namespace TensorCalculus;
using namespace VectorOperators;

// Bsp Mike
void mikeD3 (int k, int n);
void TestingWithMikeD3 ();
void random (std::vector<double> &randomNumbers);

//Bsp Stefan
void makeMPSforStefan();

void c4AOchemicalDecomp ();

void testFullMP2 ();
void testFullMP2_b();

// tests

void mpsDecompTest ();
void fullTensorTest();
void mpsTest();
void c4Test();
void c4ReadTest();
void c4GetFullMO();
void c4inverse();
void mpsTruncTest();
void vidalTest();
void qrfTest();
void matrixTest();
void svd_i1i2_i3i4();
void testEdgeAdd();


void MPSreadWrite ();
void MPSwrite ();

void testHadamard ();
void testAdd();

void testVabij2viabj();

//richtige berechnungen
void mpsAOtruncation ();
void mpsVabijTruncation ();



int main() {
   std::cout << "MPSTest..." << std::endl << std::endl;

   //mikeD3 (k,n);
   //TestingWithMikeD3 ();
   //c4AOchemicalDecomp ();
   //MPSwrite ();
   //testFullMP2 ();
   //testFullMP2_b ();

   //mpsAOtruncation ();
   //mpsVabijTruncation ();

   //mpsDecompTest ();
   //svd_i1i2_i3i4();

   //fullTensorTest();
   //mpsTest();
   //c4Test();
   //c4GetFullMO();
   //c4inverse();
   //c4ReadTest();
   //mpsTruncTest();
   //vidalTest();
   //qrfTest();
   //matrixTest();

   //testHadamard ();
   //testAdd();

   //testEdgeAdd();

   //mpsTruncTest();

   //makeMPSforStefan();

   testVabij2viabj();

   std::cout << "MPSTest...finish" << std::endl;
   return 0;
}
/* constructor entfernt, truncate mps nicht fertig
void mpsTruncTest(){
      TensorCalculus::FullTensor<double> ft("h2o.ften");
      std::vector<double> in = ft.getV();
      std::vector<int> compDim =ft.getComponentDimensions();
      using namespace VectorOperators;

      TensorCalculus::MPSRepresentation<double> mps(in, compDim, 0.0);
      TensorCalculus::MPSRepresentation<double> mps2(mps, 0.0);
      std::cout << "summations : " << mps.getSummations() << std::endl;
      //std::cout << "||A|| = " << l2_norm(copy) << std::endl;
      //std::cout << "||A|| = " << l2_norm(mps.evaluate()) << std::endl;
      //std::cout << "||A'|| = " << l2_norm(mps2.evaluate()) << std::endl;
      std::cout << "||A-A'|| = " << l2norm(mps,mps2) << std::endl;
      //std::cout << "evaluate ||A-A'|| = " << l2_norm(mps.evaluate()-mps2.evaluate()) << std::endl;

}
*/
void vidalTest(){
   std::cout<<"#### Test vidal ####"<<std::endl;
   TensorCalculus::FullTensor<double> ft("h2o.ften");
      std::vector<double> in = ft.getV();
      std::vector<double> copy = in;
      std::vector<int> compDim =ft.getComponentDimensions();
      using namespace VectorOperators;

      TensorCalculus::MPSRepresentation<double> mps(in, compDim, 0.0);
      std::cout << "summations : " << mps.getSummations() << std::endl;
      double normA = l2_norm(copy);
      std::cout << "||A|| = " << normA << std::endl;
      //std::cout << "||A|| = " << l2_norm(mps.evaluate()) << std::endl;
      double normDif = l2_norm(mps.evaluate()-copy);
      std::cout << "evaluate ||A-A'|| = " << normDif<< std::endl;
      std::cout << "evaluate ||A-A'||/||A|| = " << normDif/normA << std::endl;
      std::cout<<"#### end Test vidal ####"<<std::endl;
}



void qrfTest() {
   // funktionierende QR Faktorisierung und wieder zusammensetzung
   using namespace VectorOperators;
   int m = 4;
   int n = 3;
   int min = std::min(m,n);

   std::vector<double> T(m*n);
   for(int i=0; i<m*n; i++){
    T[i]=i;
   }

   std::cout<<"T : "<<T<<std::endl;

   std::vector<double> tau(min);
   std::vector<double> work(1);
   int lwork;

   Lapack<double>::geqrf(m,n,&T[0],m,&tau[0],&work[0], -1);
   lwork = work[0];
   work.resize(lwork);
   std::cout<<"QR Faktorisation"<<std::endl;
   Lapack<double>::geqrf(m,n,&T[0],m,&tau[0],&work[0], lwork);
   std::cout<<"QR Faktorisation FINISH"<<std::endl;
   std::vector<double> R(n*n);
   std::vector<double> Q(m*min);

   for(int i=0; i<n; i++){
    for(int j=0; j<=i && j<n; j++){
     R[j+i*n] = T[j+i*m];
    }
   }
   std::cout<<"T': "<<T<<std::endl;
   std::cout<<"R : "<<R<<std::endl;
   std::cout<<"Tau: "<<tau<<std::endl;

   Q = giveQ(m,n,&T[0],m,&tau[0]);
   std::cout<<"Q : "<<Q<<std::endl;
   Blas<double>::gemm('n','n',m, n, n, 1.0, &Q[0], m, &R[0], n, 0.0, &T[0], m);
   std::cout<<"T : "<<T<<std::endl;

   std::cout<<std::endl;

}

void mpsTest() {
   long n1 = 3;
   long n2 = 3;
   long n3 = 3;
   long n4 = 3;
   long dim = n1 * n2 * n3 * n4;

   std::vector<double> gr(dim);
   std::vector<int> compDim(4);
   compDim[0] = n1;
   compDim[1] = n2;
   compDim[2] = n3;
   compDim[3] = n4;

   for (int i = 0; i < dim; i++) {
      gr[i] = i / 7;
   }

   //TensorCalculus::FullTensor<double> ft("h2o.ften");
   //std::vector<double> in = ft.getV();
   std::vector<double> copy = gr;
   //std::vector<int> compDim =ft.getComponentDimensions();
   using namespace VectorOperators;
   std::cout << "compDim : " << compDim << std::endl;
   TensorCalculus::MPSRepresentation<double> mps(gr, compDim, 0.0);
   std::cout << "summations : " << mps.getSummations() << std::endl;
   std::cout << "||A|| = " << l2_norm(copy) << std::endl;
   std::cout << "||A-A'|| = " << l2_norm(copy - mps.evaluate()) << std::endl;

   /*
    TensorCalculus::FullTensor<double> ft("c3h8.ften");
    //TensorCalculus::FullTensor<double> ft("c3h8_a.ften");
    //TensorCalculus::FullTensor<double> ft("h2o.ften");
    using namespace VectorOperators;
    std::vector<double> v = ft.getV();

    TensorCalculus::MPSRepresentation<double> mpst(v, ft.getComponentDimensions(),0.0);
    v = ft.getV();
    //TensorCalculus::MPSRepresentation<double> mpst2(mpst.getSummations(), mpst.getComponentDimensions());
    */

}
void c4ReadTest() {
   int aoCount = 15;
   std::vector<double> ao(aoCount*aoCount*aoCount*aoCount);

   readAO(aoCount, ao); //

   std::vector<double> occAOMO;
   std::vector<double> virtAOMO;

   int occ = 10;
   int virt = 5;

   //readTransformationMatrix(virt, occ, virtAOMO, occAOMO);

   using namespace VectorOperators;
   //printMatrix(occ, aoCount, &occAOMO[0], occ);
   //printMatrix(virt, aoCount, &virtAOMO[0], virt);
   std::vector<int> compDim(4);
   for(int i=0; i<4; i++){
      compDim[i] = aoCount;
   }
   /*
   std::vector<double> ao2=ao;
   std::vector<double> copy = ao;
   double normA = l2_norm(copy);
   */
   TensorCalculus::FullTensor<double> ft (compDim, ao);
   ft.write2disk("fullAO.ften");
   /*
   TensorCalculus::MPSRepresentation<double> mps(ao, compDim, 0.0);
   std::cout<<"componentDimensions : "<<compDim<<std::endl;
   std::cout<<"summations          : "<<mps.getSummations()<<std::endl;

   TensorCalculus::MPSRepresentation<double> mps2(ao2, compDim, 1e-2);
   std::cout<<"componentDimensions : "<<compDim<<std::endl;
   std::cout<<"summations          : "<<mps2.getSummations()<<std::endl;

   std::cout<<"||A||  = "<<normA<<std::endl;
   double normDif = l2_norm(copy-mps2.evaluate());
   std::cout<<"||A-Atrunc||/||A||  = "<<normDif/normA<<std::endl;
   */

   /*
   std::vector<int> index(4);
   for(int i=0; i< aoCount*aoCount*aoCount*aoCount; i++){
      if(std::fabs(ao[i])<1e-10){
         index[0] = i%aoCount;
         index[1] = ((i-index[0])/aoCount)%aoCount;
         index[2] = ((i-index[0]-index[1]*aoCount)/aoCount/aoCount)%aoCount;
         index[3] = (i-index[0]-index[1]*aoCount-index[2]*aoCount)/aoCount/aoCount/aoCount;
         using namespace VectorOperators;
         std::cout<<index<<"   "<<ao[i]<<" || ";
      }
   }
   */

}

void c4Test() {
   int aoCount = 15;
   int occ = 5; // Ammoniak: 10 e⁻, 5 besetze Orbitale
   int virt = 10;

   std::cout<<"AOMOtransformation"<<std::endl;

   TensorCalculus::AOMOtransformation transformation(virt, occ);
   transformation.AO2MO();

   std::cout<<"AOMOtransformation int ... finish"<<std::endl;

   TensorCalculus::MPSRepresentation<double> Vabij(transformation.getVabij());

   std::cout<<"AOMOtransformation getVabij ... finish"<<std::endl;

   TensorCalculus::MPSRepresentation<double> Vabij2(transformation.AO2Vabij());

   std::cout<<"AOMOtransformation AO2Vabij ... finish"<<std::endl;
   std::cout<<"AOMOtransformation ... finish"<<std::endl;

   std::cout<<"read and evaluate v_abij"<<std::endl;

   TensorCalculus::FullTensor<double> Wabij("v_abij.ften");
   //TensorCalculus::FullTensor<double> Wiabj("v_iabj.ften");
   std::vector<double> w(Wabij.getV());

   std::vector<double> v1(Vabij.evaluate());
   std::vector<double> v2(Vabij2.evaluate());
   /*
   std::vector<double> full_v_abij = Wabij.getV();
   TensorCalculus::MPSRepresentation<double> readV(full_v_abij, Vabij.getComponentDimensions(), 0.0);
   */
   std::cout<<"read and evaluate v_abij ... finish"<<std::endl;

   std::cout<<"compute norm"<<std::endl;
   using namespace VectorOperators;
   double norm = l2_norm(v1-v2);
    std::cout<<"volle transformation - einfache transformation"<<std::endl;
    std::cout<<"norm  : "<<norm<<std::endl;

   norm = l2_norm(v1-w);
   std::cout<<"volle transformation - eingelesene, zerlegte werte"<<std::endl;
   std::cout<<"norm  : "<<norm<<std::endl;

   norm = l2_norm(v2-w);
   std::cout<<"einfache transformation - eingelesene, zerlegte werte"<<std::endl;
   std::cout<<"norm  : "<<norm<<std::endl;

/*
   std::vector<int> index(4);
   std::vector<int> index2(4);
  for(int i=0; i< occ*occ*virt*virt; i++){
     index[0] = i%virt;
     index[1] = ((i-index[0])/virt)%virt;
     index[2] = ((i-index[0]-index[1]*virt)/virt/virt)%occ;
     index[3] = (i-index[0]-index[1]*virt-index[2]*virt*virt)/virt/virt/occ;
     using namespace VectorOperators;
     index2[0] = 1+occ+index[0];
     index2[1] = 1+occ+index[1];
     index2[2] = 1+index[2];
     index2[3] = 1+index[3];
     std::cout<<index<<" || "<<index2<<" :  "<<v2[i]<<" || "<<std::endl;
  }
  */

}



void matrixTest(){

   int m = 5;
   int n = 3;

   std::vector<double> A;

   for (int j=0; j<n*m; j++){
      A.push_back(j);
   }
   printMatrix(m,n,&A[0],m);
}

void c4GetFullMO() {

   bool ausgabeMO = false;

   int occ;
   int virt;
   readOccN(occ, occ, virt, virt);
   int aoCount = virt + occ;

   TensorCalculus::AOMOtransformation transformation(occ, virt);
   transformation.AO2MO();
   TensorCalculus::MPSRepresentation<double> MO(transformation.getFullMO());
   MO.write2disk("fullMO.MPSten");
   std::vector<double> AO;
   readAO (AO);
   std::vector<int> compdim(4);
   compdim[0] = aoCount;
   compdim[1] = aoCount;
   compdim[2] = aoCount;
   compdim[3] = aoCount;
   MPSRepresentation<double> mpsAO (AO, compdim);
   compdim[0] = virt;
   compdim[1] = virt;
   compdim[2] = virt;
   compdim[3] = virt;

   MPSRepresentation<double> virtMO(mpsAO.getSummations(), compdim);

   compdim[0] = occ;
   compdim[1] = occ;
   compdim[2] = occ;
   compdim[3] = occ;

   MPSRepresentation<double> occMO(mpsAO.getSummations(), compdim);

   AO2MO(mpsAO, virtMO, occMO);
   virtMO.write2disk("virtMO.MPSten");
   occMO.write2disk("occMO.MPSten");
   //transformation.fullTensorAO2fullTensorMO();

   //std::vector<double> fullMO(MO.evaluate());
   //std::vector<double> fullMO2(transformation.getFullTensorMO());
   //printMatrix(aoCount, aoCount, &transformation.getFullTransposeAOMO()[0], aoCount);

   using namespace VectorOperators;
   //std::cout<<"norm (decomposed fullMO - fullTensorMO) = "<<l2_norm(fullMO-fullMO2)<<std::endl;

   if(ausgabeMO){
      std::cout<<"Fulltensor in d = 4"<<std::endl;
      std::cout<<"n[0] = "<<aoCount <<std::endl;
      std::cout<<"n[1] = "<<aoCount <<std::endl;
      std::cout<<"n[2] = "<<aoCount <<std::endl;
      std::cout<<"n[3] = "<<aoCount <<std::endl;

      std::vector<int> index(4);
      for(int i=0; i< aoCount*aoCount*aoCount*aoCount; i++){
         index[0] = i%aoCount;
         index[1] = ((i-index[0])/aoCount)%aoCount;
         index[2] = ((i-index[0]-index[1]*aoCount)/aoCount/aoCount)%aoCount;
         index[3] = (i-index[0]-index[1]*aoCount-index[2]*aoCount)/aoCount/aoCount/aoCount;

         index[0] += 1;
         index[1] += 1;
         index[2] += 1;
         index[3] += 1;

         using namespace VectorOperators;
         //std::cout<<index<<"   "<<fullMO2[i]<<std::endl;
      }

   }
}
/*
void c4inverse(){

   int occ = 5;
   int virt = 10;
   int aoCount = virt + occ;

   TensorCalculus::AOMOtransformation transformation(occ, virt);
   transformation.AO2MO();

   std::vector<double> occAOMO(transformation.getOccAOMO());
   std::vector<double> virtAOMO(transformation.getVirtAOMO());
   std::vector<double> fullTransposeAOMO(transformation.getFullTransposeAOMO());

   std::vector<double> copy(fullTransposeAOMO);
   std::vector<double> erg(aoCount*(aoCount));

   std::vector<double> invOccAOMO(occ*aoCount);

   int info;
   info = invertA(aoCount, aoCount, fullTransposeAOMO, aoCount);
   std::cout<<"invert fullTransposeAOMO : "<<info<<std::endl;
   Blas<double>::gemm('n','n',aoCount, aoCount, aoCount, 1.0, &fullTransposeAOMO[0], aoCount, &copy[0], aoCount, 0.0, &erg[0], aoCount);
   printMatrix(aoCount, aoCount, &erg[0], aoCount, 1e-12);

   copy = occAOMO;

   transpose(occ, aoCount, &copy[0], &occAOMO[0]);
   copy = occAOMO;

   info = invertA(aoCount, occ, occAOMO, aoCount);
   std::cout<<"invert occTransposeAOMO : "<<info<<std::endl;

   Blas<double>::gemm('n','n',aoCount, aoCount, occ, 1.0, &copy[0], aoCount, &occAOMO[0], occ, 0.0, &erg[0], aoCount);
   printMatrix(aoCount, aoCount, &erg[0], aoCount, 1e-12);
}
*/

void mikeD3 (int k, int n){

   std::vector<double> zufall (2*k*n + n);
   random(zufall);

   std::vector<double> v1(k*n);
   std::vector<double> v2(n);
   std::vector<double> v3(k*n);
   for(int i=0; i<n; i++){
      v2[i] = zufall[i];
   }
   for(int j=0; j<k*n; j++){
      v1[j] = zufall[n+j];
      v3[j] = zufall[n+n*k+j];
   }
   /*
   std::cout<<"v1 :"<<std::endl;
   printMatrix(n,k,&v1[0],n);

   std::cout<<"v2 :"<<std::endl<<v2<<std::endl;
   std::cout<<"v3 :"<<std::endl;
   printMatrix(n,k,&v3[0],n);
   */
   std::vector<double> fullV(n*n*n);
   std::vector<double> fullVPerm(n*n*n);
   double sum;
   for(int i1=0; i1<n; i1++){
      for(int i2=0; i2<n; i2++){
         for(int i3=0; i3<n; i3++){
            sum = 0.0;
            for(int j=0; j<k; j++){
               sum += v1[i1+n*j]*v2[i2]*v3[i3+n*j];
            }
            fullV[i1+n*i2+n*n*i3] = sum;
            fullVPerm[i1+n*i3+n*n*i2] = sum;
         }
      }
   }

   std::vector<int> compDim(3);
   compDim[0] = n;
   compDim[1] = n;
   compDim[2] = n;

   double normV = l2_norm(fullV);
   fullV /= normV;

   TensorCalculus::MPSRepresentation<double> mpsV (fullV, compDim, 1e-12);
   TensorCalculus::MPSRepresentation<double> mpsVPerm (fullVPerm, compDim, 1e-12);

   //std::cout<< "||fullV - mpsV||" << l2_norm(fullV - mpsV.evaluate())<<std::endl;
   std::cout<<"summations     : "<<mpsV.getSummations()<<std::endl;
   std::cout<<"summationsPerm : "<<mpsVPerm.getSummations()<<std::endl;

}

void random (std::vector<double> &randomNumbers){
   int l = randomNumbers.size();
   time_t t;
   time(&t);
   srand((unsigned int)t);

   for (int i=0; i<l; i++){
      randomNumbers[i] = (double)rand();
      randomNumbers[i] /= (double)rand();
   }
}

void TestingWithMikeD3 (){
   int n1 = 10;
   int n2 = 100;
   //int n3 = 300;

   for(int k = 1 ; k < 2*n1; k++){
      std::cout<<"k = "<<k<<" n = "<<n1<<"  ";
      mikeD3(k, n1);
   }

   for(int k = 10 ; k < 2*n2; k+=20){
      std::cout<<"k = "<<k<<" n = "<<n2<<"  ";
      mikeD3(k, n2);
   }
   /*
   for(int k = 50 ; k < 2*n3; k+=100){
      std::cout<<"k = "<<k<<" n = "<<n3<<"  ";
      mikeD3(k, n3);
   }
   */
}


void fullTensorTest(){
   int n = 4;
   int d = 3;
   std::vector<int> v(n*n*n);
   std::vector<int> compDim(d);
   compDim[0] = n;
   compDim[1] = n;
   compDim[2] = n;

   TensorCalculus::FullTensor<int> ft(compDim, v);

   std::cout<<"v : "<<ft.getV()<<std::endl;

   std::vector<int> index (d);
   index [0] = 2;
   index [1] = 3;
   index [2] = 0;
   ft.setVofIndex(index, 1);

   std::cout<<"v : "<<ft.getV()<<std::endl;

   ft.write2disk("ft.ften");

}

void mpsAOtruncation () {
   double eps = 1e-4;

   int occ;
   int virt;
   int occB;
   int virtB;
   readOccN(occ, occB, virt, virtB);

   int aoCount  = occ + virt;
   std::vector<double> fullao;

   readAO(aoCount, fullao);

   std::vector<int> compDim(4);
   compDim[0] = aoCount;
   compDim[1] = aoCount;
   compDim[2] = aoCount;
   compDim[3] = aoCount;

   TensorCalculus::MPSRepresentation<double> mpsAO(fullao, compDim, eps);
   std::cout<<"truncate AO"<<std::endl;
   std::cout<<"eps : "<<eps<<"  summations : "<<mpsAO.getSummations()<<std::endl;
}

void mpsVabijTruncation (){
   double eps = 1e-4;
   TensorCalculus::FullTensor<double> ftAO("v_abij.ften");
   TensorCalculus::MPSRepresentation<double> mpsAO(ftAO.getV(), ftAO.getComponentDimensions(), eps);
   std::cout<<"truncate v_abij"<<std::endl;
   std::cout<<"eps : "<<eps<<"  summations : "<<mpsAO.getSummations()<<std::endl;
}


void mpsDecompTest (){
   std::ifstream fin ("v_abij.ften");
   if(fin.fail()){
      throw std::invalid_argument("no input file");
   }
   TensorCalculus::FullTensor<double> ftAO("v_abij.ften");
   TensorCalculus::MPSDecomposition<double> decompAO (ftAO.getV(), ftAO.getComponentDimensions());
   //std::cout<<"computeVidalDecomposition"<<std::endl;
   decompAO.computeVidalDecomposition();

   //std::cout<<"setMPS"<<std::endl;
   TensorCalculus::MPSRepresentation<double> mpsAO(decompAO.getSummations(), decompAO.getV(), decompAO.getComponentDimensions());


   TensorCalculus::MPSDecomposition<double> decompAO2 (ftAO.getV(), ftAO.getComponentDimensions());
   std::vector <int> sum(3);
   sum[0] = 10;
   sum[1] = 15;
   sum[2] = 5;
   decompAO2.computeVidalDecomposition(sum);

   TensorCalculus::MPSRepresentation<double> mpsAO2(decompAO2.getSummations(), decompAO2.getV(), decompAO2.getComponentDimensions());


   TensorCalculus::MPSDecomposition<double> decompAO4 (ftAO.getV(), ftAO.getComponentDimensions());

   sum[0] = 5;
   sum[1] = 10;
   sum[2] = 5;
   decompAO4.computeVidalDecomposition(sum);

   TensorCalculus::MPSRepresentation<double> mpsAO4(decompAO4.getSummations(), decompAO4.getV(), decompAO4.getComponentDimensions());



   std::vector<double> evalMPS(mpsAO.evaluate());
   std::cout<<"summations : "<<mpsAO.getSummations()<<std::endl;
   std::cout<<"norm = "<<l2_norm(evalMPS-ftAO.getV())<<std::endl;

   std::vector<double> evalMPS2(mpsAO2.evaluate());
   std::cout<<"summations : "<<mpsAO2.getSummations()<<std::endl;
   std::cout<<"norm = "<<l2_norm(evalMPS2-ftAO.getV())<<std::endl;


   std::cout<<"norm 12 = "<<l2_norm(evalMPS2-evalMPS)<<std::endl;

   mpsAO2.performALS(mpsAO);
   evalMPS2 = mpsAO2.evaluate();
   std::cout<<"norm 12 = "<<l2_norm(evalMPS2-evalMPS)<<std::endl;


   // dmrg
   evalMPS2 = mpsAO2.evaluate();
   std::cout<<"norm 12 = "<<l2_norm(evalMPS2-evalMPS)<<std::endl;

   std::vector<double> evalMPS3(mpsAO2.evaluate());
   std::cout<<"summations : "<<mpsAO2.getSummations()<<std::endl;
   std::cout<<"norm = "<<l2_norm(evalMPS3-ftAO.getV())<<std::endl;

}

void c4AOchemicalDecomp (){

   int occa, occb, virta, virtb;
   readOccN(occa, occb, virta, virtb);
   int aoCount = virta + occa;
   std::vector<double> fullAO(aoCount*aoCount*aoCount*aoCount);

   readAOchemical (aoCount, fullAO);

   std::vector<int> compDim(4);
   compDim[0] = aoCount;
   compDim[1] = aoCount;
   compDim[2] = aoCount;
   compDim[3] = aoCount;
   double eps = 1e-6;
   TensorCalculus::MPSRepresentation<double> mps(fullAO, compDim, eps);
   std::cout<<" eps = "<< eps << std::endl;
   std::cout<<"summations : "<< mps.getSummations()<<std::endl;

   double norm = l2_norm(fullAO-mps.evaluate())/l2_norm(fullAO);
   std::cout<<"norm  || A - MPS(A) || / ||A|| : "<< norm<<std::endl;
}

void svd_i1i2_i3i4(){
   int occa, occb, virta, virtb;
      readOccN(occa, occb, virta, virtb);
      int aoCount = virta + occa;
      std::vector<double> fullAO(aoCount*aoCount*aoCount*aoCount);

      readAOchemical (aoCount, fullAO);

      std::vector<int> compDim(4);
      compDim[0] = aoCount;
      compDim[1] = aoCount;
      compDim[2] = aoCount;
      compDim[3] = aoCount;

      int aoCount2= aoCount * aoCount;

      std::vector<double> U(1);
      std::vector<double> VT(1);
      std::vector<double> S(aoCount2);
      std::vector<double> work(1);
      int lwork = -1;

      Lapack<double>::gesvd('n', 'n', aoCount2, aoCount2, &fullAO[0], aoCount2, &S[0], &U[0], 1, &VT[0], 1, &work[0], lwork);

      lwork = work[0];
      work.resize(lwork);

      Lapack<double>::gesvd('n', 'n', aoCount2, aoCount2, &fullAO[0], aoCount2, &S[0], &U[0], 1, &VT[0], 1, &work[0], lwork);

      //std::cout<<"S :"<<S<<std::endl;
      int r = 0;
      for(int i=0; i<aoCount2; i++){
         if(S[i] > 1e-14){
            r++;
         }
      }
      std::cout<<"rang = "<<r<<std::endl;
}

void MPSreadWrite (){
   double eps = 1e-2;

   int occ;
   int virt;
   int occB;
   int virtB;
   readOccN(occ, occB, virt, virtB);

   int aoCount  = occ + virt;
   std::vector<double> fullao;

   readAOchemical(aoCount, fullao);

   std::vector<int> compDim(4);
   compDim[0] = aoCount;
   compDim[1] = aoCount;
   compDim[2] = aoCount;
   compDim[3] = aoCount;

   TensorCalculus::MPSRepresentation<double> mpsAO(fullao, compDim, eps);
   std::cout<<"truncate AO"<<std::endl;
   std::cout<<"eps : "<<eps<<"  summations : "<<mpsAO.getSummations()<<std::endl;
   std::cout<<"write to disk"<<std::endl;
   mpsAO.write2disk("mpsChemAO.MPSten");
   std::cout<<"read from disk"<<std::endl;
   TensorCalculus::MPSRepresentation<double> readMPS;
   readMPS.read("mpsChemAO.MPSten");

   std::cout<<"norm write-read : "<<l2_norm(mpsAO.evaluate()-readMPS.evaluate())<<std::endl;
}

void MPSwrite (){
      double eps = 1e-12;

      int occ;
      int virt;
      int occB;
      int virtB;
      readOccN(occ, occB, virt, virtB);

      int aoCount  = occ + virt;
      std::vector<double> fullao;

      readAOchemical(aoCount, fullao);

      std::vector<int> compDim(4);
      compDim[0] = aoCount;
      compDim[1] = aoCount;
      compDim[2] = aoCount;
      compDim[3] = aoCount;

      TensorCalculus::MPSRepresentation<double> mpsAO(fullao, compDim, eps);
      std::cout<<"truncate AO"<<std::endl;
      std::cout<<"eps : "<<eps<<"  summations : "<<mpsAO.getSummations()<<std::endl;
      std::cout<<"write to disk"<<std::endl;
      mpsAO.write2disk("mpsChemAO.MPSten");
}

void testFullMP2 (){
   int occB = 0;
   int virtB = 0;
   int occ = 0;
   int virt = 0;
   readOccN(occ, occB, virt, virtB);

   std::vector<double> fullao;

   int aoCount = occ + virt;
   readAO(aoCount, fullao);


   AOMOtransformation transform (virt, occ);
   transform.fullTensorAO2fullTensorMO();


   std::vector<double> fullMO(transform.getFullTensorMO());
   int size = virt*virt*occ*occ;
   std::vector<double> v_vv_oo(size);
   std::vector<double> v_vv_oo_simple(size);


   std::vector<double> spinAdMO(fullMO);
   spinAdaptedAO(aoCount, fullMO, spinAdMO);

   int index1;
   int index2;


   for(int i=0; i<occ; i++){
     for(int j=0; j<occ; j++){
        for(int k=0; k<virt; k++){
           for(int l=0; l<virt; l++){
              index1 = l + virt * (k + virt * (j + occ * i));
              index2 = l + occ + aoCount * (k + occ + aoCount * (j + aoCount * i));
              v_vv_oo_simple[index1] = fullMO[index2];
              v_vv_oo[index1] = spinAdMO[index2];
           }
        }
     }
  }
   std::cout<<"||v - v_spinAd|| = "<<l2_norm(v_vv_oo_simple - v_vv_oo)<<std::endl;
   //std::cout<<v_vv_oo<<std::endl;
   Amplitudes mp2(v_vv_oo,v_vv_oo_simple);

   std::cout<<"compute mp2 energy ..."<<std::endl;
   double energie = mp2.MP2Energie();
   std::cout<<"MP2 energie : "<<energie<<std::endl;
}


void testHadamard (){
   int d = 4;
   int n = 3;
   int j = 3;
   int k = 2;
   
   TensorCalculus::Random<double> random;

   std::cout<<"d ="<<d<<std::endl;
   std::cout<<"n ="<<n<<std::endl;
   std::cout<<"j ="<<j<<std::endl;
   std::cout<<"k ="<<k<<std::endl;
   std::vector<int> compDim(d);
   for(int mu=0; mu<d; mu++){
      compDim[mu] = n+mu;
   }
   std::vector<int> summations(d-1);
   for(int mu=0; mu<d-1; mu++){
      summations[mu] = j+mu;
   }
   std::vector<std::vector<double> > CPv(d);
   for(int i=0; i<d; i++){
      CPv[i].resize(compDim[i]*k);
      for(int o=0; o< compDim[i]*k; o++){
         CPv[i][o] = random();
      }
   }
   std::cout<<"set cp ..."<<std::endl;
   CPTensorRepresentation<double> cp(k, compDim, CPv);
   std::cout<<"set cp ... finish"<<std::endl;
   std::vector<std::vector<double> > MPSv(d);
   MPSv[0].resize(summations[0]*compDim[0]);
   MPSv[d-1].resize(summations[d-2]*compDim[d-1]);
   for(int i=0; i<summations[0]*compDim[0]; i++){
      MPSv[0][i]= random();
   }
   for(int i = 0; i<summations[d-2]*compDim[d-1]; i++){
      MPSv[d-1][i]= random();
   }
   for(int mu=1; mu<d-1; mu++){
      MPSv[mu].resize(summations[mu-1]*summations[mu]*compDim[mu]);
      for(int i=0; i<summations[mu-1]*summations[mu]*compDim[mu]; i++){
         MPSv[mu][i]= random();
      }
   }
   std::cout<<"set mps ..."<<std::endl;
   MPSRepresentation<double> mps (summations, MPSv, compDim);
   MPSRepresentation<double> mps2(mps);
   mps2.setHadamardProduct(cp);
   /*
   std::cout<<"cp : "<<std::endl;
   std::cout<<cp.evaluate()<<std::endl;
   std::cout<<"mps : "<<std::endl;
   std::cout<<mps.evaluate()<<std::endl;
   std::cout<<"mps : "<<std::endl;
   std::cout<<mps2.evaluate()<<std::endl;
   */
   std::vector<double> cpFull(cp.evaluate());
   std::vector<double> mpsFull(mps.evaluate());
   std::vector<double> hadFull(mpsFull.size());
   for(unsigned int i=0; i<cpFull.size(); i++){
      hadFull[i]=cpFull[i]*mpsFull[i];
   }
   std::cout<<"norm(hadamard) = "<<l2_norm(hadFull)<<std::endl;
   std::cout<<"dist = "<<l2_norm(hadFull-mps2.evaluate())<<std::endl;
}

void testAdd(){
      TensorCalculus::Random<double> random;
  
      int d = 4;
      int n = 7;
      int j = 4;
      int k = 3;

      std::cout<<"d ="<<d<<std::endl;
      std::cout<<"n ="<<n<<std::endl;
      std::cout<<"j ="<<j<<std::endl;
      std::cout<<"k ="<<k<<std::endl;
      std::vector<int> compDim(d);
      for(int mu=0; mu<d; mu++){
         compDim[mu] = n;
      }
      std::vector<int> summations(d-1);
      std::vector<int> summations2(d-1);
      for(int mu=0; mu<d-1; mu++){
         summations[mu] = j;
         summations2[mu] = k;
      }
      std::vector<std::vector<double> > MPSv(d);
      std::vector<std::vector<double> > MPSv2(d);
      MPSv[0].resize(j*n);
      MPSv[d-1].resize(j*n);
      MPSv2[0].resize(k*n);
      MPSv2[d-1].resize(k*n);
      for(int i=0; i<j*n; i++){
         MPSv[0][i]= random();
         MPSv[d-1][i]= random();
      }
      for(int i=0; i<k*n; i++){
         MPSv2[0][i]= random();
         MPSv2[d-1][i]= random();
      }
      for(int mu=1; mu<d-1; mu++){
         MPSv[mu].resize(j*j*n);
         MPSv2[mu].resize(k*k*n);
         for(int i=0; i<j*j*n; i++){
            MPSv[mu][i]= random();
         }
         for(int i=0; i<k*k*n; i++){
            MPSv2[mu][i]= random();
         }
      }
      std::cout<<"set mps ..."<<std::endl;
      MPSRepresentation<double> mps (summations, MPSv, compDim);
      MPSRepresentation<double> mps2 (summations2, MPSv2, compDim);


      std::cout<<"summations mps : "<<mps.getSummations()<<std::endl;
      std::cout<<"summations mps2: "<<mps2.getSummations()<<std::endl;

      std::vector<double> sum(mps.evaluate()+mps2.evaluate());
      mps2.add(mps);

      std::cout<<"summations mps2: "<<mps2.getSummations()<<std::endl;
      std::cout<<"norm summe = "<<l2_norm(sum)<<std::endl;
      std::cout<<"dist       = "<<l2_norm(mps2.evaluate()-sum)<<std::endl;
}

void testEdgeAdd(){
   TensorCalculus::Random<double> random;
  
   int d = 4;
   int n = 7;
   int i = 2;
   int j = 4;
   int k = 3;

   int j1 = d-3;

   std::cout<<"d ="<<d<<std::endl;
   std::cout<<"n ="<<n<<std::endl;
   std::cout<<"i ="<<i<<std::endl;
   std::cout<<"j ="<<j<<std::endl;
   std::cout<<"k ="<<k<<std::endl;
   std::cout<<"j1 ="<<j1<<std::endl;
   std::vector<int> compDim(d);
   for(int mu=0; mu<d; mu++){
      compDim[mu] = n;
   }
   std::vector<int> summations(d-1);
   std::vector<int> summations2(d-1);
   for(int mu=0; mu<d-1; mu++){
      if(j1 != mu){
         summations[mu] = i;
         summations2[mu] = i;
      } else {
         summations[mu] = j;
         summations2[mu] = k;
      }
   }
   std::vector<std::vector<double> > MPSv(d);
   std::vector<std::vector<double> > MPSv2(d);
   MPSv[0].resize(summations[0]*n);
   MPSv[d-1].resize(summations[d-2]*n);
   MPSv2[0].resize(summations2[0]*n);
   MPSv2[d-1].resize(summations2[d-2]*n);
   for(unsigned int i=0; i<MPSv[0].size(); i++){
      MPSv[0][i]= random();
   }
   for(unsigned int i=0; i<MPSv[d-1].size(); i++){
      MPSv[d-1][i]= random();
   }
   for(unsigned int i=0; i<MPSv2[0].size(); i++){
      MPSv2[0][i]= random();
   }
   for(unsigned int i=0; i<MPSv2[d-1].size(); i++){
      MPSv2[d-1][i]= random();
   }
   for(int mu=1; mu<d-1; mu++){
      MPSv[mu].resize(summations[mu-1]*summations[mu]*n);
      MPSv2[mu].resize(summations2[mu-1]*summations2[mu]*n);
      for(unsigned int i=0; i<MPSv[mu].size(); i++){
         MPSv[mu][i]= random();
      }
      for(unsigned int i=0; i<MPSv2[mu].size(); i++){
         MPSv2[mu][i]= random();
      }
   }
   for(int mu=0; mu<d; mu++){
      if(mu != j1 && mu != j1+1){
         MPSv[mu] = MPSv2[mu];
      }
   }
   std::cout<<"set mps ..."<<std::endl;
   MPSRepresentation<double> mps (summations, MPSv, compDim);
   std::cout<<"set mps2 ..."<<std::endl;
   MPSRepresentation<double> mps2 (summations2, MPSv2, compDim);
   MPSRepresentation<double> AddMps (mps);
   std::vector<double> sum (mps.evaluate()+mps2.evaluate());
   AddMps.add(mps2);
   mps.edgeAdd(j1,mps2);

   std::cout<<"mps.summations  = "<<mps.getSummations()<<std::endl;
   std::cout<<"Add.summations  = "<<AddMps.getSummations()<<std::endl;
   std::cout<<"dist add,mps = "<<l2_norm(AddMps.evaluate()-mps.evaluate())<<std::endl;
   std::cout<<"dist add,sum = "<<l2_norm(AddMps.evaluate()-sum)<<std::endl;
   std::cout<<"dist sum,mps = "<<l2_norm(sum-mps.evaluate())<<std::endl;

}


void mpsTruncTest(){
   std::string file = "compute/AOchem2Vaibj.mps.1e-12";
   std::ifstream fin(&file[0]);
   if(fin.fail()){
      throw std::invalid_argument("no input file");
   }
   MPSRepresentation<double> in(&file[0]);
   MPSRepresentation<double> in2(&file[0]);
   std::cout<<"#############################"<<std::endl;
   in.truncateMPS(1e-2,'e',true);
   std::cout<<"#############################"<<std::endl;
   std::vector<double> in2Eval (in2.evaluate());
   double norm = l2_norm(in.evaluate() - in2Eval);
   std::cout<<"dist : "<<norm<<std::endl;
   std::cout<<"rel dist : "<<norm / l2_norm(in2Eval)<<std::endl;
   std::cout<<"summations alt: "<<in2.getSummations()<<std::endl;
   std::cout<<"summations neu: "<<in.getSummations()<<std::endl;
   /*
   in.truncateMPS(0.5,true);
   std::cout<<"#############################"<<std::endl;

   norm = l2_norm(in.evaluate() - in2Eval);
   std::cout<<"dist : "<<norm<<std::endl;
   std::cout<<"rel dist : "<<norm / l2_norm(in2Eval)<<std::endl;
   std::cout<<"summations alt: "<<in2.getSummations()<<std::endl;
   std::cout<<"summations neu: "<<in.getSummations()<<std::endl;
   */
}

void testFullMP2_b(){
   FullTensor<double> v_abij("v_abij.ften");
   FullTensor<double> t_abij("t_abij.ten");

   std::cout<<"MP2 Energie = "<< innerProduct(v_abij.getV(), t_abij.getV())<<std::endl;
}

void makeMPSforStefan(){
   std::vector<double> aoC;
   std::vector<double> aoP;
   readAO(aoP);
   readAOchemical(aoC);
   FullTensor<double> chemAO(makeFullTensorAO(aoC));
   FullTensor<double> physAO(makeFullTensorAO(aoP));

   MPSRepresentation<double> chemMPS(chemAO, 1e-12);
   MPSRepresentation<double> physMPS(physAO, 1e-12);

   chemMPS.writeToFile("chemMPS.ten");
   physMPS.writeToFile("physMPS.ten");
}

void testVabij2viabj(){
   FullTensor<double> ft (Vabij2Viabj());

   ft.write2disk("v2v_iabj.ften");
}
