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
#include "CC/C4Interface.hpp"
#include "BlasInterface.hpp"
#include "Vector/VectorOperators.hpp"


namespace TensorCalculus {

   void AO2MO (const TensorRepresentation<double> &in, TensorRepresentation<double> &virt, TensorRepresentation<double> &occ){
      //virt and occ must have been initialized properly
      //
      int occCount, virtCount;
      readOccN(occCount, occCount, virtCount, virtCount);
      int aoCount = occCount + virtCount;
      std::vector<double> virtAOMO, occAOMO;
      readTransformationMatrix(virtCount, occCount, virtAOMO, occAOMO);
      int prodOfRank;
      for(int mu=0; mu<4; mu++){
         if(mu == 0){
            prodOfRank = in.getSummation(0);
         }else if(mu == 3){
            prodOfRank = in.getSummation(2);
         }else {
            prodOfRank = in.getSummation(mu)*in.getSummation(mu-1);
         }
         Blas<double>::gemm('n','n', virtCount, prodOfRank, aoCount, 1.0, &virtAOMO[0],
               virtCount, &in.getV()[mu][0], aoCount, 0.0, &virt(mu,0), virtCount);
         Blas<double>::gemm('n','n', occCount, prodOfRank, aoCount, 1.0, &occAOMO[0],
               occCount, &in.getV()[mu][0], aoCount, 0.0, &occ(mu,0), occCount);
      }
   }

   void makeV (const std::vector<char> occState, TensorRepresentation<double> &V_mnop,
         const TensorRepresentation<double> &virt, const TensorRepresentation<double> &occ){
      std::vector<int> copyLength(4);
      copyLength[0] = V_mnop.getSummation(0);
      copyLength[1] = V_mnop.getSummation(0) * V_mnop.getSummation(1);
      copyLength[2] = V_mnop.getSummation(1) * V_mnop.getSummation(2);
      copyLength[3] = V_mnop.getSummation(2);
      for(int mu=0; mu<4; mu++){
         copyLength[mu] *= V_mnop.getComponentDimension(mu);
         if (occState[mu] == 'v'){
            //std::cout<<"mu = "<<mu<<" : "<<occState[mu]<<std::endl;
            Blas<double>::copy(copyLength[mu], &virt.getV()[mu][0], 1, &V_mnop(mu,0), 1);
         }else if (occState[mu] == 'o'){
            //std::cout<<"mu = "<<mu<<" : "<<occState[mu]<<std::endl;
            Blas<double>::copy(copyLength[mu], &occ.getV()[mu][0], 1, &V_mnop(mu,0), 1);
         }else{
            throw std::invalid_argument("undefined occumpation State");
         }
      }
   }

   MPSRepresentation<double> makeV (const std::vector<char> &occState, const MPSRepresentation<double> &AO){
      int o, v;
      readOccN(o, o, v, v);
      std::vector<int> virtCount(4);
      std::vector<int> occCount(4);
      virtCount[0] = v;
      virtCount[1] = v;
      virtCount[2] = v;
      virtCount[3] = v;
      occCount[0] = o;
      occCount[1] = o;
      occCount[2] = o;
      occCount[3] = o;
      std::vector<int> compDim(4);
      for(int mu=0; mu < 4; mu++){
         if(occState[mu] == 'v'){
            compDim[mu] = v;
         }else if (occState[mu] == 'o'){
            compDim[mu] = o;
         }else{
            throw std::invalid_argument("invalid occumpation state");
         }
      }

      MPSRepresentation<double> V_mnop(AO.getSummations(), compDim);

      MPSRepresentation<double> virt (AO.getSummations(), virtCount);
      MPSRepresentation<double> occ (AO.getSummations(), occCount);
      AO2MO(AO, virt, occ);
      makeV(occState, V_mnop, virt, occ);
      return(V_mnop);
   }

   void AOMOtransformation::init(const int virt, const int occ){
      //todo: init
      (*this).virtCount = virt;
      (*this).occCount = occ;
      (*this).aoCount = virtCount+occCount;

      readAO(aoCount, fullTensorAO);
      fullTensorMO.resize(fullTensorAO.size());

      std::vector<int> componentDimensions(4);
      componentDimensions[0] = aoCount;
      componentDimensions[1] = aoCount;
      componentDimensions[2] = aoCount;
      componentDimensions[3] = aoCount;

      //tested
      AO = MPSRepresentation<double> (fullTensorAO, componentDimensions, 0.0); //todo: wie ao definieren?
      //test
      //using namespace VectorOperators;
      //std::cout<<"init transform : norm : "<<l2_norm(AO.evaluate()-fullAO)<<std::endl;

      MOv.resize(4);
      for (int mu=0; mu<4; mu++){
         MOv[mu].resize(AO.getV()[mu].size());
      }

      readTransformationMatrix(virtCount, occCount, virtAOMO, occAOMO);
      readTransformationMatrix(aoCount, fullTransposeAOMO);

      prodOfRank.resize(4);
      prodOfRank[0]=AO.getSummation(0);
      prodOfRank[1]=AO.getSummation(0) * AO.getSummation(1);
      prodOfRank[2]=AO.getSummation(1) * AO.getSummation(2);
      prodOfRank[3]=AO.getSummation(2);

      //using namespace VectorOperators;
      //std::cout<<"prodOfRank : "<<prodOfRank<<std::endl;

   }

   void AOMOtransformation::AO2MO (){
      char aomo = 't'; //standard 't'
      for (int mu=0; mu<4; mu++){

         /*//AO in physical notation
         if(mu == 2){
            aomo = 'n';
         }
         */

         /*//AO in chemical notation --> achtung t und n richtig setzen
         if(mu == 1 || mu == 3){
            aomo = 't';
         }else{
            aomo = 'n';
         }
         */
         Blas<double>::gemm(aomo,'n', aoCount, prodOfRank[mu], aoCount, 1.0, &fullTransposeAOMO[0],
                            aoCount, &AO.getV()[mu][0], aoCount, 0.0, &MOv[mu][0], aoCount);

      }
      /*
      std::vector< std::vector <double> > temp;
      temp = MOv;
      for (int mu=0; mu<4; mu++){
         Blas<double>::gemm('n','n', aoCount, prodOfRank[mu], aoCount, 1.0, &fullTransposeAOMO[0],
                            aoCount, &temp[mu][0], aoCount, 0.0, &MOv[mu][0], aoCount);

      }
      */
   }

   MPSRepresentation<double> AOMOtransformation::getVabij (){
      //ab ... virtual Orbitals
      //ij ... occupied Orbitals
      std::vector<int> componentDimensions(4);
      componentDimensions[0] = virtCount;
      componentDimensions[1] = virtCount;
      componentDimensions[2] = occCount;
      componentDimensions[3] = occCount;

      std::vector< std::vector<double> > v(4);

      int n      = virtCount;
      int start  = occCount;

      for(int mu=0; mu<4; mu++){
         if(mu==2){
            n = occCount;
            start = 0.0;
         }
         v[mu].resize(n*prodOfRank[mu]);

//#pragma omp parallel for --> fehler!
         for(int j=0; j < prodOfRank[mu]; j++){
            //std::cout<<"mu = "<< mu <<"  j = "<< j <<std::endl;
            Blas<double>::copy(n, &MOv[mu][start+j*aoCount], 1, &v[mu][j*n], 1);
         }
      }
      std::cout<<"end copy"<<std::endl;
      return(MPSRepresentation<double> (AO.getSummations(), v, componentDimensions));
   } /*end getVabij ()*/

   MPSRepresentation<double> AOMOtransformation::AO2Vabij (){

      std::vector< std::vector<double> > v(4);

      std::vector<int> componentDimensions(4);
      componentDimensions[0] = virtCount;
      componentDimensions[1] = virtCount;
      componentDimensions[2] = occCount;
      componentDimensions[3] = occCount;

      int  m = virtCount;
      double *a = &virtAOMO[0];

      for(int mu=0; mu<4; mu++){
         std::cout<<"mu : "<<mu<<std::endl;
         if(mu==2){
            m = occCount;
            a = &occAOMO[0];
         }
         v[mu].resize(m*prodOfRank[mu]);
         Blas<double>::gemm('n', 'n', m, prodOfRank[mu], aoCount, 1.0, a, m,
                           &AO.getV()[mu][0], aoCount, 0.0, &v[mu][0], m);
      }

      return(MPSRepresentation<double> (AO.getSummations(), v, componentDimensions));
   }/*end AO2Vabij*/

   MPSRepresentation<double> AOMOtransformation::getFullMO (){
      return(MPSRepresentation<double> (AO.getSummations(), MOv, AO.getComponentDimensions()));
   }

   void AOMOtransformation::fullTensorAO2fullTensorMO (){

      int aoCount3 = aoCount * aoCount * aoCount;

      std::vector<double> temp(fullTensorAO);

      for(int mu=0; mu<4; mu++){

         //<mu | nu, lam, sig> --> <a, nu, lam, sig>
         //<nu | lam, sig, a>  --> <b, lam, sig, a>
         //<lam| sig, a, b>    --> <c, sig, a,b>
         //<sig| a, b, c>      --> <d, a, b, c>
         Blas<double>::gemm('t', 'n', aoCount, aoCount3, aoCount, 1.0, &fullTransposeAOMO[0],
                        aoCount, &temp[0], aoCount, 0.0, &fullTensorMO[0], aoCount);
         //<a | nu, lam, sig>  --> <nu, lam, sig, a>
         //<b | lam, sig, a>   --> <lam, sig, a, b>
         //<c | sig, a,b>      --> <sig, a, b, c>
         //<d | a, b, c>       --> <a, b, c, d>
         transpose(aoCount, aoCount3, &fullTensorMO[0], &temp[0]);
      }
      fullTensorMO = temp;
   } //end fullTensorAO2fullTensorMO

/*
void partialAOMO (const int aoCount, const int moCount, const int productOfRanks,
                  const std::vector<double> &AOMO, const std::vector<double> &AO,
                  std::vector<double> &MO){

   /////////////////////////////////////////////////////////////////////////////////////////////////////
   //input : aoCount      : number of atomic Orbitals (AO)
   //        moCount      : number of molecular Orbitals (MO)
   //                     : MO can be occupied (occ), virtual (virt), or full (occ+virt)
   //        productOfRank: product of all ranks of the mu-th edge
   //        AOMO         : transformationMatrix, can be virtAOMO, occAOMO (read by
   //                       readTransformationMAtrix)or fullAOMO (NEWMOS)
   //        AO           : one index (mu) of the decomposed two-elektron-integrals
   //                       AO-basis (TWOINT)
   //                       v[mu]
   //output: MO           : one index (mu) of the decomposed two-elektron-integrals transformed to
   //                       the MO-basis
   //                       new v[mu]
   /////////////////////////////////////////////////////////////////////////////////////////////////////

   if(MO.size() < moCount * productOfRanks){
      MO.resize(moCount * productOfRanks);
   }

   // compute AOMO * AO = MO ,
   //AOMO = Matrix(moCount x aoCount)
   //AO   = Matrix(aoCount x productOfRanks)
   //MO   = Matrix(moCount x productOfRanks)

   const int &m = moCount;
   const int &n = productOfRanks;
   const int &k = aoCount;
   const int &lda = m;
   const int &ldb = k;
   const int &ldc = m;

   TensorCalculus::Blas<double>::gemm('n','n', m, n, k, 1.0, &AOMO[0], lda, &AO[0], ldb, 0.0, &MO[0], ldc);

}

//void AOMOtransformation (const TensorCalculus::MPSRepresentation<double> &AO, int moCount, std::vector<double> &MO){

*/
} /*end namespace TensorCalculus */
