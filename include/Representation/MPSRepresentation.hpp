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

#ifndef __MPSREPRESENTATION_HPP
#define __MPSREPRESENTATION_HPP

#include <iostream>
#include "Matrix/MatrixOperators.hpp"
#include "Representation/TensorRepresentation.hpp"
#include "Representation/CPTensorRepresentation.hpp"
#include "Tensor/FullTensor.hpp"
#include "LapackInterface2.hpp"
#include "time.h"
#include "Representation/MPSDecomposition.hpp"

namespace TensorCalculus {

template<typename T> class MPSRepresentation: public TensorRepresentation<T> {
private:
	std::vector<std::vector<int> > setIncidenceMatrix(const int d) {
		std::vector<std::vector<int> > incidenceMatrix(d);
		incidenceMatrix[0].resize(1);
		incidenceMatrix[0][0] = 0;
		incidenceMatrix[d - 1].resize(1);
		incidenceMatrix[d - 1][0] = d - 2;

		for (int n = 1; n < d - 1; n++) {
			incidenceMatrix[n].resize(2);
			incidenceMatrix[n][0] = n - 1;
			incidenceMatrix[n][1] = n;
		}
		return (incidenceMatrix);
	}

	void MPSinit(const std::vector<T> &fullTensor, const std::vector<int> &componentDimensions,
	                     const T eps = 0.0){
	   T relEps = eps * l2_norm(fullTensor);
	   TensorCalculus::MPSDecomposition<T> decomposition(fullTensor, componentDimensions);
     decomposition.computeVidalDecomposition(relEps);
     init(decomposition.getSummations(), decomposition.getV(), componentDimensions,
           setIncidenceMatrix(componentDimensions.size()));
	}


	void mpsqlf(std::vector <std::vector <T> > &vTemp){

	   //computes QL-factorisation using MPSTensor-format.
	   //A = A(i1...id) = sum over j: L(i1|j)*Q(i2...id|j) and
	   //Q(i2...id|j) = sum over k2...k{d-1}: G2(i2, j, k2) *...*G{d-1}(i{d-1}, k{d-2}, k{d-1})
     //                                      *Gd(id, k{d-1})
	   //
	   //Q satisfies the following orthogonality property:
	   //<Q,Q> = sum over i2...id: Q(i2...id|j)Q(i2...id|j')
	   //      = delta(j,j')
	   //
	   //vTemp saves elementary reflectors of T1...Td-1,
	   // Tk is a matrix(j'k : ik, j'k-1) and coresponds to the Tensor Gk (ik, j'k-1, j'k)
	   //
	   //todo: rückgabe von v ermöglichen?
	   //--> bool computeV = false übergeben und an den test unten anhängen

	    unsigned int d = (*this).componentDimensions.size();

	    if(vTemp.size() < d){
	       vTemp.resize(d);
	    }
      for(unsigned int mu=0; mu<d; mu++){
         if(vTemp[mu].size() < (*this).v[mu].size()){
            vTemp[mu].resize((*this).v[mu].size());
         }
      }

	    vTemp[d-1]=(*this).v[d-1]; //startvalue for A;

	    //parameters for geqrf
	    int m1; //--> zum Array machen
	    int n1;
	    unsigned int min;
	    int &lda1 = m1;
	    std::vector< std::vector <double> > tau(d);    //saves tau from dgeqrf first,
                                                    //overwriten by s from dgesvd

	    //parameters for trmm
	    int m2;
	    int &n2    = n1;
      int &ldb   = m2;

      //workspace for Lapack
	    std::vector<double> work(1);
	    const int getWork = -1;
	    int lwork = -1;

	    //left to right sweep
	    //std::cout<<"left to right sweep"<<std::endl;
	    for(unsigned int mu = d-1; mu > 0; mu--){
	       //std::cout<<"mu = "<<mu<<std::endl;
	       //what happens if m1 < n1? --> fixed
	       if(mu==d-1){
	          m1 = (*this).componentDimensions[mu];
	       }else{
	          m1 = (*this).componentDimensions[mu]*(*this).summations[mu];
         }
	       n1 = (*this).summations[mu-1]; //caution: this is right, cause we reshape v[mu] at the end of the for-loop
	       //lda = m;
	       min = std::min(m1,n1);

	       //std::cout<<"m1 = "<<m1<<std::endl;
	       //std::cout<<"n1 = "<<n1<<std::endl;

	       // vTemp is copy of v, from the second loop we need vTemp because of transposition/reshape
	       std::vector<T> &A = vTemp[mu];
	       //std::vector<T> &B = vTemp[mu-1];
	       //std::vector<T> &A = v[mu];
	       // first compute v[mu-1] * R-factor, then transpose/reshape to vTemp[mu-1]
	       std::vector<T> &B = (*this).v[mu-1];
	       if(tau[mu].size()<min){
	          tau[mu].resize(min);
	       }
	       //set workspace for QR-fact.
	       Lapack<T>::geqrf(m1, n1, &A[0], lda1, &tau[mu][0], &work[0], getWork);
	       if(lwork < work[0]){
	          lwork = work[0];
	          work.resize(lwork);
	       }
	       //compute QR fact.
	       Lapack<T>::geqrf(m1, n1, &A[0], lda1, &tau[mu][0], &work[0], lwork);
	       //test; Q als reflektor stehen lassen

	       //compute v[mu-1] = v[mu-1]*R
	       if(mu != 1){
	          m2 = (*this).componentDimensions[mu-1]*(*this).summations[mu-2];
	       }else{
	          m2 = (*this).componentDimensions[mu-1]; //mu=1: there is no summations[-1] for v[0]
	       }
	       //n2 = n;

	       //need lda1 because because its the same A as for QR-factorisation
	       if(m1 >= n1){
	          Blas<T>::trmm('r','u','t','n',m2, n2, 1.0, &A[0], lda1/*lda2*/, &B[0], ldb);
	       } else {
	         // std::cout<<"m1 < n1"<<std::endl;
	          Blas<T>::trmm('r','u','t','n',m2, m1, 1.0, &A[0], lda1/*lda2*/, &B[0], ldb);

	          //there is no memory conflict::
	          //B = (B1 | B2)
	          // alpha * B2 * A + beta B1 --> B1
	          Blas<T>::gemm('n','t',m2, m1, (n1-m1), 1.0, &B[m2*m1], ldb, &A[m1*m1], lda1, 1.0, &B[0], ldb);
	          B.resize(m2*m1);
	          (*this).summations[mu-1] = m1;
	          n1 = m1;
	       }

	       (*this).v[mu] = giveQ(m1, n1, &A[0], lda1, &tau[mu][0]);

	       //reshape (i_mu, j_mu-1, j_mu) --> (j_mu | i_mu, j_mu-1)
	       //--> for mu=1 auch machen?
	       if(mu != 1){
	          transpose(m2, n2, &B[0], &vTemp[mu-1][0]);
	       }
	    }
	}

	void mpssvd(const T eps, const char epsWeight ='e', const bool display = false){
	   T relEps = eps * l2norm((*this));
	   std::vector<T> epsMu;
	   setEps((*this).d, epsMu, relEps, epsWeight, display);

	   T tempSize = *max_element((*this).summations.begin(), (*this).summations.end());
	   tempSize *= tempSize;
	   tempSize *= *max_element((*this).componentDimensions.begin(), (*this).componentDimensions.end());
	   std::vector<T> temp(tempSize);
	   int vSize;

	   const char jobu  = 'S';
	   const char jobvt = 'S';
	   int d = (*this).d;
	   int m;
	   int n;
	   int min;
	   int rTrunc;
	   //std::vector<T> U(tempSize);
	   //std::vector<T> VT(tempSize);
	   std::vector<T> U;
	   std::vector<T> VT;
	   std::vector<T> S;
	   int lWork = -1;
	   std::vector<T> work(1);
	   for (int mu = 0; mu<d-1; mu++){
	      //std::cout<<"mu = "<<mu<<std::endl;
	      if(mu == 0){
	         // truncation checked
	         m = (*this).componentDimensions[mu];
	         n = (*this).summations[mu];
	         min = std::min(m,n);

	         setSizesForSVD(jobu, jobvt, m, n, (*this).v[mu], m, S, U, m, VT, min, work, lWork);
           Lapack<T>::gesvd(jobu, jobvt, m, n, &(*this).v[mu][0], m, &S[0], &U[0], m, &VT[0], min, &work[0], lWork);
           //std::cout<<" S: "<< S <<std::endl;
           rTrunc = min;
           evaluateEps(epsMu, S, mu, rTrunc);
           //rTrunc = min;

           //std::cout<<"rTrunc : "<<rTrunc<<std::endl;

           (*this).v[mu].resize(m*rTrunc);
           // change U to v[0] with jobu = 'S';
           memcpy(&(*this).v[mu][0], &U[0], m * rTrunc * (*this).TSIZE);

           for(int i=0; i< rTrunc; i++){
              Blas<T>::scal(n, S[i], &VT[i], min);
           }
           Blas<T>::gemm('n', 't', (*this).summations[mu+1]*(*this).componentDimensions[mu+1], rTrunc, (*this).summations[mu], 1.0, &(*this).v[mu+1][0], (*this).summations[mu+1]*(*this).componentDimensions[mu+1], &VT[0], min, 0.0, &temp[0], (*this).summations[mu+1]*(*this).componentDimensions[mu+1]);
           vSize = (*this).summations[mu+1]*(*this).componentDimensions[mu+1]*rTrunc;
	      } else {
	         //truncation checked
	         m = (*this).summations[mu];
	         n = (*this).componentDimensions[mu] * (*this).summations[mu-1];
	         min = std::min(m,n);

	         setSizesForSVD(jobu, jobvt, m, n, (*this).v[mu], m, S, U, m, VT, min, work, lWork);
	         Lapack<T>::gesvd(jobu, jobvt, m, n, &(*this).v[mu][0], m, &S[0], &U[0], m, &VT[0], min, &work[0], lWork);
           //std::cout<<" S: "<< S <<std::endl;
	         rTrunc = min;
	         evaluateEps(epsMu, S, mu, rTrunc);

	         //std::cout<<"rTrunc : "<<rTrunc<<std::endl;

	         (*this).v[mu].resize(n*rTrunc);
	         transpose(rTrunc, n, &VT[0], min, &(*this).v[mu][0], n);

	         for(int i=0; i< rTrunc; i++){
	            Blas<T>::scal(m, S[i], &U[i*m], 1);
           }

	         if (mu != d-2){
              Blas<T>::gemm('n', 'n', (*this).summations[mu+1]*(*this).componentDimensions[mu+1], rTrunc, (*this).summations[mu], 1.0, &(*this).v[mu+1][0], (*this).summations[mu+1]*(*this).componentDimensions[mu+1], &U[0], (*this).summations[mu], 0.0, &temp[0], (*this).summations[mu+1]*(*this).componentDimensions[mu+1]);
              vSize = (*this).summations[mu+1]*(*this).componentDimensions[mu+1]*rTrunc;
	         } else {
	            // mu == d-2
              std::vector<T> v1((*this).componentDimensions[mu+1]*rTrunc);
              Blas<T>::gemm('n', 'n', (*this).componentDimensions[mu+1], rTrunc, (*this).summations[mu], 1.0, &(*this).v[mu+1][0], (*this).componentDimensions[mu+1], &U[0], (*this).summations[mu], 0.0, &temp[0], (*this).componentDimensions[mu+1]);
              //test
              vSize = (*this).componentDimensions[mu+1]*rTrunc;
	         }
	      }
	      (*this).v[mu+1].resize(vSize);
	      memcpy(&(*this).v[mu+1][0], &temp[0], vSize*(*this).TSIZE);
	      (*this).summations[mu]=rTrunc;
	   }
	} //end   void mpssvd(std::vector <std::vector <T> > &vTemp)

	MPSRepresentation<T>& partialHadamardProduct (const int mu, const CPTensorRepresentation<T> &CP){
	     int k = CP.getSummation(0);
	     int length = (*this).componentDimensions[mu];
	     int prodOfSummations;
	     int sizeReduced = length;
	     int sizeFull;
	     if(mu == 0){
	        prodOfSummations = (*this).summations[0];
	        sizeReduced *= prodOfSummations * k;
	        sizeFull = 0;
	     }else if(mu == (*this).d-1){
	        prodOfSummations = (*this).summations[(*this).d-2];
	        sizeReduced *= prodOfSummations * k;
	        sizeFull = 0;
	     }else {
	        prodOfSummations = (*this).summations[mu-1] * (*this).summations[mu];
	        sizeReduced *= prodOfSummations * k;
	        sizeFull = sizeReduced * k;
	     }
	     std::vector<T> vReduced(sizeReduced);
	     std::vector<T> vFull(sizeFull);
	     for(int i=0; i< (*this).componentDimensions[mu]; i++){
	        /*
	        //falsch inkrementiert
	        Blas<T>::ger(prodOfSummations, CP.getSummation(0), 1.0, &(*this).v[mu][i], length, &CP.getV()[mu][i], length, &vReduced[i], length*prodOfSummations);
	        */
	        for(int i=0; i< (*this).componentDimensions[mu]; i++){
	           for(int j1=0; j1< prodOfSummations; j1++){
	              for(int j=0; j<k; j++){
	                 vReduced[i + (*this).componentDimensions[mu]*(j1 + prodOfSummations * j)]
	                       = (*this).v[mu][i+(*this).componentDimensions[mu]*j1]
	                       * CP.getV()[mu][i+(*this).componentDimensions[mu]*j];
	              }
	           }
	        }
	     }
	     if(mu == 0 || mu == (*this).d-1){
	        (*this).v[mu] = vReduced;
	        return(*this);
	     } else {
	        int copyLength = length * (*this).summations[mu-1];
	        int sliceLength = copyLength * k;
	        for(int j=0; j< CP.getSummation(0); j++){
	           for(int j2 = 0; j2 < (*this).summations[mu]; j2++){
	              //make diagonal block matrix (ij1j;j2j)
	              memcpy(&vFull[sliceLength *(j2 + j*(*this).summations[mu])+copyLength*j],&vReduced[copyLength * (j2 + j*(*this).summations[mu])],(*this).TSIZE*copyLength);
	           }
	        }
	        //print diagonal matrix
	        //printMatrix(length*(*this).summations[mu-1]*k, (*this).summations[mu]*k, &vFull[0]);
	        (*this).v[mu].resize(sizeFull);
	        memcpy(&(*this).v[mu][0], &vFull[0], sizeFull*(*this).TSIZE);
	        return(*this);
	     }

	  } // end partialHadamardProduct (const int mu, const CPTensorRepresentation<T> &CP)

  std::vector<T> partialAdd (const int mu, const MPSRepresentation<T> &summand){
     //tested for componentDimensions = n1 for all mu for (*this) and
     // componentDimensions = n2 for all mu for summand;
     int j1, j2;
     int k1, k2;
     int n = (*this).componentDimensions[mu];
     std::vector<T> vMu;

     if(mu == 0 || mu == (*this).d-1){
        // nodes with 1 edge
        if(mu == 0){
           j1 = (*this).summations[mu];
           k1 = summand.summations[mu];
        } else {
           j1 = (*this).summations[mu-1];
           k1 = summand.summations[mu-1];
        }
        vMu.resize(n*(j1+k1));

        memcpy(&vMu[0], &(*this).v[mu][0], (*this).TSIZE*n*j1);
        memcpy(&vMu[n*j1], &summand.v[mu][0], (*this).TSIZE*n*k1);
     } else {
        // mu = 1 ... d-2 (nodes with 2 edges)
        j1 = (*this).summations[mu-1];
        k1 = summand.summations[mu-1];
        j2 = (*this).summations[mu];
        k2 = summand.summations[mu];

        vMu.resize(n*(j1+k1)*(j2+k2));
        int sliceLength = n*(j1+k1);
        int copyLength = n*j1;
        int offset = 0;
        for(int i=0; i<j2; i++){
           //first block (n*(j1+k1); j2)
           //
           //            j2
           //       (          |
           //  j1  (   v(j1,j2)|
           //     (            |
           //     -------------|
           //     (            |
           //  k1  (     0     |
           //       (          |
           memcpy(&vMu[offset+i*sliceLength], &(*this).v[mu][i*copyLength], (*this).TSIZE*copyLength);
        }
        copyLength = n*k1;
        offset = n*j1+j2*sliceLength;
        for(int i=0; i<k2; i++){
           //second block (n*(j1+k1); k2)
           //
           //            j2         k2
           //       (          |          )
           //  j1  (   v(j1,j2)|     0     )
           //     (            |            )
           //     -------------|-------------
           //     (            |            )
           //  k1  (     0     | v(k1,k2)  )
           //       (          |          )
           memcpy(&vMu[offset+i*sliceLength], &summand.v[mu][i*copyLength], (*this).TSIZE*copyLength);
        }
     }
     return(vMu);
  } // end partialAdd (const int mu, const MPSRepresentation<T> &summand)




public:
	/* findet init nicht
	 MPSRepresentation(const std::vector<int> &summations, const std::vector<int> &componentDimensions){
	 int summationsCount = summations.size();
	 int d = componentDimensions.size();

	 std::vector< std::vector< int > > incidenceMatrix = setIncidenceMatrix(d);

	 init(summations, componentDimensions, incidenceMatrix);
	 }
	 */
	MPSRepresentation(const std::vector<T> &fullTensor, const std::vector<int> &componentDimensions,
                     const T eps = 0.0/*, bool display = false*/) {
	  MPSinit(fullTensor, componentDimensions, eps);
	}

	MPSRepresentation(const FullTensor<T> fullTensor, const T eps = 0.0){
	   MPSinit(fullTensor.getV(), fullTensor.getComponentDimensions(), eps);
	}

	MPSRepresentation(const TensorCalculus::CPTensorRepresentation<T> &CPTensor) {
		//triviale Umformung
		const int d = CPTensor.getD();
		const int k = CPTensor.getSummation(0);

		std::vector<std::vector<int> > incidenceMatrix = setIncidenceMatrix(d);
		std::vector<int> summations(d - 1);
		for (int i = 0; i < d - 1; i++) {
			summations[i] = k;
		}

		std::vector<int> componentDimensions(d);
		for (int i = 0; i < d; i++) {
			componentDimensions[i] = CPTensor.getComponentDimensions()[i];
		}

		std::vector<std::vector<T> > v(d);
		v[0] = CPTensor.getV()[0];
		v[d - 1] = CPTensor.getV()[d - 1];
		int n;
		for (int mu = 1; mu < d - 1; mu++) {
			n = componentDimensions[mu];
			v[mu].resize(k * n * k);
			for (int j1 = 0; j1 < k; j1++) {
				const int inc = 1;
				Blas<T>::copy(n, &CPTensor.getV()[mu][j1 * n], inc, &v[mu][n
						* (j1 + k * j1)], inc);
			}
			std::cout << v[mu] << std::endl;
		}
		init(summations, v, componentDimensions, incidenceMatrix);
	}//end MPSRepresentation(CPTensor)

	MPSRepresentation(const MPSRepresentation<T> &MPSTensor) {
	    init(MPSTensor.getSummations(), MPSTensor.getV(), MPSTensor.getComponentDimensions(), MPSTensor.getIncidenceMatrix());
	  }


	MPSRepresentation(std::vector<int> summations, std::vector< std::vector<T> > v, std::vector<int> componentDimensions) {
	   int d = componentDimensions.size();
	   init(summations, v, componentDimensions, setIncidenceMatrix(d));
 }//end MPSRepresentation(summations, v, componentDimensions)

	MPSRepresentation(std::vector<int> summations, std::vector<int> componentDimensions) {
     int d = componentDimensions.size();
     std::vector<std::vector<T> > v(d);
     int size;
     for(int mu=0; mu<d; mu++){
        size = componentDimensions[mu];
        if(mu == 0){
           size *= summations[mu];
        }else if (mu == d-1){
           size *= summations[d-2];
        }else {
           size *= summations[mu-1]*summations[mu];
        }
        v[mu].resize(size);
     }
	   init(summations, v, componentDimensions, setIncidenceMatrix(d));
	 }//end MPSRepresentation(summations, v, componentDimensions)

	MPSRepresentation(const char *filename){
	   read(filename);
	}

	MPSRepresentation(){

	}

	void write2disk (const char* MPSten){
	   std::ofstream fout(&MPSten[0]);
	   fout.setf(std::ios::scientific, std::ios::floatfield);
	   fout.precision(20);
	   fout<<"MPSTensor in d = "<<(*this).d<<std::endl;

	   for(int mu=0; mu<(*this).d; mu++){
       fout<<"n["<<mu<<"] = "<<(*this).componentDimensions[mu]<<std::endl;
	   }
	   for(int mu=0; mu<(*this).d-1; mu++){
	      fout<<"j["<<mu<<"] = "<<(*this).summations[mu]<<std::endl;
	   }
	   int n;
	   for(int mu=0; mu<(*this).d; mu++){
	      n = (*this).componentDimensions[mu];
	      if(mu == 0){
	         n *= (*this).summations[mu];
	      } else if (mu == (*this).d -1){
	         n *= (*this).summations[mu-1];
	      } else {
	         n *= (*this).summations[mu-1];
	         n *= (*this).summations[mu];
	      }
	      for(int i=0; i<n; i++){
	         fout<<(*this).v[mu][i]<<std::endl;
	      }
	   }
	}

	void read (const char* MPSten){
      std::ifstream from(&MPSten[0]);
      if(from.fail()){
         throw std::invalid_argument("There is no file to read.");
      }
      char t[64];
      int d;
      from.setf(std::ios::scientific, std::ios::floatfield);
      from >> t >> t >> t >> t;
      from >> d;
      std::vector<int> componentDimensions(d);
      for (int mu=0; mu<d; mu++){
       from >> t;
       from >> t;
       from >> componentDimensions[mu];
      }
      std::vector<int> summations(d-1);
      for (int mu=0; mu<d-1; mu++){
         from >> t;
         from >> t;
         from >> summations[mu];
      }
      std::vector<std::vector<T> > v(d);
      int n;
      for (int mu=0; mu<d; mu++){
         n = componentDimensions[mu];
         if(mu == 0){
            n *= summations[mu];
         } else if (mu == d -1){
            n *= summations[mu-1];
         } else {
            n *= summations[mu-1];
            n *= summations[mu];
         }
         v[mu].resize(n);
         for(int i=0; i<n; i++){
            from >> v[mu][i];
         }
      }
      init(summations, v, componentDimensions, setIncidenceMatrix(d));


      std::cout<<"read file: "<<MPSten<<std::endl;
      std::cout<<"d = "<<d<<std::endl;
      std::cout<<"componentDimensions = ";
      using namespace VectorOperators;
      std::cout<<componentDimensions<<std::endl;
      std::cout<<"summations = ";
      std::cout<<summations<<std::endl;
      //std::cout<<v <<std::endl;

	} //end read

	MPSRepresentation<T>& setHadamardProduct (const CPTensorRepresentation<T> &CP){
	   //tested for componentDimensions = n for all mu;
	   (*this).checkCompatibility(CP);

	   int d = (*this).d;
	   int k = CP.getSummation(0);
	   /*
	   std::vector<int> summations(d-1);
	   for(int mu=0; mu<d-1; mu++){
	      summations[mu] = (*this).summations[mu] * k;
	   }
	   std::vector<std::vector<T> > v(d);
	   v[0].resize((*this).componentDimensions[0]*summations[0]);
	   v[d-1].resize((*this).componentDimensions[d-1]*summations[d-2]);
	   for(int mu = 1; mu < d-1; mu++){
	      v[mu].resize((*this).componentDimensions[mu]*summations[mu-1]*summations[mu]);
	   }
	   */
	   for(int mu=0; mu<d; mu++){
	      (*this).partialHadamardProduct(mu, CP);
	      //memcpy(&v[mu][0],&had[0], v[mu].size()*(*this).TSIZE);
	   }
	   for(int mu=0; mu<d-1; mu++){
	       (*this).summations[mu] *= k;
	   }
	   return(*this);
	} //end setHadamardProduct (CPTensorRepresentation<T> CP)



	MPSRepresentation<T>& add (const MPSRepresentation<T> &summand){
	   //tested for componentDimensions = n1 for all mu for (*this) and
     // componentDimensions = n2 for all mu for summand;
	   (*this).checkCompatibility(summand);
	   int d = (*this).d;
	   int copyLength = (*this).componentDimensions[0]*summand.summations[0];
	   int offset     = (*this).componentDimensions[0]*(*this).summations[0];
	   (*this).v[0].resize(copyLength + offset);

	   memcpy(&(*this).v[0][offset], &summand.v[0][0], (*this).TSIZE * copyLength);

	   copyLength = (*this).componentDimensions[d-1]*summand.summations[d-2];
	   offset     = (*this).componentDimensions[d-1]*(*this).summations[d-2];
	   (*this).v[d-1].resize(copyLength + offset);

	   memcpy(&(*this).v[d-1][offset], &summand.v[d-1][0], (*this).TSIZE * copyLength);

	   for(int mu=1; mu<d-1; mu++){
	      /*
	      copyLength  = (*this).componentDimensions[mu];
	      copyLength *= summand.summations[mu-1] + (*this).summations[mu-1];
	      copyLength *= summand.summations[mu] + (*this).summations[mu];
	      (*this).v[mu].resize(copyLength);
	      memcpy(&(*this).v[mu][0], &partialAdd(mu,summand)[0],(*this).TSIZE * copyLength);
	      */
	      (*this).v[mu] = partialAdd(mu, summand);
	   }
	   for(int mu=0; mu<d-1; mu++){
        (*this).summations[mu] += summand.getSummation(mu);
	   }
	   return(*this);
	} // end add

	MPSRepresentation& edgeAdd (const int edge, const MPSRepresentation<T> &summand){
     //
	   // Addition for MPSRepresentations that differ in one edge and the corresponding nodes;
	   // sum_j'i^r' [sum_j1...ji-1 ji+1...jn v(j1) ... v'(ji-1,j'i) v'(j'i, ji+1) ... v(jn)]
	   //    + sum_j''i^r'' [sum_j1...ji-1 ji+1...jn v(j1) ... v''(ji-1,j''i) v''(j''i, ji+1) ... v(jn)]
	   //    = sum_ji^(r=r'+r'') [sum_j1...ji-1 ji+1...jn v(j1) ... v(ji-1,ji) v(ji, ji+1) ... v(jn)]
	   //
	   // mit v(ji-1,ji) = v'(ji-1, ji) für ji = 1...r'
	   //     v(ji-1,ji) = v''(ji-1, ji) für ji = r'+1...r''
	   //     v(ji,ji+1) = v'(ji, ji+1) für ji = 1...r'
     //     v(ji,ji+1) = v''(ji, ji+1) für ji = r'+1...r''
	   //
	   // i.e.: abcd + abef = ab(cd+ef)

	   (*this).checkCompatibility(summand);
     if(edge != 0 ){
        if(summand.summations[edge - 1] != (*this).summations[edge - 1]){
           throw std::invalid_argument("The tensor summations[edge - 1] is not equal.");
        }
     }
     if(edge != (*this).d - 2 ){
        if(summand.summations[edge + 1] != (*this).summations[edge + 1]){
           throw std::invalid_argument("The tensor summations[edge + 1] is not equal.");
        }
     }

     int size = (*this).componentDimensions[edge];

	   if(edge != 0){
	     size *= (*this).summations[edge - 1];
	   }
	   int copyLength = size * summand.summations[edge];
	   int offset     = size * (*this).summations[edge];
	   /*
	   std::cout<<"offset = "<<offset<<std::endl;
	   std::cout<<"oldsize = "<<(*this).v[edge].size()<<std::endl;
	   std::cout<<"copyL = "<<copyLength<<std::endl;
	   std::cout<<"sumSize = "<<summand.v[edge].size()<<std::endl;
	   */
	   size = copyLength + offset;
	   (*this).v[edge].resize(size);

	   memcpy(&(*this).v[edge][offset],&summand[edge][0], copyLength * (*this).TSIZE);

	   size = (*this).componentDimensions[edge+1];
     if(edge != (*this).d - 2){
        size *= (*this).summations[edge + 1];
        throw std::invalid_argument("not implemented");

     } else {
        offset     = size  * (*this).summations[edge];
        copyLength = size  * summand.summations[edge];
        size       = offset + copyLength;
        /*
        std::cout<<"offset = "<<offset<<std::endl;
        std::cout<<"oldsize = "<<(*this).v[edge+1].size()<<std::endl;
        std::cout<<"copyL = "<<copyLength<<std::endl;
        std::cout<<"sumSize = "<<summand.v[edge+1].size()<<std::endl;
        */
        (*this).v[edge+1].resize(size);

        memcpy(&(*this).v[edge+1][offset],&summand[edge+1][0], copyLength * (*this).TSIZE);
     }
     (*this).summations[edge]+=summand.summations[edge];

     return (*this);
	} // end edgeAdd

	MPSRepresentation& edgeAdd (const int edge, const MPSRepresentation<T> &summand, bool truncate, T eps = 0.0){
	   (*this).edgeAdd(edge, summand);
	   if(truncate = true){
	      //todo
	   }
	   return(*this);
	}


  MPSRepresentation& truncateMPS(const T eps, const char epsWeight = 'e', const bool display = false) {

     //truncation for a MPSTensor as discribed by Osoledets;
     //nicht complet implementiert
     if(display){std::cout<<"truncateMPS"<<std::endl;}

     time_t start, end;
     time(&start);

     std::vector< std::vector <T> > vTemp;

     mpsqlf(vTemp);

     time(&end);
     if(display){std::cout<<"LQ-factorisation took "<<difftime(end,start)<<"s"<<std::endl;}
     time(&start);

     mpssvd(eps, epsWeight, display);

     time(&end);
     if(display){std::cout<<"SVD took "<<difftime(end,start)<<"s"<<std::endl;}
     return(*this);
  } // end truncateMPS

  T tensorScalarProduct(MPSRepresentation<T> mps2){
     //throw std::invalid_argument("not implemented");
     T result;
     int d = (*this).d;
     (*this).checkCompatibility(mps2);
     int wSize = *( std::max_element ( (*this).summations.begin(), (*this).summations.end() ) );
     wSize    *= *( std::max_element ( mps2.summations.begin(), mps2.summations.end() ) );
     int vSize = wSize * *( std::max_element ( (*this).componentDimensions.begin(), (*this).componentDimensions.end() ) );
     std::vector<double> v (vSize);
     std::vector<double> w (wSize);

     Blas<T>::gemm('t', 'n', (*this).summations[0], mps2.summations[0], (*this).componentDimensions[0],
           1.0, &(*this).v[0][0], (*this).componentDimensions[0], &mps2.v[0][0], mps2.componentDimensions[0],
           0.0, &w[0], (*this).summations[0]);
     for(int j=0; j<(*this).summations[1]; j++){
        Blas<T>::gemm('n', 'n', (*this).componentDimensions[1], mps2.summations[0], (*this).summations[0],
              1.0, &(*this).v[1][j * (*this).componentDimensions[1] * (*this).summations[0]],
              (*this).componentDimensions[1], &w[0], (*this).summations[0], 0.0,
              &v[j * (*this).componentDimensions[1] * mps2.summations[0]], (*this).componentDimensions[1]);
     }
     for(int mu=1; mu<d-1; mu++){
        int k = (*this).componentDimensions[mu]*mps2.summations[mu-1];
        Blas<T>::gemm('t', 'n', (*this).summations[mu], mps2.summations[mu], k, 1.0, &v[0],
              k , &mps2.v[mu][0], k, 0.0, &w[0], (*this).summations[mu]);
        int summation;
        if(mu != d-2){
           summation = (*this).summations[mu+1];
        } else {
           summation = 1;
        }
        for(int j=0; j < summation; j++){
           Blas<T>::gemm('n', 'n', (*this).componentDimensions[mu+1], mps2.summations[mu], (*this).summations[mu],
                         1.0, &(*this).v[mu+1][j * (*this).componentDimensions[mu+1] * (*this).summations[mu]],
                         (*this).componentDimensions[mu+1], &w[0], (*this).summations[mu], 0.0,
                         &v[j * (*this).componentDimensions[mu+1] * mps2.summations[mu]], (*this).componentDimensions[mu+1]);
        }
     }
     result = Blas<T>::dot((*this).componentDimensions[d-1]*mps2.summations[d-2], &v[0], 1,
                          &mps2.v[d-1][0], 1);
     return result;
  }

};
}

#endif /* __MPSREPRESENTATION_HPP */
