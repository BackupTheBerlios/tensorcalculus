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

#ifndef MPSDECOMPOSITION_HPP_
#define MPSDECOMPOSITION_HPP_

#include "Representation/MPSRepresentation.hpp"
#include <vector>
#include "BlasInterface.hpp"
#include "Vector/VectorOperators.hpp"


namespace TensorCalculus{
   template <typename T>
   class MPSDecomposition{
   private:
      int d;
      std::vector<int> componentDimensions;
      std::vector<int> summations;
      std::vector<int> n;

      //output for MPSRepresentation
      std::vector< std::vector<T> > v;

      // on input: Vector (i_1,...,i_d)
      // Vector with multidirections (j,i_n,...,i_d),
      std::vector<T> tempV;
      std::vector<T> VT; // saves VT if JOBVT = 'S'
      std::vector<T> U;  // saves U if JOBU = 'S'
      std::vector<T> svdWork;
      int            svdLwork;
      std::vector<T> singularValues;

      //for truncation
      T epsNull;
      std::vector<T> epsMu;

      void init(const std::vector<T> &fullTensor,const std::vector<int> &componentDimensions){
         (*this).componentDimensions = componentDimensions;
         tempV = fullTensor;
         d = componentDimensions.size();
         v.resize(d);
         summations.resize(d-1);
         n.resize(d);
         n[d-1] = 1;
         for(int mu=d-1; mu>0; mu--){
            n[mu-1] = n[mu] * (*this).componentDimensions[mu];
         }
         svdWork.resize(1);
         svdLwork = 1;
         epsNull = 1e-12;
         //todo: größe richtig setzen
      };

      void vidalSVD(int mu){
         char jobvt = 'S';
         char jobu  = 'S';

         int m = componentDimensions[mu];
         if(mu != 0){
            m *= summations[mu-1];
         }
         int min = std::min(m, n[mu]);

         singularValues.resize(min);

         int Usize = min*m;
         if(U.size() < Usize){
            U.resize(Usize);
         }

         int VTsize = min * n[mu];
         if(VT.size() < VTsize){
            VT.resize(VTsize);
         }
         //std::cout<<"set work"<<std::endl;
         int setWork = -1;

        Lapack<T>::gesvd(jobu, jobvt, m, n[mu],&tempV[0] , m, &singularValues[0], &U[0], m, &VT[0], min, &svdWork[0], setWork);
        if (svdWork[0] > svdLwork) {
          svdLwork = svdWork[0];
          svdWork.resize(svdLwork);
        }
        Lapack<T>::gesvd(jobu, jobvt, m, n[mu],&tempV[0] , m, &singularValues[0], &U[0], m, &VT[0], min, &svdWork[0], svdLwork);
        //std::cout<<"end vidalSVD"<<std::endl;
      } //end vidalSVD

      void truncateVT (int mu, int rFull, int rTrunc){
         for (int i = n[mu]-1; i >= 0; i--){
             VT.erase(VT.begin() + rTrunc + i * rFull, VT.begin() + rFull
                 + i * rFull);
         }
      }

      void setVidalV (int mu, int rTrunc){
         int rOld = 1;
         if(mu != 0){
            rOld = summations[mu-1];
         }
         v[mu].resize(rOld * rTrunc * componentDimensions[mu]);
         Blas<T>::copy(componentDimensions[mu]*rTrunc*rOld, &U[0], 1, &v[mu][0], 1);
         //std::cout<<"set tempV"<<std::endl;
         if(mu != d-1){
          for(int j=0; j<n[mu+1]; j++){
            transpose(rTrunc, componentDimensions[mu+1],
                  &VT[j*rTrunc*componentDimensions[mu+1]],
                  &tempV[j*rTrunc*componentDimensions[mu+1]]);
            for(int i=0; i<rTrunc; i++){
               Blas<T>::scal(componentDimensions[mu+1],singularValues[i],
                     &tempV[j*rTrunc*componentDimensions[mu+1]+i*componentDimensions[mu+1]], 1);
            }
          }
         }else{
           v[d-1].resize(rTrunc * componentDimensions[d-1]);
           Blas<T>::copy(rTrunc * componentDimensions[d-1], &tempV[0], 1, &v[d-1][0], 1);
         }
      } //end setVidalV



   public:
      MPSDecomposition(const std::vector<T> &fullTensor,const std::vector<int> &componentDimensions){
         init(fullTensor, componentDimensions);
      }

      std::vector< std::vector <T> > getV (){return(v);}
      std::vector< int >             getSummations (){return(summations);}
      std::vector< int >             getComponentDimensions (){return(componentDimensions);}

      void computeVidalDecomposition (){
         computeVidalDecomposition (epsNull);
      }
      void computeVidalDecomposition (const T eps){
         setEps(d, epsMu, eps, 'e', true);
         int rTrunc;
         int rFull;
         for (int mu=0; mu<d; mu++){
            //std::cout<<"mu = "<<mu<<std::endl;
            if(mu != d-1){
               //std::cout<<"computeSVD"<<std::endl;
               vidalSVD(mu);
               rFull = singularValues.size();
               rTrunc = rFull;
               //std::cout<<"evalaute eps"<<std::cout;
               evaluateEps(epsMu, singularValues, mu, rTrunc);
               if (rTrunc < 1){
                  rTrunc = 1;
               }
               summations[mu]=rTrunc;
               if(rTrunc < rFull){
                  truncateVT(mu, rFull, rTrunc);
               }
            }
            //std::cout<<"setV"<<std::endl;
            setVidalV(mu, rTrunc);
         }
      }//end computeVidalDecompostition

      void computeVidalDecomposition (std::vector<int> maxSummations){
         if(maxSummations.size() != d-1){
            std::cout<<"maxSummations has wrong size"<<std::endl;
         }

         setEps(d, epsMu, epsNull);
         int rTrunc;
         int rFull;
         for (int mu=0; mu<d; mu++){
            //std::cout<<"mu = "<<mu<<std::endl;
            if(mu != d-1){
               //std::cout<<"computeSVD"<<std::endl;
               vidalSVD(mu);
               rFull = singularValues.size();
               rTrunc = rFull;
               //std::cout<<"evaluateEps"<<std::endl;
               evaluateEps(epsMu, singularValues, mu, rTrunc);

               if(maxSummations[mu] < rTrunc){
                  rTrunc = maxSummations[mu];
               }
               if (rTrunc < 1){
                  rTrunc = 1;
               }
               summations[mu] = rTrunc;

               if(rTrunc < rFull){
                  truncateVT(mu, rFull, rTrunc);
               }
            }
            //std::cout<<"setV"<<std::endl;
            setVidalV(mu, rTrunc);
         }
      }//end computeVidalDecompostition

   }; //end class MPSDecomposition
   template <typename T>
   void setEps(const int d, std::vector<T> &epsMu, const T eps = 1e-15, const char weight = 'e', bool display = false){
      //set EpsMu[0]...Eps[d-2] = epsMu, sqrt(sum(i=rTrunc to rFull-1)  S[i]*S[i]) < epsMu
      //set EpsMu[d-1] = epsMu*epsMu, overwrite after evaluation with squared border for next iteration

      // e...all epsMu are equal
      // m...bigger epsMu in the middle, implemented for d=4 : 1 k 1 in epsMu, k>4
      // f...whole epsilon for the first
      /*
      if (epsMu.size() < d){
         epsMu.resize(d);
      }
      */
      if(weight == 'e'){
         if(display){std::cout<<"use EpsE"<<std::endl;}
         setEpsE(d, eps, epsMu);
      }
      if(weight == 'm'){
         if(display){std::cout<<"use EpsM"<<std::endl;}
         setEpsM(d, eps, epsMu);
      }
      if(weight == 'f'){
        if(display){std::cout<<"use EpsM"<<std::endl;}
        setEpsF(d, eps, epsMu);
      }
      //std::cout<<"set Eps finish"<<std::endl;
   }

   template <typename T>
   void setEpsE(const int d, const T eps, std::vector <T> &epsMu){
      const T eps1Mu = eps / (T)(d-1);
      const T eps2Mu = eps1Mu * eps1Mu;
      if(epsMu.size() < d){
         epsMu.resize(d);
      }
      epsMu[d-1] = eps2Mu;
      for(unsigned int mu=0; mu<d-1; mu++){
         epsMu[mu] = eps1Mu;
      }
   }

   template <typename T>
   void setEpsM(const int d, const T eps, std::vector <T> &epsMu){
        //for d=4
        if(d == 4){
           double k = 1000.0; // epsMu = (1, k, 1) in eps1Mu
           const T eps1Mu = eps / (T)(k+2);
           const T eps2Mu = eps1Mu * eps1Mu;
           if(epsMu.size() < d){
              epsMu.resize(d);
           }
           epsMu[0] = eps1Mu;
           epsMu[1] = eps1Mu *k;
           epsMu[2] = eps1Mu;
           epsMu[d-1] = eps2Mu;
        }else{
           std::cout<<"setEps not implemented"<<std::endl;
        }
     }

   template <typename T>
   void setEpsF(const int d, const T eps, std::vector <T> &epsMu){
        const T eps1Mu = eps;
        const T eps2Mu = eps1Mu * eps1Mu;
        if(epsMu.size() < d){
           epsMu.resize(d);
        }
        epsMu[d-1] = eps2Mu;
        epsMu[0] = eps1Mu;
     }

   template <typename T>
   void evaluateEps(std::vector<T> &epsMu, const std::vector<T> &singularValues, const int mu, int &rTrunc){
      //input: rTrunc == min (m,n)
      int d = epsMu.size();
      T s2 = singularValues[rTrunc-1]*singularValues[rTrunc-1];
      T s2sum = s2;
      T s2sumOld = 0;
      int i = 0;
      while(s2sum < epsMu[d-1]){
         i++;
         s2sumOld = s2sum;
         s2 = singularValues[rTrunc-i-1]*singularValues[rTrunc-i-1];
         s2sum += s2;
      }
      rTrunc -= i;
      if(mu != d-2){
         epsMu[mu+1] += epsMu[mu] - sqrt(s2sumOld);
         epsMu[d-1] = epsMu[mu+1]*epsMu[mu+1];
      }
   }


} //end namespace TensorCalculus

#endif /* MPSDECOMPOSITION_HPP_ */
