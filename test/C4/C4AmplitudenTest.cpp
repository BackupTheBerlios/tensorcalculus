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

#include <vector>
#include "CC/SpinAdapt.hpp"
#include "CC/C4Interface.hpp"
#include "Tensor/FullTensor.hpp"
#include "Representation/MPSRepresentation.hpp"
#include "CC/AOMOtransformation.hpp"
#include "CC/Amplitudes.hpp"

using namespace TensorCalculus;
using namespace VectorOperators;

void test ();
void test2 ();

void makeAmplitudesPhys(const char *MOin, const char *Tout, const double eps); // uses EpsInv from DKTS (udo)
//void makeAmplitudesChem(const char *MOin, const char *Tout, const double eps);

void makeAmplitudes (const char *MOin, const char *Tout, std::vector<char> &occState, const double eps);

int main () {
   //test();
   //test2();

   std::vector<char> occState(4);
   occState[0] = 'v';
   occState[1] = 'v';
   occState[2] = 'o';
   occState[3] = 'o';

   makeAmplitudes("compute/AOphys2Vabij.mps.1e-12", "compute/AOphys2V2Tabij.mps.1e-12", occState, 1e-12);
   makeAmplitudes("compute/AOphys2Vabij.mps.1e-6", "compute/AOphys2V2Tabij.mps.1e-6", occState, 1e-6);
   makeAmplitudes("compute/AOphys2Vabij.mps.1e-4", "compute/AOphys2V2Tabij.mps.1e-4", occState, 1e-4);
   makeAmplitudes("compute/AOphys2Vabij.mps.1e-2", "compute/AOphys2V2Tabij.mps.1e-2", occState, 1e-2);

   occState[0] = 'v';
   occState[1] = 'o';
   occState[2] = 'v';
   occState[3] = 'o';

   makeAmplitudes("compute/AOchem2Vaibj.mps.1e-12", "compute/AOchem2V2Taibj.mps.1e-12", occState, 1e-12);
   makeAmplitudes("compute/AOchem2Vaibj.mps.1e-6", "compute/AOchem2V2Taibj.mps.1e-6", occState, 1e-6);
   makeAmplitudes("compute/AOchem2Vaibj.mps.1e-4", "compute/AOchem2V2Taibj.mps.1e-4", occState, 1e-4);
   makeAmplitudes("compute/AOchem2Vaibj.mps.1e-2", "compute/AOchem2V2Taibj.mps.1e-2", occState, 1e-2);
}

void test(){
   std::vector<double> fullSpinAdAO;

   //readAOSpinAdapted(fullSpinAdAO);

   readAO(fullSpinAdAO);
   //Blas<double>::scal(fullSpinAdAO.size(), 4.0, &fullSpinAdAO[0],1);

   FullTensor<double> ftSpinAdAO(makeFullTensorAO(fullSpinAdAO));
   std::cout<<ftSpinAdAO.getComponentDimensions()<<std::endl;
   MPSRepresentation<double> mpsSpinAdAO (ftSpinAdAO, 1e-12);

   std::vector<char> occState(4);
   occState [0] = 'v';
   occState [1] = 'v';
   occState [2] = 'o';
   occState [3] = 'o';

   //MPSRepresentation<double>  t_abij(makeV(occState, mpsSpinAdAO));
   FullTensor<double> ft_v ("v_abij.ften");
   MPSRepresentation<double> t_abij(ft_v, 1e-12);
   std::vector<double> t_vor(t_abij.evaluate());
   V2T(t_abij);
   t_abij.truncateMPS(1e-12);
   std::vector<double> t_nach(t_abij.evaluate());
   std::cout<<"|| t_vor - t_nach|| "<<l2_norm(t_vor-t_nach)<<std::endl;

   FullTensor<double> read_t("t_abij.ten");

   double dist = l2_norm(t_abij.evaluate()-read_t.getV());
   double rel_dist = dist / l2_norm(read_t.getV());

   std::cout<<"dist = "<<dist<<"   rel dist = "<<rel_dist<<std::endl;
   FullTensor<double> t_test (t_abij.getComponentDimensions(), t_abij.evaluate());
   t_test.write2disk("t_test");
}

void makeAmplitudesPhys(const char *MOin, const char *Tout, const double eps){
   std::cout<<"make Amplitudes Physical ...";
   std::ifstream fin (&MOin[0]);
   if(fin.fail()){
      throw std::invalid_argument("no input file");
   }
   std::ifstream fin2 (&Tout[0]);
   if(!fin2.fail()){
      std::cout<<"output file already exist"<<std::endl;
      return;
   }

   MPSRepresentation<double> t_abij(&MOin[0]);
   V2T(t_abij);
   t_abij.truncateMPS(eps);
   t_abij.write2disk(&Tout[0]);

   std::cout<<"make Amplitudes Physical ... finshed"<<std::endl;
}

void makeAmplitudes (const char *MOin, const char *Tout, std::vector<char> &occState, const double eps){
   std::cout<<"make Amplitudes ...";
   std::ifstream fin (&MOin[0]);
   if(fin.fail()){
      throw std::invalid_argument("no input file");
   }
   std::ifstream fin2 (&Tout[0]);
   if(!fin2.fail()){
      std::cout<<"output file already exist"<<std::endl;
      return;
   }

   MPSRepresentation<double> t_abij(&MOin[0]);
   V2T(t_abij, occState);
   t_abij.truncateMPS(eps);
   t_abij.write2disk(&Tout[0]);

   std::cout<<"make Amplitudes ... finshed"<<std::endl;
}

void test2 (){
   std::vector<double> epsFromDKTS (readEpsilonInvDKTS().evaluate());
   //std::cout<<epsFromDKTS<<std::endl;
   std::vector<char> occState(4);
   occState[0] = 'v';
   occState[1] = 'v';
   occState[2] = 'o';
   occState[3] = 'o';
   std::vector<double> epsFromEPS (makeEpsilonInv(occState).evaluate());
   //std::cout<<"DKTS"<<epsFromDKTS<<std::endl;
   //std::cout<<"EPSILON"<<epsFromEPS<<std::endl;
   std::cout<<"DKTS: rel dist : "<<l2_norm(epsFromEPS - epsFromDKTS) / l2_norm(epsFromDKTS)<<std::endl;

   std::vector<double> f;
   readOrbitalEnergies(f);

   int occ, virt;
   readOccN(occ, occ, virt, virt);
   std::vector<double> fullEpsInv (occ*occ*virt*virt);
   for(int k=0; k<occ; k++){
      for(int l=0; l<occ; l++){
         for(int m=0; m<virt; m++){
            for(int n=0; n<virt; n++){
               fullEpsInv[n + virt*(m + virt * (l + occ * k))] = 1 / (-1.0 * f[k] - f[l] - f[occ+m] - f[occ+n]);
            }
         }
      }
   }
   std::cout<<"full: rel dist : "<<l2_norm(epsFromEPS - fullEpsInv) / l2_norm(fullEpsInv)<<std::endl;
   std::cout<<"DKTS -full: rel dist : "<<l2_norm(epsFromDKTS - fullEpsInv) / l2_norm(fullEpsInv)<<std::endl;
}
