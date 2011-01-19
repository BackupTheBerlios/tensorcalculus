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

#include "Representation/MPSRepresentation.hpp"
#include <vector>
#include "CC/MP2.hpp"
#include "CC/C4Interface.hpp"
#include "Tensor/FullTensor.hpp"

using namespace TensorCalculus;

int test1 () {
   std::ifstream fin ("compute/MP2");
   if (!fin.fail()){
      //throw std::invalid_argument("MP2 file already exist");
   }
   std::ofstream fout ("compute/MP2");
   fout.precision(15);
   double E2;
   E2 = mp2_ft ("v_abij.ften","t_abij.ten");
   std::cout<<" from canonical quantities, full MP2 ... "<< E2 <<std::endl;
   fout<<" from canonical quantities, full MP2 ... "<< E2 <<std::endl;

   /*
   FullTensor<double> Vft ("v_abij.ften");
   FullTensor<double> Tft ("t_abij.ten");
   MPSRepresentation<double> Vmps (Vft, 1e-12);
   MPSRepresentation<double> Tmps (Tft, 1e-12);
   */


   MPSRepresentation<double> Vmps ("compute/AOchem2Vaibj.mps.1e-12");
   MPSRepresentation<double> Tmps ("compute/AOchem2V2Taibj.mps.1e-12");

   E2 = mp2_mps(Vmps, Tmps, true);

   //E2 = mp2_mps ("compute/AO2Vabij.mps.1e-12", "compute/AO2V2Tabij.mps.1e-12");
   return 0;
}

int main (){
   //todo: E2 schritt für chemical verdrahten

   //mp2_mps(1e-12, "compute/mp2_chem_1e-12", true);
   mp2_mps(1e-2, "compute/mp2_chem_1e-2", true);
   mp2_mps(1e-4, "compute/mp2_chem_1e-4", true);
   mp2_mps(1e-6, "compute/mp2_chem_1e-6", true);

   mp2_mps(1e-2, "compute/mp2_phys_1e-2", false);
   mp2_mps(1e-4, "compute/mp2_phys_1e-4", false);
   mp2_mps(1e-6, "compute/mp2_phys_1e-6", false);
   return 0;
}

