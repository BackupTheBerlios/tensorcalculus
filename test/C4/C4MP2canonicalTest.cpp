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

#include "CC/MP2.hpp"

using namespace TensorCalculus;

int main (){
   std::ofstream fout ("compute/mp2_canonical_log");
   fout.precision(12);
   double E2;
   E2 = mp2_ft ("v_abij.ften", "t_abij.ten");
   fout<<"MP2, conventionel   "<< E2<<std::endl;
   fout<<"### chemical ###"<<std::endl;
   E2 = mp2_mps_canonical(1.0, "compute/mp2_canonical_chem_1e+0", true);
   fout<<"1.00E+00   "<<E2<<std::endl;
   E2 = mp2_mps_canonical(0.8, "compute/mp2_canonical_chem_8e-1", true);
   fout<<"8.00E-01   "<<E2<<std::endl;
   E2 = mp2_mps_canonical(0.6, "compute/mp2_canonical_chem_6e-1", true);
   fout<<"6.00E-01   "<<E2<<std::endl;
   E2 = mp2_mps_canonical(0.4, "compute/mp2_canonical_chem_4e-1", true);
   fout<<"4.00E-01   "<<E2<<std::endl;
   E2 = mp2_mps_canonical(0.2, "compute/mp2_canonical_chem_2e-1", true);
   fout<<"2.00E-01   "<<E2<<std::endl;
   E2 = mp2_mps_canonical(0.1, "compute/mp2_canonical_chem_1e-1", true);
   fout<<"1.00E-01   "<<E2<<std::endl;
   E2 = mp2_mps_canonical(0.08, "compute/mp2_canonical_chem_8e-2", true);
   fout<<"8.00E-02   "<<E2<<std::endl;
   E2 = mp2_mps_canonical(0.06, "compute/mp2_canonical_chem_6e-2", true);
   fout<<"6.00E-02   "<<E2<<std::endl;
   E2 = mp2_mps_canonical(0.04, "compute/mp2_canonical_chem_4e-2", true);
   fout<<"4.00E-02   "<<E2<<std::endl;
   E2 = mp2_mps_canonical(0.02, "compute/mp2_canonical_chem_2e-2", true);
   fout<<"2.00E-02   "<<E2<<std::endl;
   E2 = mp2_mps_canonical(0.01, "compute/mp2_canonical_chem_1e-2", true);
   fout<<"1.00E-02   "<<E2<<std::endl;
   E2 = mp2_mps_canonical(0.008, "compute/mp2_canonical_chem_8e-3", true);
   fout<<"8.00E-03   "<<E2<<std::endl;
   E2 = mp2_mps_canonical(0.006, "compute/mp2_canonical_chem_6e-3", true);
   fout<<"6.00E-03   "<<E2<<std::endl;
   E2 = mp2_mps_canonical(0.004, "compute/mp2_canonical_chem_4e-3", true);
   fout<<"4.00E-03   "<<E2<<std::endl;
   E2 = mp2_mps_canonical(0.002, "compute/mp2_canonical_chem_2e-3", true);
   fout<<"2.00E-03   "<<E2<<std::endl;
   E2 = mp2_mps_canonical(0.001, "compute/mp2_canonical_chem_1e-3", true);
   fout<<"1.00E-03   "<<E2<<std::endl;
   E2 = mp2_mps_canonical(0.0001, "compute/mp2_canonical_chem_1e-4", true);
   fout<<"1.00E-04   "<<E2<<std::endl;
   E2 = mp2_mps_canonical(0.00001, "compute/mp2_canonical_chem_1e-5", true);
   fout<<"1.00E-05   "<<E2<<std::endl;

   return 0;
}

