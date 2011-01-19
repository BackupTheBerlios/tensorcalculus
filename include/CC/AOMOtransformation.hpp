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

#ifndef AOMOTRANSFORMATION_HPP_
#define AOMOTRANSFORMATION_HPP_

#include <vector>
#include "Representation/MPSRepresentation.hpp"
#include "Representation/TensorRepresentation.hpp"

namespace TensorCalculus{


void AO2MO (const TensorRepresentation<double> &in, TensorRepresentation<double> &virt, TensorRepresentation<double> &occ);
void makeV (const std::vector<char> occState, TensorRepresentation<double> &V_mnop,
      const TensorRepresentation<double> &virt, const TensorRepresentation<double> &occ);

MPSRepresentation<double> makeV (const std::vector<char> &occState, const MPSRepresentation<double> &AO);

//for testing
class AOMOtransformation {

   protected:
      int   aoCount;
      int  occCount;
      int virtCount;
      std::vector<double>  occAOMO;
      std::vector<double> virtAOMO;
      std::vector<double> fullTransposeAOMO;
      std::vector<double> fullTensorAO;
      std::vector<double> fullTensorMO;
      MPSRepresentation<double> AO;
      std::vector< std::vector <double> > MOv;

      std::vector< int > prodOfRank;

      void init(const int virt, const int occ);

   private:

   public:
      AOMOtransformation (const int virtCount, const int occCount){
         init(virtCount, occCount);
      };

      //~AOMOtransformation ();

      //void AO2vvooMO ();
      void AO2MO ();
      void fullTensorAO2fullTensorMO (); //tested --> M(a, mu) M(b, nu) M(c, lam) M(d, sig)
      MPSRepresentation<double> getVabij ();
      MPSRepresentation<double> AO2Vabij ();
      MPSRepresentation<double> getFullMO ();

      std::vector<double> getOccAOMO (){return(occAOMO);};
      std::vector<double> getVirtAOMO (){return(virtAOMO);};
      std::vector<double> getFullTransposeAOMO (){return(fullTransposeAOMO);};
      std::vector<double> getFullTensorMO (){return(fullTensorMO);};

      /////////////////////////////////////////////////////////////////////////////////////////


}; /* end class AOMOtransformation */

/*
void partialAOMO (const int aoCount, const int moCount, const int productOfRanks,
                  const std::vector<double> &AOMO, const std::vector<double> &AO,
                  std::vector<double> &MO);

void AOMOtransformation (const TensorCalculus::MPSRepresentation<double> &AO, std::vector<double> &MO);
*/

} /*end namespace TensorCalculus */
#endif /* AOMOTRANSFORMATION_HPP_ */
