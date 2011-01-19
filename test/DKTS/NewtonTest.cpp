/*
 * Copyright (C) Mike Espig
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

#include "DKTSTruncation2Eps.hpp"

#include "SimpleException.hpp"

#include <iostream>
#include <fstream>

#include <stdio.h>

int truncate2Eps(LongReal& eps, LongReal& epsN, LongReal& preC, LongInt& kMin)
 {
    const IString fileName   ("tensor_input.ten");
    const IString coFileName ("coefficientSystem.ten");
    const IString ortFileName("orthogonalBasis.ten");

    Timer time;

    IString date(time.date());

    cout << "Tensor Decompositions       : " << endl;
    cout << "Date                        : " << date;

    cout << endl;
    cout << "Truncate2Eps " << endl;

    cout << "kMin         = ";
    cin  >> kMin;

    cout << "eps          = ";
    cin  >> eps;

    cout << "precision    = ";
    cin  >> preC;

    cout << endl;


    eps = MAX(1.0e-7, eps);

    cout << "Reading Initial-Tensor from : " << fileName << endl << endl;


    DKTSTruncation2Eps truncator(fileName, kMin);
    //DKTSTruncation2Eps truncator(fileName, coFileName, ortFileName, kMin);

    truncator.truncateF2(eps, epsN, preC);

/*
    IString name("DKTSTruncation");

    truncator.attr_truncationLog.plotAll("truncationLog", name);
    truncator.attr_truncationLogR1.plotAll("truncationLogR1", name);
*/

   return 0;
 }


int computeCoefficientTensor()
 {
    const IString fileName   ("tensor_input.ten");
    const IString coFileName ("coefficientSystem.ten");
    const IString ortFileName("orthogonalBasis.ten");

    DKTS A, Z, a;

    A.readDataFrom(fileName);

    Z.setOrthogonal2(A);
    a.setCoefficientsSystemOf(A, Z);

    a.writeDataTo(coFileName);
    Z.writeDataTo(ortFileName);

   return 0;
 }


int main(int argc, char* argv[]) {
	 LongInt kMin = 1;
	 LongInt kMax = 9;
	 LongInt r    = 1;

	 LongReal eps   = 1.0e-3;
    LongReal preC  = 1.0e-3;
    LongReal epsN  = 1.0e-4;

    if(argc==2)
     {
        IString temp(argv[1]);
        eps   = MAX(temp.asDouble(), 1.0e-12);
        kMax  = kMin;
     }
    else if(argc==3)
     {
        IString temp1(argv[1]);
        eps   = MAX(temp1.asDouble(), 1.0e-12);
        kMax  = kMin;

        IString temp2(argv[2]);
        preC   = MAX(temp2.asDouble(), 1.0e-10);
     }
    else if(argc==4)
     {
        IString temp1(argv[1]);
        eps   = MAX(temp1.asDouble(), 1.0e-12);
        kMax  = kMin;

        IString temp2(argv[2]);
        preC   = MAX(temp2.asDouble(), 1.0e-10);

        IString temp3(argv[3]);
        epsN   = MAX(temp3.asDouble(), 1.0e-12);
     }


    try
     {
        //computeCoefficientTensor();
        truncate2Eps(eps, epsN, preC, r);
     }catch(SimpleException e)
      {
         cout << e.text() << endl;
      }

}

