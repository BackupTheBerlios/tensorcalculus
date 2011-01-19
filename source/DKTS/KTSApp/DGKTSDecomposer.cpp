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

// DGKTSDecomposer.cpp: Implementierung der Klasse DGKTSDecomposer.
//
//////////////////////////////////////////////////////////////////////

#include "DGKTSDecomposer.hpp"



DGKTSDecomposer::DGKTSDecomposer()
 {
 }


DGKTSDecomposer::~DGKTSDecomposer()
 {

 }

DGKTSDecomposer& DGKTSDecomposer::setDefaultParameter(const LongReal& normA)
 {
    const LongInt d = DGKTSDDataBase::d();
    const LongInt r = DGKTSDDataBase::k();
 
    const LongReal nA = normA*normA;
    const LongReal nD = pow(nA*nA, 1.0/((LongReal)d))*(LongReal)(d*r);

    (*this)
     .setMu(1.0)
     .setGamma(0.9)
     .setDelta(1.0e-1)
     .setKappa(0.0)
     .setKappa1(1.0/nA)
     .setKappa2(1.0/nD)
     .setKappa3(1.0/nA)
     .setLambda1(1.0e-3)
     .setLambda2(0.0)
     .setEpsilonCR(1.0e-10)
     .setMaxStepsCR(350)
     .setMaxStepsCRL1(30)
     .setMaxStepsCRL2(10)
    ;
        
   return (*this);
 }


DKTSDIterationInfo DGKTSDecomposer::decompose(DKTS& a, DKTS& xi, const LongReal& normA, DKTSDDataBlock& dataBlock)
 {
    DKTSDIterationInfo infoEntry;
            
    if(resize(a, xi)==true)
     {
       if(20<=a.k())
        {
           (*this)
            .setUseNewtonMethode(false)
           ;
        }
       else
        {
           (*this)
            .setUseNewtonMethode(true)
           ;
        }

        xi *= 1.0/normA;
        a  *= 1.0/normA;        
/*
        if(attr_printCout)
         {
            cout << "|A| = " << normA << '\t';
         }
*/
        (*this)
         .setDefaultParameter(1.0)
        ;
                
        infoEntry = startIteration(a, xi, 1.0, dataBlock);
        
        xi *= normA;
        a  *= normA;
     }
    else
     {
        throw SimpleException(IString("Warning In DGKTSD::decompose(const DKTS& a, DKTS& x, const LongReal& normA, DKTSDDataBlock& dataBlock), Error in resize !!!"));
     }

   return infoEntry;
 }

 
DKTSDIterationInfo DGKTSDecomposer::decompose(DKTS& a, DKTS& xi, DKTSDDataBlock& dataBlock)
 {
    DKTSDIterationInfo infoEntry;
            
    if(resize(a, xi)==true)
     {
       if(20<=a.k())
        {
           (*this)
            .setUseNewtonMethode(false)
           ;
        }
       else
        {
           (*this)
            .setUseNewtonMethode(true)
           ;
        }

        (*this)
         .setDefaultParameter(1.0)
        ;
                
        infoEntry = startIteration(a, xi, 1.0, dataBlock);
        
     }
    else
     {
        throw SimpleException(IString("Warning In DGKTSD::decompose(const DKTS& a, DKTS& x, DKTSDDataBlock& dataBlock), Error in resize !!!"));
     }

   return infoEntry;
 } 
