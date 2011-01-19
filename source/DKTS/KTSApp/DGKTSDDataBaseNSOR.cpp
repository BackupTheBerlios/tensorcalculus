/*
 * Copyright (C) 2011 Mike Espig
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

// DGKTSDDataBaseNSOR.cpp: Implementierung der Klasse DGKTSDDataBaseNSOR.
//
//////////////////////////////////////////////////////////////////////

#include "DGKTSDDataBaseNSOR.hpp"


DGKTSDDataBaseNSOR::DGKTSDDataBaseNSOR()
 {
    (*this)
     .setDefaultParameter(1.0)
     .allocateWorkSpace(0, 0, 0, 0)
    ;
 }


DGKTSDDataBaseNSOR& DGKTSDDataBaseNSOR::setDefaultParameter(const LongReal& normA)
 {
    (*this)
     .setNormA2(normA*normA)
     .setEpsilon(1.0e-1)
     .setAccuracy(1.0e-4)
     .setMaxIterations(1000)
     .setError(1.0)
     .setAlpha(1.0)
    ;              
   
   return (*this);
 }


DGKTSDDataBaseNSOR& DGKTSDDataBaseNSOR::allocateWorkSpace(const LongInt& d, const LongInt& k, const LongInt& r, const LongInt& m)
 {    
    (*this)
     .setD(d)
     .setR(r)
     .setK(k)
     .setM(m)
    ;    
    
    attr_gradient.resize(d, r, m);

    attr_xa.resize(d, r, k);
    attr_xx.resize(d, r, r);
    
    
    const LongInt SIPA = d*r*k;
    const LongInt SIPX = d*r*r;
    
    attr_valuesAu = (LongRealPointer) new LongReal[SIPA];
    attr_valuesAo = (LongRealPointer) new LongReal[SIPA];
    attr_valuesXu = (LongRealPointer) new LongReal[SIPX];
    attr_valuesXo = (LongRealPointer) new LongReal[SIPX];
    
   return (*this);
 }


DGKTSDDataBaseNSOR::~DGKTSDDataBaseNSOR()
 {
    (*this)    
     .deleteWorkSpace()
    ;
 }


DGKTSDDataBaseNSOR& DGKTSDDataBaseNSOR::deleteWorkSpace()
 {
    attr_gradient.resize(0, 0, 0);
    
    attr_xa.resize(0, 0, 0);
    attr_xx.resize(0, 0, 0);
    
    
    delete [] attr_valuesAu;
    delete [] attr_valuesAo;
    delete [] attr_valuesXu;
    delete [] attr_valuesXo;
    
    
   return (*this);
 }


bool DGKTSDDataBaseNSOR::resize(const DKTS& x, const DKTS& a)
 {
    bool  value = true;
/* 
    (*this)
     .setNormA2(innerProduct(a, a))
    ;
*/
    //! Parameter for x
    const LongInt d = x.d();
    const LongInt r = x.k();
    const LongInt m = x.n();
    
    
    //! Parameter for a
    const LongInt k = a.k();

    // old    
    //! Parameter for xOld
    const LongInt dOld = DGKTSDDataBaseNSOR::d();
    const LongInt rOld = DGKTSDDataBaseNSOR::r();
    const LongInt mOld = DGKTSDDataBaseNSOR::m();
    
    //! Parameter for a
    const LongInt kOld = DGKTSDDataBaseNSOR::k();
    
    if(d!=dOld || r!=rOld || m!=mOld || k!=kOld)
     {
        (*this) 
         .deleteWorkSpace()
         .allocateWorkSpace(d, k, r, m)
        ;
     }
    
   return value;
 }
   
