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

// DGKTSDDataBaseCG.cpp: Implementierung der Klasse DGKTSDDataBaseCG.
//
//////////////////////////////////////////////////////////////////////

#include "DGKTSDDataBaseCG.hpp"


DGKTSDDataBaseCG::DGKTSDDataBaseCG()
 {
    (*this)
     .setDefaultParameter(1.0)
     .allocateWorkSpace(0, 0, 0, 0)
    ;
 }


DGKTSDDataBaseCG& DGKTSDDataBaseCG::setDefaultParameter(const LongReal& normA)
 {
    (*this)
     .setNormA2(normA*normA)
     .setEpsilon(1.0e-2)
     .setAccuracy(1.0e-1)
     .setMaxIterations(500)
     .setMaxIterationsCR(250)
     .setEpsilonStepSize(1.0e-7)
     .setC(0.600)
     .setDD(0.600)
     .setError(1.0)
     .setAlpha(1.0e-15)
    ;              
   
   return (*this);
 }


DGKTSDDataBaseCG& DGKTSDDataBaseCG::allocateWorkSpace(const LongInt& d, const LongInt& k, const LongInt& r, const LongInt& m)
 {    
    (*this)
     .setD(d)
     .setR(r)
     .setK(k)
     .setM(m)
    ;    
    
    attr_gradient.resize(d, r, m);
    attr_gradientS.resize(d, r, m);
    attr_gradientOld.resize(d, r, m);
    attr_direction.resize(d, r, m);
    attr_work.resize(d, r, m);

    attr_xa.resize(d, r, k);    
    attr_da.resize(d, r, k);

    attr_xx.resize(d, r, r);
    attr_xd.resize(d, r, r);
    attr_dd.resize(d, r, r);
    
    attr_v1.resize(k);
    attr_v2.resize(r);

    const LongInt dim_G = d*r*r;
    
    attr_inversOfA = (LongRealPointer)  new LongReal[dim_G];
    
    const LongInt SIPA = d*r*k;
    const LongInt SIPX = d*r*r;
    
    attr_valuesAu = (LongRealPointer) new LongReal[SIPA];
    attr_valuesAo = (LongRealPointer) new LongReal[SIPA];
    attr_valuesXu = (LongRealPointer) new LongReal[SIPX];
    attr_valuesXo = (LongRealPointer) new LongReal[SIPX];
    
   return (*this);
 }


DGKTSDDataBaseCG::~DGKTSDDataBaseCG()
 {
    (*this)    
     .deleteWorkSpace()
    ;
 }


DGKTSDDataBaseCG& DGKTSDDataBaseCG::deleteWorkSpace()
 {
    attr_gradient.resize(0, 0, 0);
    attr_gradientOld.resize(0, 0, 0);
    attr_direction.resize(0, 0, 0);
    attr_gradientS.resize(0, 0, 0);
    attr_work.resize(0, 0, 0);        
    
    attr_xa.resize(0, 0, 0);
    attr_da.resize(0, 0, 0);
    attr_xx.resize(0, 0, 0);
    attr_xd.resize(0, 0, 0);
    attr_dd.resize(0, 0, 0);
    
    attr_v1.resize(0);
    attr_v2.resize(0);
    
    delete [] attr_inversOfA;
    
    delete [] attr_valuesAu;
    delete [] attr_valuesAo;
    delete [] attr_valuesXu;
    delete [] attr_valuesXo;
    
    
   return (*this);
 }


bool DGKTSDDataBaseCG::resize(const DKTS& x, const DKTS& a)
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
    const LongInt dOld = DGKTSDDataBaseCG::d();
    const LongInt rOld = DGKTSDDataBaseCG::r();
    const LongInt mOld = DGKTSDDataBaseCG::m();
    
    //! Parameter for a
    const LongInt kOld = DGKTSDDataBaseCG::k();
    
    if(d!=dOld || r!=rOld || m!=mOld || k!=kOld)
     {
        (*this) 
         .deleteWorkSpace()
         .allocateWorkSpace(d, k, r, m)
        ;
     }
    
   return value;
 }
   
