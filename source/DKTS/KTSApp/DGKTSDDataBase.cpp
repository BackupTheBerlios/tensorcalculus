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

// DGKTSDDataBase.cpp: Implementierung der Klasse DGKTSDDataBase.
//
//////////////////////////////////////////////////////////////////////

#include "DGKTSDDataBase.hpp"


DGKTSDDataBase::DGKTSDDataBase()
:attr_gradient(0, 0, 0), attr_direction(0, 0, 0)
 {
    (*this)
     .setDefaultParameter()
     .allocateWorkSpace(0, 0, 0, 0)
    ;
 }


DGKTSDDataBase& DGKTSDDataBase::setDefaultParameter()
 {
    (*this)
     .setEpsilon(1.0e-4)
     .setMaxSteps(350)
     .setEpsilonCR(1.0e-10)
     .setMaxStepsCR(250)
     .setMaxStepsAmijo(32)
     .setSigma(5.0e-4)
     .setPrecision(1.0e-3)
     .setMu(1.0)
     .setKappa(0.0)
     .setKappa1(1.0)
     .setKappa2(1.0)
     .setKappa3(1.0)
     .setLambda1(1.0e-3)
     .setLambda2(0.0)
     .setUpperBound(1.0e20)
     .setPreC(6)
     .setPrintCout(false)
     .setUseNewtonMethode(true)
     .setGamma(0.9)
     .setDelta(1.0e-1)     
     .setMaxStepsCR(350)
     .setMaxStepsCRL1(20)
     .setMaxStepsCRL2(5)     
    ;

    attr_error = 1.0;
   
   return (*this);
 }


DGKTSDDataBase& DGKTSDDataBase::allocateWorkSpace (const LongInt& d, const LongInt& k, const LongInt& l, const LongInt& m)
 {
    (*this)
     .setD(d)
     .setK(k) 
     .setL(l)
     .setM(m)
    ; 

    attr_dl = d*l;

    const LongInt dim_a  = d*k*l;
    const LongInt dim_x  = (LongInt)((d*k*(k+1))*0.5);
    const LongInt dim_xd = k*k*d;

    
    attr_workK = (LongRealPointer) new LongReal [k];
    attr_workL = (LongRealPointer) new LongReal [l];

    attr_a  = (LongRealPointer) new LongReal [dim_a];
    attr_x  = (LongRealPointer) new LongReal [dim_x];

    attr_ad = (LongRealPointer) new LongReal [dim_a];
    attr_xd = (LongRealPointer) new LongReal [dim_xd];
    attr_dd = (LongRealPointer) new LongReal [dim_x];

    const LongInt  dd      = (LongInt)(d*(d+1)*0.5);
    const LongInt dim_phi = k*l*dd;
    const LongInt dim_psi = k*k*dd;

    attr_phi = (LongRealPointer) new LongReal [dim_phi];
    attr_psi = (LongRealPointer) new LongReal [dim_psi];


    const LongInt lwork = MAX(1, 3*k-1);

    attr_eigWork  = (LongRealPointer) new LongReal [lwork];
    attr_eigValue = (LongRealPointer) new LongReal [k];
    attr_eigA     = (LongRealPointer) new LongReal [k*k];


    const LongInt dim_G = d*k*k;
    attr_A = (LongRealPointer) new LongReal [dim_G];

    const LongInt dim_W = d*k;
    attr_W = (LongRealPointer) new LongReal [dim_W];

    attr_gradient.resize(d, k, m);
    attr_gradientV.resize(d, k, m);
    attr_direction.resize(d, k, m);
    attr_work.resize(d, k, m);

    attr_w.resize(d, k, m);
    attr_r.resize(d, k, m);

    attr_ap.resize(d, k, m);
    attr_apo.resize(d, k, m);

    attr_b.resize(d, k, m);
    attr_bo.resize(d, k, m);

    attr_p.resize(d, k, m);
    attr_po.resize(d, k, m);

    const LongInt mem = 2*(dim_a + dim_x) + dim_xd + dim_phi + dim_psi + dim_G + 12*d*k*m + dim_W + 2*k + lwork + k*k + l;

    (*this)
     .setMemory(mem)
    ;

   return (*this);
 }


DGKTSDDataBase::~DGKTSDDataBase()
 {
    deleteWorkSpace();
 }


DGKTSDDataBase& DGKTSDDataBase::deleteWorkSpace()
 {
 
    delete [] attr_workK;
    delete [] attr_workL;

    delete [] attr_eigWork;
    delete [] attr_eigValue;
    delete [] attr_eigA;

    delete [] attr_a;
    delete [] attr_x;

    delete [] attr_ad;
    delete [] attr_dd;
    delete [] attr_xd;

    delete [] attr_phi;
    delete [] attr_psi;
 
    delete [] attr_A;
    delete [] attr_W;

   return (*this);
 }


bool DGKTSDDataBase::resize(const DKTS& a, DKTS& xi)
 {
    bool  value = false;

    const LongInt da = a.d();
    const LongInt dx = xi.d();

    const LongInt na = a.n();
    const LongInt nx = xi.n();

    if(da==dx && na==nx)
     {
        value = true;

        const LongInt ds  = dx;
        const LongInt ks  = xi.k();
        const LongInt ls  = a.k();
        const LongInt ms  = nx;

        const LongInt di  = d();
        const LongInt ki  = k();
        const LongInt li  = l();
        const LongInt mi  = m();

        if(ds!=di || ks!=ki || ls!=li || ms!=mi)
         {
            (*this) 
             .deleteWorkSpace()
             .allocateWorkSpace(ds, ks, ls, ms)
            ;
         }
     }

   return value;
 }


LongInt DGKTSDDataBase::isIndexOffsetOfa(const LongInt& j, const LongInt& i) const
 {
    LongInt index = j*attr_dl+i*attr_d;

   return index;
 }


LongInt DGKTSDDataBase::isIndexOffsetOfx(const LongInt& j1, const LongInt& j2) const
 {
    const LongInt j = MAX(j1,j2);

    LongInt index = attr_d*( (LongInt)(0.5*(j*(j+1))) + MIN(j1,j2) );
 
   return index;
 }


LongInt DGKTSDDataBase::isIndexOffsetOfxd(const LongInt& j1, const LongInt& j2) const
 {
    LongInt index = attr_d*( attr_k*j1 + j2);
 
   return index;
 }

LongReal& DGKTSDDataBase::A(const LongInt& j1, const LongInt& j2, const LongInt& mu) const
 {
    const LongInt index = (MAX(j1,j2) + mu*attr_k)*attr_k + MIN(j1,j2);
 
   return attr_A[index];
 }


LongReal& DGKTSDDataBase::W(const LongInt& j,  const LongInt& mu) const
 {
    const LongInt index = j + mu*k();
 
   return attr_W[index];
 }


LongReal& DGKTSDDataBase::A(const LongInt& mu) const
 {
    const LongInt index = mu*attr_k*attr_k;
 
   return attr_A[index];
 }


LongReal& DGKTSDDataBase::ad(const LongInt& j, const LongInt& i, const LongInt& mu) const
 {
    LongInt index = isIndexOffsetOfa(j, i);
    
    index += mu;

   return attr_ad[index];
 }


LongReal& DGKTSDDataBase::dd(const LongInt& j1, const LongInt& j2, const LongInt& mu) const
 {
    LongInt index = isIndexOffsetOfx(j1, j2);
    
    index += mu;

   return attr_dd[index];
 }


LongReal& DGKTSDDataBase::xd(const LongInt& j1, const LongInt& j2, const LongInt& mu) const
 {
    LongInt index = isIndexOffsetOfxd(j1, j2);
    
    index += mu;

   return attr_xd[index];
 }


LongReal& DGKTSDDataBase::x(const LongInt& j1, const LongInt& j2, const LongInt& mu) const
 {
    LongInt index = isIndexOffsetOfx(j1, j2);
    
    index += mu;

   return attr_x[index];
 }


LongReal& DGKTSDDataBase::a(const LongInt& j, const LongInt& i, const LongInt& mu) const
{
    LongInt index = isIndexOffsetOfa(j, i);
    
    index += mu;

   return attr_a[index];
 }


LongReal& DGKTSDDataBase::psi(const LongInt& j1, const LongInt& j2, const LongInt& mu1, const LongInt& mu2) const
 {
    const LongInt index_jj = j1*attr_k + j2;
    const LongInt index2   = isIndexOffsetOfmuX(mu1, mu2);

   return attr_psi[index_jj+index2];
 }


LongInt DGKTSDDataBase::isIndexOffsetOfmuX(const LongInt& mu1, const LongInt& mu2) const
 {
    const LongInt mu = MAX(mu1,mu2);

    LongInt index = ( (LongInt)(0.5*((mu)*(mu+1))) + MIN(mu1,mu2))*attr_k*attr_k;
 
   return index;
 }


LongInt DGKTSDDataBase::isIndexOffsetOfmuA(const LongInt& mu1, const LongInt& mu2) const
 {
    const LongInt mu = MAX(mu1,mu2);

    LongInt index = ( (LongInt)(0.5*((mu)*(mu+1))) + MIN(mu1,mu2))*attr_l*attr_k;
 
   return index;
 }


LongReal& DGKTSDDataBase::phi(const LongInt& j, const LongInt& i, const LongInt& mu1, const LongInt& mu2) const
 {
    const LongInt index_ji = (j*attr_l+i);
    const LongInt index2   = isIndexOffsetOfmuA(mu1, mu2);

   return attr_phi[index_ji+index2];
 }
