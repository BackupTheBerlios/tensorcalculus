/*
 * Copyright (C) Mike Espig, Henry Auer
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

// DKTSDSimpleInnerProducts.cpp: Implementierung der Klasse DKTSDSimpleInnerProducts.
//
//////////////////////////////////////////////////////////////////////

#include "DKTSDSimpleInnerProducts.hpp"
#include "BlasInterface.hpp"

DKTSDSimpleInnerProducts::DKTSDSimpleInnerProducts(const DKTS& x, const DKTS& y)
 {      
    (*this)
     .allocateDataSpace(MIN(x.d(), y.d()), x.k(), y.k())
    ; 
        
 }


DKTSDSimpleInnerProducts::DKTSDSimpleInnerProducts(const LongInt& d, const LongInt& r1, const LongInt& r2)
 {
    (*this)
     .allocateDataSpace(d, r1, r2)
    ;
 }

     
DKTSDSimpleInnerProducts::~DKTSDSimpleInnerProducts()
 {
    (*this)
     .deleteDataSpace()
    ; 
 }


DKTSDSimpleInnerProducts& DKTSDSimpleInnerProducts::resize(const DKTS& x, const DKTS& y)
 {    
    const LongInt d  = MIN(x.d(), y.d());
    const LongInt r1 = x.k();
    const LongInt r2 = y.k();

    (*this)
     .resize(d, r1, r2)
    ;                    

   return (*this);
 }


DKTSDSimpleInnerProducts& DKTSDSimpleInnerProducts::resize(const LongInt& d, const LongInt& r1, const LongInt& r2)
 {
    const LongInt dOld  = DKTSDSimpleInnerProducts::d();
    const LongInt r1Old = DKTSDSimpleInnerProducts::r1();
    const LongInt r2Old = DKTSDSimpleInnerProducts::r2();

    if(dOld!=d || r1Old!=r1 || r2Old!=r2)
     {
        (*this)
         .deleteDataSpace()
         .allocateDataSpace(d, r1, r2)
        ;
     }
     
    (*this)
     .setAllEntrys()
    ;
 
   return (*this);
 }


DKTSDSimpleInnerProducts& DKTSDSimpleInnerProducts::allocateDataSpace(const LongInt& d, const LongInt& r1, const LongInt& r2)
 {          
     (*this)
      .setD(d)
      .setR1(r1)
      .setR2(r2)
     ;

     const LongInt m = r1*r2;               
     const LongInt n = d*m;
                                   
     attr_valuesS = (LongRealPointer) new LongReal[n];     
     attr_valuesG = (LongRealPointer) new LongReal[n];

     attr_valuesAu = (LongRealPointer) new LongReal[n];     
     attr_valuesAo = (LongRealPointer) new LongReal[n];


     (*this)
      .setAllEntrys()
     ;

   return (*this);
 }


DKTSDSimpleInnerProducts& DKTSDSimpleInnerProducts::deleteDataSpace()
 {
    (*this)
     .setD(0)
     .setR1(0)
     .setR2(0)
    ;

    delete [] attr_valuesG;    
    delete [] attr_valuesS;    

    delete [] attr_valuesAu;    
    delete [] attr_valuesAo;    

   return (*this);
 }


DKTSDSimpleInnerProducts& DKTSDSimpleInnerProducts::setAllEntrys()
 {
    
    
   return (*this);
 }


DKTSDSimpleInnerProducts& DKTSDSimpleInnerProducts::computeInnerProducts(const DKTS& x, const DKTS& y) 
 {
    (*this)
     .resize(x, y)
    ;
    
    const LongInt d  = DKTSDSimpleInnerProducts::d();
    const LongInt r1 = DKTSDSimpleInnerProducts::r1();
    const LongInt r2 = DKTSDSimpleInnerProducts::r2();
    
    
    for(LongInt i1=0; i1<r1; i1++)
     {
        for(LongInt i2=0; i2<r2; i2++)
         {
            for(LongInt mu=0; mu<d; mu++)
             {
                (*this).innerProductAt(i1, i2, mu) = innerProduct(x(i1,mu), y(i2,mu));
             }
         }
     }
    
    
   return (*this);
 }


DKTSDSimpleInnerProducts& DKTSDSimpleInnerProducts::computeSimpleInnerProducts()
 {
    const LongInt d = DKTSDSimpleInnerProducts::d();
    const LongInt r1 = DKTSDSimpleInnerProducts::r1();
    const LongInt r2 = DKTSDSimpleInnerProducts::r2();

    for(LongInt mu=0; mu<d; mu++)
     {
        for(LongInt i1=0; i1<r1; i1++)
         {
            for(LongInt i2=0; i2<r2; i2++)
             {
                LongReal& value = tensorInnerProductAt(i1, i2, mu);

                value = 1.0;

                for(LongInt nu=0; nu<d; nu++)
                 {
                    if(mu!=nu)
                     {
                        value *=innerProductAt(i1, i2, nu);
                     }                    
                 }
             }
         }
     }

   return (*this);
 }
 

DKTSDSimpleInnerProducts& DKTSDSimpleInnerProducts::computeAllInnerProducts (const DKTS& x, const DKTS& y)
 {
    (*this)    
     .computeInnerProducts(x, y)     
    // .computeSimpleInnerProducts()
     .computeSimpleLeaveOut()
    ;
   return (*this);
 }

DKTSDSimpleInnerProducts& DKTSDSimpleInnerProducts::computeaU ()
 {
    LongInt d  = DKTSDSimpleInnerProducts::d();
    LongInt r1 = DKTSDSimpleInnerProducts::r1();
    LongInt r2 = DKTSDSimpleInnerProducts::r2();
    
    
    for (LongInt i1=0; i1 < r1; i1++)
     {
        for (LongInt i2=0; i2 < r2; i2++)
         {
            (*this).aU(i1,i2,0) = 1;
            for (LongInt mu=1; mu < d; mu++)
             {
                LongInt nu = mu-1;
                (*this).aU(i1,i2,mu) = (*this).aU(i1,i2,nu)*(*this).innerProductAt(i1,i2,nu);
             }
         }
     }
   return(*this);
 }

DKTSDSimpleInnerProducts& DKTSDSimpleInnerProducts::computeaO ()
 {
    LongInt d  = DKTSDSimpleInnerProducts::d();
    LongInt r1 = DKTSDSimpleInnerProducts::r1();
    LongInt r2 = DKTSDSimpleInnerProducts::r2();
    
    for (LongInt i1=0; i1 < r1; i1++)
     {
        for (LongInt i2=0; i2 < r2; i2++)
         {
            //Achtung for-Schleife r�ckw�rts
            (*this).aO(i1,i2,d-1) = 1;
            for (LongInt mu = d-2; 0 <= mu; mu--)
             {
                LongInt nu = mu+1;
                (*this).aO(i1,i2,mu) = (*this).aO(i1,i2,nu) * (*this).innerProductAt(i1,i2,nu);
             }
         }
     }
   return(*this);
 }
 
DKTSDSimpleInnerProducts& DKTSDSimpleInnerProducts::computeSimpleLeaveOut ()
 {
    LongInt  d = DKTSDSimpleInnerProducts::d();
    LongInt r1 = DKTSDSimpleInnerProducts::r1();
    LongInt r2 = DKTSDSimpleInnerProducts::r2();
    
    (*this)
    .computeaO()
    .computeaU();
    
    for(LongInt i1=0; i1 < r1; i1++)
     {
        for (LongInt i2=0; i2 < r2; i2++)
         {
            for (LongInt mu=0; mu < d; mu++)
             {
                (*this).tensorInnerProductAt(i1,i2,mu) = (*this).aU(i1,i2,mu) * (*this).aO(i1,i2,mu);
             }
         }
     }
   return (*this);
 }
 
DKTSDSimpleInnerProducts& DKTSDSimpleInnerProducts::operator *= (const LongReal& alpha)
 {
    //skaliert attr_valuesS[]
    const LongInt r1 = DKTSDSimpleInnerProducts::r1();
    const LongInt r2 = DKTSDSimpleInnerProducts::r2();
    const LongInt d  = DKTSDSimpleInnerProducts::d();
    
    const LongInt n          = r1*r2*d;
    const LongInt const_inc  = 1;
    
    LongReal& vec = (*this).innerProductAt(0,0,0);
    
    // dscal(&n, &alpha, &vec, &const_inc);
    TensorCalculus::Blas<double>::scal(n, alpha, &vec, const_inc);
    
   return(*this);
 }
