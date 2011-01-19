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

// IntegerArray.cpp: Implementierung der Klasse IntegerArray.
//
//////////////////////////////////////////////////////////////////////

#include "IntegerArray.hpp"


IntegerArray::IntegerArray(const LongInt& d)
 {
    (*this)
     .allocateDataSpace(d)
    ; 
    
    (*this)
     .setAllEntrys()
    ;
    
 }


IntegerArray::IntegerArray (const IntegerArray& v)
 {
    (*this)
     .allocateDataSpace(v.order())
    ;
    
    (*this) = v; 
 }

IntegerArray::IntegerArray(const LongInt& d, const LongInt& nMax)
 {
    (*this)
     .allocateDataSpace(d)
    ; 
    
    (*this)
     .setRand(nMax)
    ;
    
 }
 
     
IntegerArray::~IntegerArray()
 {
    (*this)
     .deleteDataSpace()
    ; 
 }


bool IntegerArray::operator == (const IntegerArray& a) const
 {
   return !((*this)!=a);
 }


bool IntegerArray::operator != (const IntegerArray& a) const
 {
    const LongInt n1 = (*this).order();
    const LongInt n2 = a.order();
 
    bool value = false;
 
    if(n1==n2)
     {
        for(LongInt mu=0; mu<n1 && value==false; mu++)
         {
            if((*this)(mu)!=a(mu))
             {
                value = true;
             }
         }
     }
    else
     {
        value = true;
     }
  
   return value;
 } 


IntegerArray& IntegerArray::operator = (const IntegerArray& a)
 {
    const LongInt n = a.order();
    
    (*this)
     .resize(n)
    ;
    
    for(LongInt i=0; i<n; i++)
     {
        (*this)(i) = a(i);
     }
    
   return (*this);
 }


IntegerArray& IntegerArray::operator *= (const LongInt& multiplier)
 {
    const LongInt d = order();
  
    for(LongInt mu=0; mu<d; mu++)
     {
        (*this)(mu) *= multiplier;
     }
  
   return (*this);
 }


LongInt IntegerArray::sumOfAllEntry2(const LongInt& nu) const
 {
    LongInt sum = 0;
    
    for(LongInt mu=0; mu<nu; mu++)
     {
        sum += (*this)(mu);
     }
     
   return sum;
 }


IntegerArray& IntegerArray::resize(const LongInt& d)
 {    
    const LongInt orderOld = IntegerArray::order();
                
    if(orderOld!=d)
     {
        (*this)
         .deleteDataSpace()
         .allocateDataSpace(d)
        ;
     }
     
    (*this)
     .setAllEntrys()
    ;

   return (*this);
 }


IntegerArray& IntegerArray::allocateDataSpace(const LongInt& d)
 {
    (*this)
     .setOrder(d)     
    ; 

     attr_elements = (LongIntPointer) new LongInt[attr_order];

   return (*this);
 }


IntegerArray& IntegerArray::deleteDataSpace()
 {
    (*this)
     .setOrder(0)
    ;

    delete [] attr_elements;

   return (*this);
 }


IntegerArray& IntegerArray::setAllEntrys(const LongInt& value)
 {
    const LongInt order = IntegerArray::order();
    
    for(LongInt l=0; l<order; l++)
     {
        (*this)(l) = value;
     }
 
   return (*this);
 }


IntegerArray& IntegerArray::setRand(const LongInt& vMax, const LongInt& vMin)
 {
    const LongInt order = IntegerArray::order();
    Random rand;
    
    for(LongInt l=0; l<order; l++)
     {
        (*this)(l) = rand.randomLongInt(vMax, vMin);
     }    
 
   return (*this);
 }


IntegerArray& IntegerArray::setPwProductOf (const LongInt& alpha, const IntegerArray& v1, const IntegerArray& v2)
 {
    const LongInt d = MIN(v1.order(), v2.order());

    (*this).resize(d);

    for(LongInt mu=0; mu<d; mu++)
     {
        (*this)(mu) = alpha*v1(mu)*v2(mu);
     }

   return (*this);
 }


ostream& operator << (ostream& s, const IntegerArray& A)
 {
    const LongInt order = A.order();
    const LongInt lMax  = order-1;
    
    s << "(";        
    for(LongInt l=0; l<order; l++)
     {
        s << A(l);
        
        if(l<lMax)
         {
            s << ", ";
         }
        else
         {
            s << ")";
         }        
     }
   return s;
 }
 
 
LongInt l1(const IntegerArray& v) 
 {
   LongInt value = 0;
 
   const LongInt d = v.order();
   
   for(LongInt i=0; i<d; i++)
    {
       value += v(i);
    }
 
   return value;
 }


LongInt maxEntryOf(const IntegerArray& v)
 {
    const LongInt d = v.order();
    
    LongInt value = v(0);
 
    for(LongInt mu=1; mu<d; mu++)
     {
        const LongInt& v_mu = v(mu);
        
        if(value<v_mu)
         {
            value = v_mu;
         }
     }
 
   return value;
 }
