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

// DKTVector.cpp: Implementierung der Klasse DKTVector.
//
//////////////////////////////////////////////////////////////////////

#include "DKTVector.hpp"
#include "BlasInterface.hpp"


LongInt DKTVector::inc = 1;


DKTVector::DKTVector()
:KTVector(NULL), attr_element(NULL)
 {
 
 }

DKTVector::DKTVector(LongIntPointer n, LongRealPointer element)
:KTVector(n), attr_element(element)
 {
 
 }
 
 
DKTVector::~DKTVector()
 {
 
 }



DKTVector& DKTVector::set(LongInt& n, LongReal& element)
 {
    attr_n       = &n;
    attr_element = &element;
    
   return (*this);
 }
 

DKTVector& DKTVector::operator = (const DKTVector& x)
 {  
    const LongInt& n = DKTVector::n();   
    
    TensorCalculus::Blas<double>::copy(n, &x(), inc, &(*this)(), inc);
    
   return (*this);
 }


DKTVector& DKTVector::setCopyOf(const LongReal& alpha, const DKTVector& x)
 {
    (*this)  = x;
    (*this) *= alpha;
    
   return (*this); 
 }
 
 
DKTVector& DKTVector::operator += (const DKTVector& x)
 {
    LongReal e = 1.0; 
    
    const LongInt& n = DKTVector::n();    
    TensorCalculus::Blas<double>::axpy(n, e, &x(), inc, &(*this)(), inc);
    
   return (*this);
 } 

 
DKTVector& DKTVector::operator -= (const DKTVector& x)
 {
    LongReal e = -1.0;
    const LongInt& n = DKTVector::n();    
    
    TensorCalculus::Blas<double>::axpy(n, e, &x(), inc, &(*this)(), inc);
    
   return (*this);
 } 

 
DKTVector& DKTVector::operator *= (const LongReal&  alpha)
 {  
    const LongInt& n = DKTVector::n();      
    TensorCalculus::Blas<double>::scal(n, alpha, &(*this)(), inc);        
    
   return (*this);
 }


DKTVector& DKTVector::operator /= (const LongReal&  alpha)
 {
    if(fabs(alpha)<EPS_NULL)
     {
        cout << "alpha = " << alpha << endl;
        throw SimpleException(IString("Warning In DKTVector::operator /= (const LongReal&  alpha), fabs(alpha)<EPS_NULL !!!"));
     }
    else
     {
        (*this) *= 1.0/alpha;
     }
 
   return (*this);
 }

//
// y<--alpha x + beta y
DKTVector& DKTVector::add(const LongReal& alpha, const DKTVector& x, const LongReal& beta)
 {
     const LongInt& n = DKTVector::n();    
     
     TensorCalculus::Blas<double>::scal(n, beta,  &(*this)(), inc);
     TensorCalculus::Blas<double>::axpy(n, alpha, &x(), inc, &(*this)(), inc);  
     
   return (*this);
 }

//
//y  <--  alpha*x + y
DKTVector& DKTVector::update(const LongReal& alpha, const DKTVector& x)
 {
    const LongInt& n = DKTVector::n();    
    TensorCalculus::Blas<double>::axpy(n, alpha, &x(), inc, &(*this)(), inc);
     
   return (*this);
 }


DKTVector& DKTVector::pwProduct(const LongReal& alpha, const DKTVector& x, const DKTVector& y)
 {
    const LongInt n = DKTVector::n(); 

    for(LongInt i=0; i<n; i++)
     {
        (*this)(i) = alpha*x(i)*y(i);
     }
    
   return (*this);
 }


DKTVector& DKTVector::setNull()
 {
    (*this)
     .setAllEntrys(0.0)
    ;
       
   return (*this);
 }


DKTVector& DKTVector::setAllEntrys(const LongReal& v)
{
    const LongInt n = DKTVector::n();
    
    for(LongInt i=0; i<n; i++)
     {
        (*this)(i) = v;
     }
   
   return (*this);
 }

DKTVector& DKTVector::setRand(const LongReal& eps)
 {
    const LongInt n = DKTVector::n();

    Random rand;

    for(LongInt i=0; i<n; i++)
     {
        (*this)(i) = rand.randomLongReal(eps);
     }
 
   return (*this);
 }


LongReal DKTVector::normalized()
 {
    LongReal value = l2((*this));
 
    (*this) /= value;
 
   return value;
 }


LongReal innerProduct (const DKTVector& a, const DKTVector& b)
 {
    LongInt n = MIN(a.n(), b.n());
    
    LongReal value = TensorCalculus::Blas<double>::dot(n, &a(), DKTVector::inc, &b(), DKTVector::inc);
         
   return value;
 }
 
 
LongReal l2(const DKTVector& x)
 {
   LongReal value = TensorCalculus::Blas<double>::nrm2(x.n(), &x(), DKTVector::inc);
     
   return value;
 }
 

LongInt DKTVector::indexOfMax() const
 {
    const LongInt n = DKTVector::n();
    
    LongInt index = 0;
    
    LongReal maxValue = fabs((*this)(0));  
    
    for(LongInt i=1; i<n; i++)
     {
        const LongReal temp = fabs((*this)(i));                        
        
        if(maxValue < temp)
         {
            maxValue = temp;
            index = i;
         }
     }
 
   return index;
 }


LongReal maximumNormOf (const DKTVector& x)
 {
    const LongInt i_max  = x.indexOfMax();    
    const LongReal value = fabs(x(i_max));
    
   return value;
 }

LongReal maximumValueOf (const DKTVector& x, LongInt& indexOfMaxValue)
 {
    const LongInt n = x.n();
    
    LongReal value  = x(0);
    indexOfMaxValue = 0;
    
    for(LongInt l=1; l<n; l++)
     {
        const LongReal& value_l = x(l);
        
        if(value < value_l)
         {
            value           = value_l;
            indexOfMaxValue = l;
         }
     }
    
   return value;
 }


LongReal minimumValueOf (const DKTVector& x, LongInt& indexOfMinValue)
 {
    const LongInt n = x.n();
    
    LongReal value  = x(0);
    indexOfMinValue = 0;
    
    for(LongInt l=1; l<n; l++)
     {
        const LongReal& value_l = x(l);
        
        if(value_l < value)
         {
            value           = value_l;
            indexOfMinValue = l;
         }
     }
    
   return value;
 }


ostream& operator << (ostream& os, const DKTVector& x)
 {
    const LongInt n = x.n();
  
    LongRealPointer pE = &x();
    
    os << n << endl;
    for(LongInt i=0; i<n; i++)
     {
        os << i << " " << pE[i] << endl;
     }
    
   return os;
 }
