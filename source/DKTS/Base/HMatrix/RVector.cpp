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

// RVector.cpp
//

#include "RVector.hpp"
#include "RMatrix.hpp"
#include "RkRMatrix.hpp"
#include "BlasInterface.hpp"

RVector::RVector(LongInt dim)
:AVector(dim), attr_refCounter(0)
 {
	_pkomp	= new LongReal [attr_dimension];
    for(LongInt i=0; i<dim; i++)
     {
        _pkomp[i] = 0.0;
     }
 }


RVector::RVector(const RVector &v)
:AVector(v.dimension())
 {
	_pkomp	= new LongReal [attr_dimension];

    (*this) = v;
 }


RVector::RVector(const DKTVector& v)
:AVector(v.n())
 {    
    LongInt nx = 1;
    
    _pkomp	= new LongReal [attr_dimension];
    
    // dcopy (&attr_dimension, (LongReal*)(&v(0)), &nx, (LongReal*)(&(*this)(0)), &nx);
    TensorCalculus::Blas<double>::copy (attr_dimension, (&v(0)), nx, (&(*this)(0)), nx);
 }


RVector::RVector(const LongInt& dim, const LongReal* add)
:AVector(dim), attr_refCounter(0)
 {
    attr_refCounter++;   
    _pkomp = (LongReal*)add;
 }


RVector::~RVector()
 {
    if(attr_refCounter==0)
     {
 	    delete []_pkomp;
     }
    else
     {
        attr_refCounter--;
     }
 }


RVector& RVector::update(const LongReal& a, const RVector& d)
 {
    integer n  = attr_dimension;
    integer ic = 1;

    // daxpy (&n, (LongReal*)&a, (LongReal*)&d(0), &ic, (LongReal*) &(*this)(0), &ic);
    TensorCalculus::Blas<double>::axpy (n, a, &d(0), ic,  &(*this)(0), ic);
    
   return (*this);
 }


LongReal innerProduct(const RVector &v1, const RVector &v2)
 {
    integer n    = MIN(v1.dimension(), v2.dimension());
    integer inc  = 1;
    LongReal temp = 0.0;

    // temp = ddot (&n, (LongReal*)&v1(0), &inc, (LongReal*) &v2(0), &inc);
    temp = TensorCalculus::Blas<double>::dot(n, &v1(0), inc, &v2(0), inc);

   return temp;
 }


RVector&  RVector::setRand(const LongReal& eps)
 {
    const LongInt n = dimension();

    srand((unsigned)time(NULL));

    for(LongInt i=0; i<n; i++)
     {
        (*this)(i) = (LongReal)rand()*eps; 
     }

   return (*this);
 }


LongReal RVector::smallestEntry() const
 {
    const LongInt n = dimension();
    
    LongReal value = (*this)(0);

    for(LongInt i=0; i<n; i++)
     {
        value = MIN(value, (*this)(i));
     }
 
   return value;
 }


LongReal  L1(const RVector& v)
 {
    LongReal temp = 0.0;

    integer n  = v.dimension();
    integer ic = 1;

    // temp = dasum (&n, (LongReal*) &v(0), &ic);    
    temp = TensorCalculus::Blas<double>::asum (n, &v(0), ic);    

   return temp;
 }


LongReal  L2(const RVector &v)
 {
    LongReal temp = 0.0;

    integer n  = v.dimension();
    integer ic = 1;
 
    // temp = dnrm2 (&n, (LongReal*) &v(0), &ic);
    temp = TensorCalculus::Blas<double>::nrm2 (n, &v(0), ic);

   return temp;
 }


LongReal Max(const RVector &v)
 {
    LongReal max  = -1.0;
    LongReal temp =  0.0;
    const LongInt n = v.dimension();    

    for(int i=0; i<n; i++)
     {
        temp = fabs(v(i));
        if(max<temp)
         {
            max = temp;
         }
     }

   return max;     
 }


LongReal Min(const RVector &v)
 {
    LongReal min  = 1.0e30;
    LongReal temp =  0.0;
    
    const LongInt n = v.dimension();    

    for(LongInt i=0; i<n; i++)
     {
        temp = fabs(v(i));
        
        if(temp<min)
         {
            min = temp;
         }
     }

   return min;     
 }


RVector& RVector::operator = (const RVector& v)
 { 
    integer n  = v.dimension();
    integer m  = attr_dimension;
    integer nx = 1;

    if(n!=m)
     {
        delete []_pkomp;
        _pkomp = new LongReal [n];
        attr_dimension   = n; 
     }

    attr_refCounter = v.refCounter();

    // dcopy (&n, (LongReal*)(&v(0)), &nx, (LongReal*)(&(*this)(0)), &nx);
    TensorCalculus::Blas<double>::copy (n, (&v(0)), nx, (&(*this)(0)), nx);

   return *this;
 }


RVector& RVector::operator += (const RVector& v)
 {
    update(LongReal(1.0), v);
   return (*this);
 }


RVector& RVector::operator -= (const RVector& v)
 {
    update(LongReal(-1.0), v);
   return (*this);
 }


RVector& RVector::operator *= (const LongReal& a)
 {
    integer n  = attr_dimension;
    integer ic = 1;

    // dscal (&n, (LongReal*) &a, (LongReal*) &(*this)(0), &ic);
    TensorCalculus::Blas<double>::scal (n, a, &(*this)(0), ic);
   return (*this);
 }


RVector& RVector::operator /= (const LongReal& a)
 {
    (*this) *= 1.0/a; 
   return (*this);
 }


RVector RVector :: operator + (const RVector &v) const
 {
    RVector sum((*this));
    
    sum += v;
    		
   return sum;
 }


RVector RVector :: operator - (const RVector& v) const
 {
    RVector diff((*this));
    
    diff -= v;
    		
   return diff;
 }

	
RVector  operator *  (LongReal alpha, const RVector &v)
 {
    RVector result(v);
    
    result *= alpha;
    		
   return result;
 }


RVector  operator /  (const RVector &v, LongReal alpha)
 {
    RVector result(v);
    
    result /= alpha;
    		
   return result;
 }


bool	RVector :: operator == (const RVector& v) const
 {
    if((*this).attr_dimension == v.attr_dimension)
     {
	       for (LongInt i=0 ; i < v.attr_dimension ; i++)
         {
			         if((*this)(i) != v(i)) return false;
		       }
		       
	       return true;
	    }
	   
   return false;
 }


bool RVector::operator != (const RVector& v) const
 {
	   
	 return !((*this)==v);
}


RVector& RVector::setNull ()
 {
 
    for(LongInt i=0; i<attr_dimension; i++)
     {
        _pkomp[i] = 0.0;
     }

   return (*this);
 }

ostream&  operator << (ostream & s,const RVector &v)
 {
    LongInt dim = v.dimension();
	s << dim << '\n';
	for (LongInt i=0; i < dim ; i++){
		s << v._pkomp[i]<<'\t';
	}
	s << '\n';
   return s;
 }


istream&  operator >> (istream & s, RVector &v)
 {
    LongInt dim =0;
	   s >> dim;
	   v.newDimension(dim);
	   for (LongInt i=0; i < dim ; i++)
     {
		      s >> v._pkomp[i];
	    }
	    
  return s;
 }


RVector&  RVector ::resize(LongInt newdim)
 {
	   delete []_pkomp;
	
	   attr_dimension	= newdim;
	   _pkomp	= new LongReal [attr_dimension];
	   
	   for(LongInt i=0 ; i<attr_dimension ; i++)
	    {
	       _pkomp[i] = (LongReal)0;
	    }
	    
   return (*this);
 }


RVector& RVector ::newDimension(LongInt newdim)
 {
	   (*this)
	    .resize(newdim)
	   ;
	      
   return (*this);
 }

