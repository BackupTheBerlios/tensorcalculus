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

// CVector.cpp
//

#include "CVector.hpp"
#include "CMatrix.hpp"
#include "RkCMatrix.hpp"
#include "BlasInterface.hpp"

CVector::CVector(LongInt dim)
:AVector(dim), attr_refCounter(0)
 {
	_pkomp	= new std::complex<double> [attr_dimension];
 }


CVector::CVector(const CVector &v)
:AVector(v.dimension())
 {
	_pkomp	= new std::complex<double> [attr_dimension];

    (*this) = v;
 }


CVector::CVector(const RVector &v)
:AVector(v.dimension())
 {
	_pkomp	= new std::complex<double> [attr_dimension];

    const LongInt dim = v.dimension();

    for(LongInt i=0; i<dim; i++)
     {
        (*this)(i) = v(i);
     }

 }


CVector::CVector(const LongInt& dim, const std::complex<double>* add)
:AVector(dim), attr_refCounter(0)
 {
    attr_refCounter++;   
    _pkomp = (std::complex<double>*)add;
 }


CVector::~CVector()
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


CVector& CVector::update(const std::complex<double>& a, const CVector& d)
 {
    integer n  = attr_dimension;
    integer ic = 1;

    // zaxpy (&n, (_MKL_Complex16*)&a, (_MKL_Complex16*)&d(0), &ic, (_MKL_Complex16*) &(*this)(0), &ic);
    TensorCalculus::Blas< std::complex<double> >::axpy (n, a, &d(0), ic, &(*this)(0), ic);
    
   return (*this);
 }

std::complex<double> innerProduct(const CVector &v1, const CVector &v2)
 {
    integer n    = MIN(v1.dimension(), v2.dimension());
    integer inc  = 1;
    std::complex<double> temp = 0.0;

    // zdotc ((_MKL_Complex16*) &temp, &n, (_MKL_Complex16*)&v1(0), &inc, (_MKL_Complex16*) &v2(0), &inc);
    temp = TensorCalculus::Blas< std::complex<double> >::dotc(n, &v1(0), inc, &v2(0), inc);

   return temp;
 }


LongReal  L1(const CVector &v)
 {
    LongReal temp = 0.0;

    integer n  = v.dimension();
    integer ic = 1;

    // temp = dzasum (&n, (_MKL_Complex16*) &v(0), &ic);
    temp = TensorCalculus::Blas< std::complex<double> >::asum(n, &v(0), ic);

   return temp;
 }


LongReal  L2(const CVector &v)
 {
    LongReal temp = 0.0;

    integer n  = v.dimension();
    integer ic = 1;
 
    // temp = dznrm2 (&n, (_MKL_Complex16*) &v(0), &ic);
    temp = TensorCalculus::Blas< std::complex<double> >::nrm2(n, &v(0), ic);

   return temp;
 }


LongReal Max(const CVector &v)
 {
    LongReal max  = -1.0;
    LongReal temp =  0.0;
    const LongInt n = v.dimension();    

    for(int i=0; i<n; i++)
     {
        temp = std::norm(v(i));
        if(max<temp)
         {
            max = temp;
         }
     }

   return sqrt(max);     
 }


CVector& CVector :: operator = (const CVector& v)
 { 
    integer n  = v.dimension();
    integer m  = attr_dimension;
    integer nx = 1;

    if(n!=m)
     {
        delete []_pkomp;
        _pkomp = new std::complex<double> [n];
        attr_dimension   = n; 
     }

    attr_refCounter = v.refCounter();

    // zcopy (&n, (_MKL_Complex16*)(&v(0)), &nx, (_MKL_Complex16*)(&(*this)(0)), &nx);
    TensorCalculus::Blas< std::complex<double> >::copy (n, (&v(0)), nx, (&(*this)(0)), nx);

   return *this;
 }


CVector& CVector::operator += (const CVector& v)
 {
    update(std::complex<double>(1.0,0.0), v);
   return (*this);
 }


CVector& CVector::operator -= (const CVector& v)
 {
    update(std::complex<double>(-1.0,0.0), v);
   return (*this);
 }


CVector& CVector::operator *= (const std::complex<double>& a)
 {
    integer n  = attr_dimension;
    integer ic = 1;

    // zscal (&n, (_MKL_Complex16*) &a, (_MKL_Complex16*) &(*this)(0), &ic);
    TensorCalculus::Blas< std::complex<double> >::scal (n, a, &(*this)(0), ic);
   return (*this);
 }


CVector& CVector::operator /= (const std::complex<double>& a)
 {
    (*this) *= 1.0/a; 
   return (*this);
 }


CVector CVector :: operator + (const CVector &v) const
 {
    CVector sum((*this));
    
    sum += v;
    		
   return sum;
 }


CVector CVector :: operator - (const CVector& v) const
 {
    CVector diff((*this));
    
    diff -= v;
    		
   return diff;
 }

	
CVector  operator *  (std::complex<double> alpha, const CVector &v)
 {
    CVector result(v);
    
    result *= alpha;
    		
   return result;
 }


CVector  operator /  (const CVector &v, std::complex<double> alpha)
 {
    CVector result(v);
    
    result /= alpha;
    		
   return result;
 }


bool	CVector :: operator == (const CVector& v) const
 {
    if((*this).attr_dimension == v.attr_dimension)
     {
	    for (LongInt i=0 ; i < v.attr_dimension ; i++)
         {
			if((*this)(i) != v(i))return false;
		 }
	   return true;
	 }
	else
       return false;
 }

bool CVector :: operator != (const CVector& v) const
{
	if((*this).attr_dimension != v.attr_dimension)
     {
		for (LongInt i=0 ; i < v.attr_dimension ; i++)
         {
		    if((*this)(i) != v(i))return true;
		 }
	  return false;
	 }
	else
	  return true;
}


CVector& CVector::setNull ()
 {
    (*this) *= COMPLEX_NULL;

   return (*this);
 }

ostream&  operator << (ostream & s,const CVector &v)
 {
    LongInt dim = v.dimension();
	s << dim << '\n';
	for (LongInt i=0; i < dim ; i++){
		s << v._pkomp[i]<<'\t';
	}
	s << '\n';
   return s;
 }


istream&  operator >> (istream & s, CVector &v)
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


CVector& CVector ::newDimension(LongInt newdim)
 {
	delete []_pkomp;
	attr_dimension	= newdim;
	_pkomp	= new std::complex<double> [attr_dimension];
	for(LongInt i=0 ; i<attr_dimension ; i++)_pkomp[i] = (std::complex<double>)0;
   return (*this);
 }



RVector CVector::realPart () const
 {
    const LongInt dim = dimension();
    RVector value (dim);

    for(LongInt i=0 ; i<dim; i++)
     {
        value(i) = _pkomp[i].real();	
     }

   return value;
 }


RVector CVector::imagPart () const
 {
    const LongInt dim = dimension();
    RVector value (dim);
    
    for(LongInt i=0 ; i<dim; i++)
     {
        value(i) = _pkomp[i].imag();	
     }

   return value;
 }
