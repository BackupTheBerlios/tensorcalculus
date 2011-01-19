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

// ComplexSvd.cpp: Implementierung der Klasse ComplexSvd.
//
//////////////////////////////////////////////////////////////////////

#include "ComplexSvd.hpp"
#include "LapackInterface2.hpp"

const char   ComplexSvd::job[] = "S";

ComplexSvd::ComplexSvd(const CMatrix& A, const LongInt k)
: attr_U(A.numberOfRows(), MIN(k, MIN(A.numberOfRows(), A.numberOfColumns()))), 
  attr_V(A.numberOfColumns(), MIN(k, MIN(A.numberOfRows(), A.numberOfColumns()))), 
  attr_sigma(MIN(A.numberOfColumns(), A.numberOfRows()))
 {
    (*this)
    ;   
    svd(A); 
 }


ComplexSvd::~ComplexSvd()
 {
 }

ComplexSvd& ComplexSvd::svd (const CMatrix& A)
 {
    integer m      = A.numberOfRows();
    integer n      = A.numberOfColumns();
    integer k      = MIN(n, m);
    integer maxDim = MAX(n, m);
    integer lwork  = maxDim + 5*k + k*k;
    integer info   = 0;


    CVector work(lwork);
    RVector rwork(5*maxDim);
    CMatrix U(m, k);
    CMatrix V(k, n);

    // zgesvd ( job, job, &m, &n, 
    //          (_MKL_Complex16*) &A(0,0), &m, &attr_sigma(0), 
    //          (_MKL_Complex16*) &U(0,0), &m, 
    //          (_MKL_Complex16*) &V(0,0), &k, 
    //          (_MKL_Complex16*) &work(0), &lwork, 
    //          &rwork(0), &info);
    info = TensorCalculus::Lapack< std::complex<double> >::gesvd (*job, *job, m, n,
				   &A(0,0), m, &attr_sigma(0), 
				   &U(0,0), m, &V(0,0), k, &work(0), lwork, &rwork(0));

    generateUby(U);
    generateVby(V);


   return (*this);
 }

ComplexSvd& ComplexSvd::multiplyUSigma()
 {
    LongInt n = attr_U.numberOfRows();
    LongInt m = attr_U.numberOfColumns();

    for(LongInt i=0; i<n; i++)
     {
        for(LongInt j=0; j<m; j++)
         {
            attr_U(i,j) *=attr_sigma(j);
         }
     }

   return (*this);
 }


ComplexSvd& ComplexSvd::generateUby (const CMatrix& U)
 {
    LongInt m = attr_U.numberOfRows();
    LongInt n = attr_U.numberOfColumns();

    for(LongInt i=0; i<m; i++)
     {
        for(LongInt j=0; j<n; j++)
         {
            attr_U(i,j) = U(i,j);
         }
     }

   return (*this);
 }


ComplexSvd& ComplexSvd::generateVby (const CMatrix& V)
 {
    LongInt m = attr_V.numberOfRows();
    LongInt n = attr_V.numberOfColumns();
    
    for(LongInt i=0; i<m; i++)
     {
        for(LongInt j=0; j<n; j++)
         {
            attr_V(i,j) = std::conj(V(j,i));
         }
     }

   return (*this);
 }
