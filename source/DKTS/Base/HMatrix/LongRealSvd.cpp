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

// LongRealSvd.cpp: Implementierung der Klasse LongRealSvd.
//
//////////////////////////////////////////////////////////////////////

#include "LongRealSvd.hpp"
#include "LapackInterface2.hpp"

const char   LongRealSvd::job[] = "S";

LongRealSvd::LongRealSvd(const RMatrix& A, const LongInt k)
: attr_U(A.numberOfRows(), MIN(k, MIN(A.numberOfRows(), A.numberOfColumns()))), 
  attr_V(A.numberOfColumns(), MIN(k, MIN(A.numberOfRows(), A.numberOfColumns()))), 
  attr_sigma(MIN(A.numberOfColumns(), A.numberOfRows()))
 {
    (*this)
    ;   
    svd(A); 
 }


LongRealSvd::LongRealSvd(const RMatrix& A)
: attr_U(A.numberOfRows(), MIN(A.numberOfRows(), A.numberOfColumns())), 
  attr_V(A.numberOfColumns(), MIN(A.numberOfRows(), A.numberOfColumns())), 
  attr_sigma(MIN(A.numberOfColumns(), A.numberOfRows()))
 {
    (*this)
    ;   
    svd(A); 
 }

LongRealSvd::~LongRealSvd()
 {

 }

LongRealSvd& LongRealSvd::svd (const RMatrix& A)
 {
    integer m      = A.numberOfRows();
    integer n      = A.numberOfColumns();
    integer k      = MIN(n, m);
    integer maxDim = MAX(n, m);
    integer lwork  = maxDim + 5*k + k*k;
    integer info   = 0;


    RVector work(lwork);
    RMatrix U(m, k);
    RMatrix V(k, n);

    // dgesvd ( job, job, &m, &n, 
    //          (LongReal*) &A(0,0), &m, &attr_sigma(0), 
    //          (LongReal*) &U(0,0), &m, 
    //          (LongReal*) &V(0,0), &k, 
    //          (LongReal*) &work(0), &lwork, 
    //           &info);
    info = TensorCalculus::Lapack<double>::gesvd(*job, *job, m, n, &A(0,0), m, &attr_sigma(0), 
             &U(0,0), m, &V(0,0), k, &work(0), lwork);

    generateUby(U);
    generateVby(V);

   return (*this);
 }


LongRealSvd& LongRealSvd::multiplyUSigma()
 {
    LongInt n = attr_U.numberOfRows();
    LongInt m = attr_U.numberOfColumns();

    for(LongInt i=0; i<n; i++)
     {
        for(LongInt j=0; j<m; j++)
         {
            attr_U(i,j) *= attr_sigma(j);
         }
     }

   return (*this);
 }


LongRealSvd& LongRealSvd::generateUby (const RMatrix& U)
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


LongRealSvd& LongRealSvd::generateVby (const RMatrix& V)
 {
    LongInt m = attr_V.numberOfRows();
    LongInt n = attr_V.numberOfColumns();
    
    for(LongInt i=0; i<m; i++)
     {
        for(LongInt j=0; j<n; j++)
         {
            attr_V(i,j) = V(j,i);
         }
     }

   return (*this);
 }
