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

// LongRealQR.cpp: Implementierung der Klasse LongRealQR.
//
//////////////////////////////////////////////////////////////////////

#include "LongRealQR.hpp"
#include "LapackInterface2.hpp"

LongRealQR::LongRealQR(const RMatrix& A, const LongInt& k)
:attr_Q(A.numberOfRows(), MIN(k,A.numberOfRows()) ),
 attr_R(MIN(k,A.numberOfRows()), A.numberOfColumns())
 {    
    (*this)
     .setOrder(MIN(k,A.numberOfRows()))
    ;

    const LongInt  m = MIN(A.numberOfRows(),k);
    const LongReal e = 1.0;

    for(LongInt i=0; i<m; i++)
     {
        attr_Q(i,i) = e;
     }

    qr(A); 
 }


LongRealQR::LongRealQR (const RMatrix& A)
:attr_Q(A.numberOfRows(), A.numberOfRows()),
 attr_R(A.numberOfRows(), A.numberOfColumns())
 {
    (*this)
     .setOrder(A.numberOfColumns())
    ;

    const LongInt  m = A.numberOfRows();
    const LongReal e = 1.0;

    for(LongInt i=0; i<m; i++)
     {
        attr_Q(i,i) = e;
     }

    qr(A); 
 }


LongRealQR::~LongRealQR()
 {
 }


LongRealQR& LongRealQR::qr(const RMatrix& A)
 {
    integer m     = A.numberOfRows();
    integer n     = A.numberOfColumns(); 
    integer dim   = MIN(n,m);
    integer lda   = m;
    integer lwork = MAX(m,n);
    integer info  = 0;

    RVector work(lwork), tau(dim); 

    //dgeqrf ( &m, &n, &A(0,0), &lda, &tau(0), 
    //          &work(0), &lwork, &info);   
    info = TensorCalculus::Lapack<double>::geqrf(m, n, &A(0,0), lda, &tau(0), &work(0), lwork);   

    generateRby(A);
 
    integer mQ   = attr_Q.numberOfRows();
    integer nQ   = attr_Q.numberOfColumns();
    integer dimQ = nQ;

    // dorgqr ( &mQ, &nQ, &dimQ, &A(0,0), &lda, &tau(0), 
    //          &work(0), &lwork, &info);
    info = TensorCalculus::Lapack<double>::orgqr( mQ, nQ, dimQ, &A(0,0), lda, &tau(0), &work(0), lwork);

    generateQby(A);    
    
   return (*this);
 }


LongRealQR& LongRealQR::generateQby (const RMatrix& A)
 {
    const LongInt m  = attr_Q.numberOfRows();
    const LongInt n  = MIN(attr_Q.numberOfColumns(),A.numberOfColumns());

    for(int i=0; i<m; i++)
     {
        for(int j=0; j<n; j++)
         {
            attr_Q(i,j) = A(i,j);
         }         
     }

   return (*this);
 }


LongRealQR& LongRealQR::generateRby (const RMatrix& A)
 {
    const LongInt n = attr_R.numberOfRows();
    const LongInt m = attr_R.numberOfColumns();

    for(int i=0; i<n; i++)
     {
        for(int j=i; j<m; j++)
         {
            attr_R(i,j) = A(i,j);
         }
     }

   return (*this);
 }


ostream& operator << (ostream& s, const LongRealQR &qr)
 {
    s << qr.Q() << endl << qr.R() << endl;
  
   return s;
 }


istream& operator >> (istream& s, LongRealQR &qr)
 {
    s >> qr.attr_Q;
    s >> qr.attr_R;

   return s;
 }
