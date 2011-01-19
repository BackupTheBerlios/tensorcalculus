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

// ComplexQR.cpp: Implementierung der Klasse ComplexQR.
//
//////////////////////////////////////////////////////////////////////

#include "ComplexQR.hpp"
#include "LapackInterface2.hpp"

ComplexQR::ComplexQR(const CMatrix& A, const LongInt& k)
:attr_Q(A.numberOfRows(), MIN(k, A.numberOfRows()) ),
 attr_R(MIN(k,A.numberOfRows()), A.numberOfColumns())
 {
    
    (*this)
     .setOrder(MIN(k,A.numberOfRows()))
    ;

    const LongInt m = MIN(A.numberOfRows(),k);
    const std::complex<double> e(1.0, 0.0);

    for(LongInt i=0; i<m; i++)
     {
        attr_Q(i,i) = e;
     }

    qr(A); 
 }


ComplexQR::~ComplexQR()
 {
 }


ComplexQR& ComplexQR::qr (const CMatrix& A)
 {
    integer m     = A.numberOfRows();
    integer n     = A.numberOfColumns(); 
    integer dim   = MIN(n,m);
    integer lda   = m;
    integer lwork = MAX(m,n);
    integer info  = 0;

    CVector work(lwork), tau(dim); 

    // zgeqrf ( &m, &n, (_MKL_Complex16*) &A(0,0), &lda, (_MKL_Complex16*) &tau(0), 
    //          (_MKL_Complex16*) &work(0), &lwork, &info);   
    info = TensorCalculus::Lapack< std::complex<double> >::geqrf (m, n, &A(0,0), lda, &tau(0), &work(0), lwork);

    generateRby(A);
 
    integer mQ   = attr_Q.numberOfRows();
    integer nQ   = attr_Q.numberOfColumns();
    integer dimQ = nQ;


    // zungqr ( &mQ, &nQ, &dimQ, (_MKL_Complex16*) &A(0,0), &lda, (_MKL_Complex16*) &tau(0), 
    //          (_MKL_Complex16*) &work(0), &lwork, &info);
    info = TensorCalculus::Lapack< std::complex<double> >::ungqr(mQ, nQ, dimQ, &A(0,0), lda, &tau(0), &work(0), lwork);

    generateQby(A);    
    
   return (*this);
 }



ComplexQR& ComplexQR::generateQby (const CMatrix& A)
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


ComplexQR& ComplexQR::generateRby (const CMatrix& A)
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


ostream& operator << (ostream& s, const ComplexQR &qr)
 {
    s << qr.attr_Q << endl << qr.attr_R;
  
   return s;
 }


istream& operator >> (istream& s, ComplexQR &qr)
 {
    s >> qr.attr_Q;
    s >> qr.attr_R;

   return s;
 }
