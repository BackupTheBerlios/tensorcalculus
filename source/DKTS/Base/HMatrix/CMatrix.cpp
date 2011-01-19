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

// CMatrix.cpp: Implementierung der Klasse CMatrix.
//
//////////////////////////////////////////////////////////////////////

#include "CMatrix.hpp"
#include "RkCMatrix.hpp"
#include "BlasInterface.hpp"
#include "LapackInterface2.hpp"

CMatrix::CMatrix(LongInt m, LongInt n)
:ACMatrix(m, n)
 {
    (*this)
     .disableConjugateTranspose()     
     .setType(isCMatrix)
    ;

    _pelm =  (std::complex<double>*) new std::complex<double> [attr_numberOfRows*attr_numberOfColumns];
 }


CMatrix::~CMatrix()
 {
	delete [] _pelm;
 }


CMatrix::CMatrix(const CMatrix &A)
:ACMatrix(A.numberOfRows(), A.numberOfColumns())
 {   
    _pelm = (std::complex<double>*) new std::complex<double> [attr_numberOfRows*attr_numberOfColumns];

	(*this) = A;
 }


CMatrix::CMatrix(const RMatrix &A)
:ACMatrix(A.numberOfRows(), A.numberOfColumns())
 {   
    const LongInt m = A.numberOfRows();
    const LongInt n = A.numberOfColumns();   

    (*this)
     .disableConjugateTranspose()     
     .setType(isCMatrix)
    ;

    _pelm = (std::complex<double>*) new std::complex<double> [m*n];

    for(LongInt i=0; i<m; i++)
      {
         for(LongInt j=0; j<n; j++)
          {
             (*this)(i,j) = A(i,j);
          }
      }

 }


CMatrix& CMatrix :: operator = (const CMatrix& A)
 {	
    LongInt n = A.numberOfRows();
    LongInt m = A.numberOfColumns();   
    
    setDimension(n, m);
 
    attr_type       =  A.attr_type;
    setConjugateTranspose(A.isConjugateTranspose());
    setTranspose(A.isTranspose());

    integer dim = n*m;
    integer nx  = 1;

    // zcopy(&dim, (_MKL_Complex16*)(&A._pelm[0]), &nx, (_MKL_Complex16*)(&this->_pelm[0]), &nx);
    TensorCalculus::Blas< std::complex<double> >::copy(dim, &A._pelm[0], nx, &this->_pelm[0], nx);
    
   return (*this);
 }


CMatrix& CMatrix::update (const std::complex<double>& a, const CMatrix& B)
 {
    integer n  = attr_numberOfRows*attr_numberOfColumns;
    integer ic = 1;

    // zaxpy(&n, (_MKL_Complex16*)&a, (_MKL_Complex16*)&B._pelm[0], &ic, (_MKL_Complex16*) (&this->_pelm[0]), &ic);
    TensorCalculus::Blas< std::complex<double> >::axpy(n, a, &B._pelm[0], ic, &this->_pelm[0], ic);

   return (*this);
 }


CMatrix& CMatrix::operator += (const CMatrix& A)
 {
    update(COMPLEX_UNIT, A);

   return (*this);
 }


CMatrix& CMatrix::operator -= (const CMatrix& A)
 {
    update(-COMPLEX_UNIT, A);

   return (*this);
 }


CMatrix&  CMatrix::addIdentityScaledBy (const std::complex<double>& z)
 {
    const LongInt dim = MIN(numberOfRows(),numberOfColumns());
  
    for(LongInt i=0; i<dim; i++)
     {
        (*this)(i,i) += z;
     }

   return(*this);
 }


CMatrix& CMatrix::addMatrix (const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex)
 {
    const AMatrixType theMatrixA = A.type();

    LongInt m  = numberOfRows();
    LongInt n  = numberOfColumns();
    LongInt m1 = A.numberOfRows();
    LongInt n1 = A.numberOfColumns();
    
    if(theMatrixA==isRkCMatrix)
     {

        const RkCMatrix& R = (const RkCMatrix&)A;
        LongInt k = R.rank();

        LongInt mR = R.numberOfRows();
        LongInt nR = R.numberOfColumns(); 
      
        // zgemm ( &notConjTrans, &conjTrans, &MIN(m,mR), &MIN(n,nR), &k, (_MKL_Complex16*) &Complex::complexUnit,
	//             (_MKL_Complex16*) &(R.attr_A(rowIndex, 0)), &mR, 
        //         (_MKL_Complex16*) &(R.attr_B(colIndex, 0)), &nR, (_MKL_Complex16*) &Complex::complexUnit,
	//             (_MKL_Complex16*) &(*this)(0, 0), &m);
        TensorCalculus::Blas< std::complex<double> >::gemm(notConjTrans, conjTrans, MIN(m,mR), MIN(n,nR), k, COMPLEX_UNIT,
			    &(R.attr_A(rowIndex, 0)), mR, &(R.attr_B(colIndex, 0)), nR, COMPLEX_UNIT,
	            &(*this)(0, 0), m);

     }
    else if(theMatrixA==isCMatrix)
     {
        const CMatrix& C = (const CMatrix&)A;        
        int ic = 1;
        for(LongInt i=0; i<n; i++)
         {
	   // zaxpy (&m,  (_MKL_Complex16*) &Complex::complexUnit, (_MKL_Complex16*) &C(rowIndex, colIndex+i), 
           //        &ic, (_MKL_Complex16*) &(*this)(0, i), &ic);
	   TensorCalculus::Blas< std::complex<double> >::axpy (m,  COMPLEX_UNIT, &C(rowIndex, colIndex+i), ic, &(*this)(0, i), ic);
         }
     }
    else if(theMatrixA==isHMatrixInterface)
     {    
        CMatrix C(m1, n1);

        C.isConversionOf(A);
        addMatrix(C, rowIndex, colIndex);
     }
    else
     {
        throw SimpleException(IString("Warning In CMatrix::addMatrix (const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex), Type Not Implemented !!!"));
     }    

   return (*this);
 }


CMatrix&  CMatrix::addScalMatrix (const std::complex<double>& z, const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex)
 {
    const AMatrixType theMatrixA = A.type();

    LongInt m  = numberOfRows();
    LongInt n  = numberOfColumns();
    LongInt m1 = A.numberOfRows();
    LongInt n1 = A.numberOfColumns();
    
    if(theMatrixA==isRkCMatrix)
     {

        const RkCMatrix& R = (const RkCMatrix&)A;
        LongInt k = R.rank();

        LongInt mR = R.numberOfRows();
        LongInt nR = R.numberOfColumns(); 
      
        // zgemm ( &notConjTrans, &conjTrans, &MIN(m,mR), &MIN(n,nR), &k, (_MKL_Complex16*) &z,
	//             (_MKL_Complex16*) &(R.attr_A(rowIndex, 0)), &mR, 
        //         (_MKL_Complex16*) &(R.attr_B(colIndex, 0)), &nR, (_MKL_Complex16*) &Complex::complexUnit,
	//             (_MKL_Complex16*) &(*this)(0, 0), &m);
        TensorCalculus::Blas< std::complex<double> >::gemm (notConjTrans, conjTrans, MIN(m,mR), MIN(n,nR), k, z,
	            &(R.attr_A(rowIndex, 0)), mR, 
                &(R.attr_B(colIndex, 0)), nR, COMPLEX_UNIT,
	            &(*this)(0, 0), m);

     }
    else if(theMatrixA==isCMatrix)
     {
        const CMatrix& C = (const CMatrix&)A;        
        int ic = 1;
        for(LongInt i=0; i<n; i++)
         {
	   // zaxpy (&m,  (_MKL_Complex16*) &z, (_MKL_Complex16*) &C(rowIndex, colIndex+i), 
           //         &ic, (_MKL_Complex16*) &(*this)(0, i), &ic);
           TensorCalculus::Blas< std::complex<double> >::axpy (m,  z, &C(rowIndex, colIndex+i), ic, &(*this)(0, i), ic);
         }
     }
    else if(theMatrixA==isHMatrixInterface)
     {    
        CMatrix C(m1, n1);

        C.isConversionOf(A);
        addScalMatrix(z, C, rowIndex, colIndex);
     }
    else
     {
        throw SimpleException(IString("Warning In CMatrix::addScalMatrix (const Complex& z, const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex), Type Not Implemented !!!"));
     }    

   return (*this);
 }


CMatrix& CMatrix::operator -= (const RMatrix& A)
 {
    const LongInt n  = attr_numberOfRows;
    const LongInt m  = attr_numberOfColumns;

    for(LongInt i=0; i<n; i++)
     {
        for(LongInt j=0; j<m; j++)
         {
            (*this)(i,j) -= A(i,j);
         }
     }

   return (*this);
 }


CMatrix& CMatrix::operator *= (const std::complex<double>& alpha)
 {
    integer n  = attr_numberOfRows*attr_numberOfColumns;
    integer ic = 1;

    // zscal (&n, (_MKL_Complex16*) &alpha, (_MKL_Complex16*) (&this->_pelm[0]), &ic);
    TensorCalculus::Blas< std::complex<double> >::scal (n, alpha, (&this->_pelm[0]), ic);
   return (*this);
 }


CMatrix& CMatrix::operator /= (const std::complex<double>& alpha)
 {
    (*this) *= 1.0/alpha; 
   return (*this);
 }


CMatrix& CMatrix::update(const std::complex<double>& a, CMatrix& A, CMatrix& B, const std::complex<double>& b,
                         const LongInt& rowOffset, const LongInt& colOffset)
 {
    integer l = A.numberOfRows();
    integer n = B.numberOfColumns();
    integer m = A.numberOfColumns();
    integer lda = l;
    integer ldb = m;
             
    if(A.isConjugateTranspose() || A.isTranspose())
     {
        l   = A.numberOfColumns();
        m   = A.numberOfRows();
        lda = m;
        ldb = m;
     }

    if(B.isConjugateTranspose() || B.isTranspose())
     {
        m   = B.numberOfColumns();
        n   = B.numberOfRows();
        ldb = n;
     }

    integer ldc = l;

    (*this).setDimension(l, n);

    TensorCalculus::Blas< std::complex<double> >::gemm (*A.correctChar(), *B.correctChar(), l, n, m,
           a, &A._pelm[rowOffset], lda, 
           &B._pelm[colOffset], ldb, b, 
           (&this->_pelm[0]), ldc);

    // zgemm (A.correctChar(), B.correctChar(), &l, &n, &m, 
    //        (_MKL_Complex16*) &a, (_MKL_Complex16*) &A._pelm[rowOffset], &lda, 
    //        (_MKL_Complex16*) &B._pelm[colOffset], &ldb, (_MKL_Complex16*) &b, 
    //        (_MKL_Complex16*) (&this->_pelm[0]), &ldc);

      clearConjugateTranspose();
      A.clearConjugateTranspose();
      B.clearConjugateTranspose();

   return (*this);
 }


CMatrix& CMatrix::operator *= (CMatrix& B)
 {
    CMatrix A((*this));

    (*this).update(COMPLEX_UNIT, A, B, COMPLEX_NULL);

   return (*this);
 }


CMatrix& CMatrix::operator *= (const RkCMatrix& R)
 {
    (*this) *= (CMatrix &)R.attr_A;
    (*this) *= (((CMatrix &)R.attr_B).setConjugateTranspose());
   return (*this);
 }

CMatrix& CMatrix::leftMultiplied( CMatrix& A)
 {
    CMatrix B((*this));

    (*this).update(COMPLEX_UNIT, A, B, COMPLEX_NULL);

   return (*this); 
 }


CMatrix&  CMatrix::isConversionOf (const AMatrix& A)
 {    
    const AMatrixType theMatrixA = A.type();
    
    if(theMatrixA==isCMatrix)
     {
        (*this) = (CMatrix&) A;
     }
    else if(theMatrixA==isRkCMatrix)
     {
        (*this) = ((RkCMatrix&) A).convert2CMatrix();
     }
    else 
     {
        const LongInt m = A.numberOfRows();
        const LongInt n = A.numberOfColumns();
        
        CVector temp(n);

        setDimension(m, n);

        for(LongInt j=0; j<n; j++)
         {
            CVector helpVp_w(m, &(*this)(0,j));
            temp(j) = COMPLEX_UNIT;
            A(0, helpVp_w, 0, temp);
            temp(j) = COMPLEX_NULL;
         }           
     }

   return (*this);
 }


CMatrix& CMatrix::update(const AMatrix& A, const AMatrix& B)
 {
    const AMatrixType theMatrixA = A.type();
    const AMatrixType theMatrixB = B.type();

    if(theMatrixA==isCMatrix && theMatrixB==isCMatrix)
     {
        update(COMPLEX_UNIT, (CMatrix &)A , (CMatrix &) B, COMPLEX_UNIT);
     }
    else if (theMatrixA!=isCMatrix && theMatrixB==isCMatrix)
     {
        const LongInt m = A.numberOfRows();
        const LongInt n = B.numberOfColumns();

        CMatrix& C = (CMatrix&)B;

        for(LongInt j=0; j<n; j++)
         {
            CVector helpVp_v(m, &C(0,j)); 
            CVector helpVp_w(m, &(*this)(0,j));
  
            A.addEvaluateAt(helpVp_v, helpVp_w);
         }
     }
    else if (theMatrixA==isCMatrix && theMatrixB!=isCMatrix)
     {
        LongInt m = A.numberOfRows();
        LongInt k = A.numberOfColumns();
        LongInt n = B.numberOfColumns();
        int     e = 1;

        CVector v(k), w(n);
        CMatrix& C = (CMatrix&)A;
        for(int i=0; i<m; i++)            
         {
           // zcopy (&k, (_MKL_Complex16*) &C(i,0), &m, (_MKL_Complex16*) &v(0), &e);
           //  zcopy (&n, (_MKL_Complex16*) &(*this)(i,0), &m, (_MKL_Complex16*) &w(0), &e);
           //  B.addEvaluateTransposeAt(v, w);
           //  zcopy (&n, (_MKL_Complex16*) &w(0), &e, (_MKL_Complex16*) &(*this)(i,0), &m);

           TensorCalculus::Blas< std::complex<double> >::copy (k, &C(i,0), m, &v(0), e);
           TensorCalculus::Blas< std::complex<double> >::copy (n, &(*this)(i,0), m, &w(0), e);
	       B.addEvaluateTransposeAt(v, w);
           TensorCalculus::Blas< std::complex<double> >::copy (n, &w(0), e, &(*this)(i,0), m);
         }
     }
    else if (theMatrixA!=isCMatrix && theMatrixB!=isCMatrix)
     {
        const LongInt m = A.numberOfRows();
        const LongInt n = A.numberOfColumns();
        const LongInt o = B.numberOfColumns();
        
        CMatrix S(m,n);
        S.isConversionOf(A);

        CMatrix T(n,o);
        T.isConversionOf(B);
    
        update(COMPLEX_UNIT, S , T, COMPLEX_UNIT);
     }

   return (*this);
 }


CMatrix CMatrix :: operator + (const CMatrix& A)
 {	
    CMatrix result((*this));
 
    result += A;

   return result;
 }


CMatrix CMatrix :: operator - (const CMatrix& A)
 {	
    CMatrix result((*this));
 
    result -= A;

   return result;
 }


CMatrix CMatrix :: operator * ( CMatrix& A)
 {
    CMatrix result((*this));
 
    result *= A;

   return result;
 }


CMatrix  operator / (const CMatrix& A, const std::complex<double> alpha)
 {
    CMatrix result(A);
 
    result /= alpha;

   return result;
 }


CMatrix  operator * (const std::complex<double> alpha, const CMatrix& A)
 {	
    CMatrix result(A);
 
    result *= alpha;

   return result;
 }


CVector CMatrix :: operator * (const CVector& v)
 {
    integer incx = 1;
    integer m    = attr_numberOfRows;
    integer n    = attr_numberOfColumns;
    integer lda  = attr_numberOfRows;
    
	CVector result (attr_numberOfRows);

	TensorCalculus::Blas< std::complex<double> >::gemv (*correctChar(), m, n, COMPLEX_UNIT, (&_pelm[0]),
            lda, &v(0), incx, COMPLEX_NULL,
            &result(0), incx );

    // zgemv ( correctChar(), &m, &n, (_MKL_Complex16*) &Complex::complexUnit, (_MKL_Complex16*) (&_pelm[0]), 
    //         &lda, (_MKL_Complex16*) &v(0), &incx, (_MKL_Complex16*) &Complex::complexNull, 
    //         (_MKL_Complex16*)&result(0), &incx
    //       );
        clearConjugateTranspose();
        clearTranspose();
   return result;
 }


bool CMatrix::operator () (const LongInt& i, CVector& w, const LongInt& j, const CVector& v)
 {
    bool value = true;

    integer incx = 1;
    integer m    = numberOfRows();
    integer n    = numberOfColumns();
    integer lda  = m;

    TensorCalculus::Blas< std::complex<double> >::gemv (*correctChar(), m, n, COMPLEX_UNIT, (&_pelm[0]),
            lda, &v(j), incx, COMPLEX_UNIT,
            &w(i), incx);

    // zgemv (  correctChar(), &m, &n, (_MKL_Complex16*) &Complex::complexUnit, (_MKL_Complex16*) (&_pelm[0]), 
    //         &lda, (_MKL_Complex16*) &v(j), &incx, (_MKL_Complex16*) &Complex::complexUnit, 
    //         (_MKL_Complex16*)&w(i), &incx
    //       );

    clearConjugateTranspose();
    clearTranspose();
   return value;
 }


bool CMatrix::partialEvaluateTransposeAt (const LongInt& i, CVector& w, const LongInt& j, const CVector& v) const
{
    bool value = true;

    integer incx = 1;
    integer n    = numberOfRows();
    integer m    = numberOfColumns();
    integer lda  = n;
    

    // zgemv (  &trans, &n, &m, (_MKL_Complex16*) &Complex::complexUnit, (_MKL_Complex16*) (&_pelm[0]), 
    //         &lda, (_MKL_Complex16*) &v(j), &incx, (_MKL_Complex16*) &Complex::complexUnit, 
    //         (_MKL_Complex16*)&w(i), &incx
    //       );
    TensorCalculus::Blas< std::complex<double> >::gemv (trans, n, m, COMPLEX_UNIT, (&_pelm[0]),
            lda, &v(j), incx, COMPLEX_UNIT,
            &w(i), incx
          );

   return value;
 }
 

bool CMatrix::partialEvaluateConjugateTransposeAt (const LongInt& i, CVector& w, const LongInt& j, const CVector& v) const
{
    bool value = true;

    integer incx = 1;
    integer n    = numberOfRows();
    integer m    = numberOfColumns();
    integer lda  = n;
    

    // zgemv (  &conjTrans, &n, &m, (_MKL_Complex16*) &Complex::complexUnit, (_MKL_Complex16*) (&_pelm[0]), 
    //         &lda, (_MKL_Complex16*) &v(j), &incx, (_MKL_Complex16*) &Complex::complexUnit, 
    //         (_MKL_Complex16*)&w(i), &incx
    //       );
    TensorCalculus::Blas< std::complex<double> >::gemv (conjTrans, n, m, COMPLEX_UNIT,  (&_pelm[0]),
            lda, &v(j), incx, COMPLEX_UNIT,
            &w(i), incx
          );

   return value;
 }


bool CMatrix::operator == (const CMatrix& A) const
 {

	for( LongInt i=0 ; i < attr_numberOfRows ; i++)
     {
        for(LongInt k=0; i<attr_numberOfColumns; k++)
         {
		    if((*this)(i,k) != A(i,k))
             {
               return false;
             }
         }
	 }
	return true;
 }



CMatrix&  CMatrix::invert()
 {
    integer m   = numberOfRows();
    integer n   = numberOfColumns();
    integer lda = MAX(1,m);
    integer dim = MIN(m,n);
    integer info = 0; 

    integer* ipiv = new integer[dim];

    // zgetrf (&m, &n, (_MKL_Complex16*) &_pelm[0], &lda, ipiv, &info);
    info = TensorCalculus::Lapack< std::complex<double> >::getrf (m, n, &_pelm[0], lda, ipiv);

    integer lwork = m*n;
    CVector work(lwork);

    // zgetri (&m, (_MKL_Complex16*) &_pelm[0], &lda, ipiv, (_MKL_Complex16*) &work(0), &lwork, &info);
    info = TensorCalculus::Lapack< std::complex<double> >::getri (m, &_pelm[0], lda, ipiv, &work(0), lwork);

    delete [] ipiv;
    
   return (*this);
 }


CMatrix& CMatrix::isInverseOf(const AMatrix& s, AMatrix& w)
 {
    const AMatrixType theMatrixS = s.type();
    const AMatrixType theMatrixW = w.type();
    
    if(theMatrixS==isCMatrix && theMatrixW==isCMatrix)
     {
        const CMatrix& S = (const CMatrix&)s;
        const CMatrix& W = (const CMatrix&)w;        

        (*this) = S; 

        integer m   = numberOfRows();
        integer n   = numberOfColumns();
        integer lda = MAX(1,m);
        integer dim = MIN(m,n);
        integer info = 0; 
        integer lwork = m*n;

        integer* ipiv = new integer[dim];

        info = TensorCalculus::Lapack< std::complex<double> >::getrf(m, n, &_pelm[0], lda, ipiv);
	    info = TensorCalculus::Lapack< std::complex<double> >::getri(m, &_pelm[0], lda, ipiv, &W(0,0), lwork);
	// zgetrf (&m, &n, (_MKL_Complex16*) &_pelm[0], &lda, ipiv, &info);
	// zgetri (&m, (_MKL_Complex16*) &_pelm[0], &lda, ipiv, (_MKL_Complex16*) &W(0,0), &lwork, &info);

        delete [] ipiv;
     }
    else
     {
        throw SimpleException(IString("Warning In CMatrix::isInverseOf(const AMatrix& s, const AMatrix& w), Bad Type!!!"));     
     }

   return (*this);
 }


char* CMatrix::correctChar() const
 {
    if(isConjugateTranspose())
     {
       return &conjTrans;
     }
    else if(isTranspose())
     {
       return &trans;
     }
    else
     {
       return &notConjTrans;
     }
 }


CMatrix& CMatrix::attachColumns (const CMatrix& A)
 {
    const LongInt m  = numberOfRows();
    const LongInt n1 = numberOfColumns();
    const LongInt n2 = A.numberOfColumns();
    const LongInt n  = n1+n2;

    CMatrix temp((*this));
  
    setDimension(m, n);

    integer dim1 = m*n1;
    integer nx   = 1;

    TensorCalculus::Blas< std::complex<double> >::copy (dim1, (&temp._pelm[0]), nx, (&_pelm[0]), nx);
    // zcopy (&dim1, (_MKL_Complex16*)(&temp._pelm[0]), &nx, (_MKL_Complex16*)(&_pelm[0]), &nx);

    integer dim2 = m*n2;

    TensorCalculus::Blas< std::complex<double> >::copy (dim2, (&A._pelm[0]), nx, (&_pelm[dim1]), nx);
    // zcopy (&dim2, (_MKL_Complex16*)(&A._pelm[0]), &nx, (_MKL_Complex16*)(&_pelm[dim1]), &nx);

   return (*this);
 }


CMatrix&  CMatrix::attachRows (const CMatrix& A)
 {
    //! \todo zcopy_ verwenden
    const LongInt n  = numberOfColumns();
    const LongInt m1 = numberOfRows();
    const LongInt m2 = A.numberOfRows();
    const LongInt m  = m1+m2; 

    CMatrix temp((*this));

    setDimension(m,n);
     
    for(LongInt j=0; j<n; j++)
     {  
        LongInt i;
        for(i=0; i<m1; i++)
         {
            (*this)(i,j) = temp(i,j);
         }
        for(i=0; i<m2; i++)
         {
            (*this)(i+m1,j) = A(i,j);
         }
     }

   return (*this);
 }


RMatrix CMatrix::realPart() const
 {
    const LongInt n = attr_numberOfRows;
    const LongInt m = attr_numberOfColumns;
    
    RMatrix rePart(n,m);

    for(LongInt i=0; i<n; i++)
     {
        for(LongInt j=0; j<m; j++)
         {
            rePart(i,j) = (*this)(i,j).real(); 
         }
     }
   return rePart; 
 }


RMatrix CMatrix::imagPart() const
 {
    const LongInt n = attr_numberOfRows;
    const LongInt m = attr_numberOfColumns;
    
    RMatrix imPart(n,m);

    for(LongInt i=0; i<n; i++)
     {
        for(LongInt j=0; j<m; j++)
         {
            imPart(i,j) = (*this)(i,j).imag(); 
         }
     }
   return imPart; 
 }


ostream&  operator << (ostream & s,const CMatrix &A)
 {
    LongInt m = A.numberOfRows();
    LongInt n = A.numberOfColumns();

	s << m << '\t'<< n <<'\n';

	for (LongInt i=0; i<m; i++)
     {
	    for(LongInt k=0; k<n; k++)
         {
		    s << A(i,k)<<'\t';
		 }
		s << '\n';
	 }	

   return s; 
 }


istream&  operator >> (istream & s, CMatrix &A)
 {
    LongInt numberOfRows=0, numberOfColumns=0;
	s >> numberOfRows;
	s >> numberOfColumns;
	A.setDimension(numberOfRows, numberOfColumns);

	for (LongInt i=0; i < numberOfRows ; i++)
     {
	    for(LongInt k=0; k < numberOfColumns;k++)
         {
		    s >> A(i,k);
		 }
	}

   return s;
 }


CMatrix& CMatrix :: setDimension(LongInt n, LongInt m)
 {
    if(n!=attr_numberOfRows || m!=attr_numberOfColumns)
     {
        delete []_pelm;	

	    attr_numberOfRows    = n;
	    attr_numberOfColumns = m;

        _pelm =  (std::complex<double>*) new std::complex<double> [attr_numberOfRows*attr_numberOfColumns];
     }

   return (*this);
 }


CMatrix& CMatrix::conjugateTransposed()
 {
    LongInt i, j;

    const LongInt m = numberOfRows();
    const LongInt n = numberOfColumns();

    CMatrix value((*this));

    setDimension(n,m);

	for (i=0; i<m ; i++)
     {
	    for(j=0; j<n ; j++)
         {
            (*this)(j,i) = std::conj(value(i,j));
         }
     }
    
   return (*this);
 }


CMatrix&  CMatrix::conjugated ()
 {
    LongInt i;

    const LongInt dim = numberOfRows()*numberOfColumns();

	for (i=0; i<dim ; i++)
     {
        _pelm[i] = std::conj(_pelm[i]); // TODO: faster
     }

   return (*this);
 }


CMatrix& CMatrix :: setColumn(LongInt k, const CVector &v)
 {
    LongInt i;
	for (i=0 ; i < attr_numberOfRows ; i++)
		(*this)(i,k)=v(i);
   return (*this);
 }


CVector CMatrix :: getColumn(LongInt k)
 {
    LongInt i;
	CVector temp(attr_numberOfRows);
	for (i=0 ; i < attr_numberOfRows ; i++)
		temp(i)=(*this)(i,k);
   return temp;
 }


CMatrix& CMatrix :: setRow(LongInt i, const CVector &v)
 {
	for (LongInt k=0 ; k < attr_numberOfColumns ; k++)
		(*this)(i,k)=v(k);
   return (*this);
 }


CVector CMatrix :: getRow(LongInt i)
 {
	CVector temp(attr_numberOfColumns);
	for (LongInt k=0 ; k < attr_numberOfColumns ; k++)
		temp(k)=(*this)(i,k);
   return temp;
 }


AMatrix& CMatrix::setNullMatrix()
 {
    (*this) *= COMPLEX_NULL;
   return (*this);
 }


LongReal CMatrix::frobeniusNorm () const
 {
    LongReal temp = 0.0;

    integer dim  = numberOfRows()*numberOfColumns();
    integer ic = 1;
 
    // temp = dznrm2 (&dim, (_MKL_Complex16*) &(*this)(0,0), &ic);
    temp = TensorCalculus::Blas< std::complex<double> >::nrm2 (dim, &(*this)(0,0), ic);

   return temp;
 }
