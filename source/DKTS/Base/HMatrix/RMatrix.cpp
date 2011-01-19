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

// RMatrix.cpp: Implementierung der Klasse RMatrix.
//
//////////////////////////////////////////////////////////////////////

#include "RMatrix.hpp"
#include "RkRMatrix.hpp"
#include "BlasInterface.hpp"
#include "LapackInterface2.hpp"

RMatrix::RMatrix(LongInt m, LongInt n)
:ARMatrix(m, n)
 {
    const LongInt dim = m*n;

    (*this)
     .disableTranspose()     
     .setType(isRMatrix)
    ;

    _pelm =  (LongRealPointer) new LongReal [dim];

    LongRealPointer p = &_pelm[0];

    for(LongInt i=0; i<dim; i++)
     {
        *p = 0.0;
        p++;
     }

 }


RMatrix::~RMatrix()
 {
	   delete [] _pelm;
 }


RMatrix::RMatrix(const RMatrix &A)
:ARMatrix(A.numberOfRows(), A.numberOfColumns())
 {   
    _pelm = (LongRealPointer) new LongReal [attr_numberOfRows*attr_numberOfColumns];

	(*this) = A;
 }


RMatrix& RMatrix :: operator = (const RMatrix& A)
 {	
    LongInt n = A.numberOfRows();
    LongInt m = A.numberOfColumns();   
    
    setDimension(n, m);
 
    attr_type  =  A.attr_type;
    setTranspose(A.isTranspose());
    

    integer dim = n*m;
    integer nx  = 1;

    TensorCalculus::Blas<double>::copy (dim, (&A._pelm[0]), nx, (&this->_pelm[0]), nx);
    
   return (*this);
 }


RMatrix& RMatrix::update (const LongReal& a, const RMatrix& B)
 {
    integer n  = attr_numberOfRows*attr_numberOfColumns;
    integer ic = 1;

    // daxpy (&n, (LongReal*)&a, (LongReal*) &B._pelm[0], &ic, (LongReal*) &(this->_pelm[0]), &ic);
    TensorCalculus::Blas<double>::axpy (n, a, &B._pelm[0], ic, &(this->_pelm[0]), ic);

   return (*this);
 }


RMatrix& RMatrix::operator += (const RMatrix& A)
 {
    update(COMPLEX_UNIT.real(), A);

   return (*this);
 }


RMatrix& RMatrix::setTranspose(bool b)
 {
    AMatrix::setTranspose(b);

   return (*this);
 }


RMatrix& RMatrix::operator -= (const RMatrix& A)
 {
    update(-COMPLEX_UNIT.real(), A);

   return (*this);
 }


RMatrix&  RMatrix::addIdentityScaledBy (const LongReal& z)
 {
    const LongInt dim = MIN(numberOfRows(),numberOfColumns());
  
    for(LongInt i=0; i<dim; i++)
     {
        (*this)(i,i) += z;
     }

   return(*this);
 }


RMatrix& RMatrix::addMatrix (const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex)
 {
    const AMatrixType theMatrixA = A.type();

    LongInt m  = numberOfRows();
    LongInt n  = numberOfColumns();
    LongInt m1 = A.numberOfRows();
    LongInt n1 = A.numberOfColumns();
    
    if(theMatrixA==isRkRMatrix)
     {

        const RkRMatrix& R = (const RkRMatrix&)A;
        LongInt k = R.rank();

        LongInt mR = R.numberOfRows();
        LongInt nR = R.numberOfColumns(); 

        LongReal mike = 1.0;
      
	//  dgemm ( &notConjTrans, &conjTrans, &MIN(m,mR), &MIN(n,nR), &k, (LongReal*) &(Complex::complexUnit),
	//            (LongReal*) &(R.attr_A(rowIndex, 0)), &mR, 
        //         (LongReal*) &(R.attr_B(colIndex, 0)), &nR, (LongReal*) &(Complex::complexUnit),
	//             (LongReal*) &(*this)(0, 0), &m);
	TensorCalculus::Blas<double>::gemm (notConjTrans, conjTrans, MIN(m,mR), MIN(n,nR), k, 1.0,
	            &(R.attr_A(rowIndex, 0)), mR, 
                &(R.attr_B(colIndex, 0)), nR, 1.0,
	            &(*this)(0, 0), m);

     }
    else if(theMatrixA==isRMatrix)
     {
        const RMatrix& C = (const RMatrix&)A;        
        int ic = 1;
        for(LongInt i=0; i<n; i++)
         {
	   // daxpy (&m,  (LongReal*) &(Complex::complexUnit), (LongReal*) &C(rowIndex, colIndex+i), 
           //         &ic, (LongReal*) &(*this)(0, i), &ic);
           TensorCalculus::Blas<double>::axpy (m, 1.0, &C(rowIndex, colIndex+i), ic, &(*this)(0, i), ic);
          }
     }
    else if(theMatrixA==isHMatrixInterface)
     {    
        RMatrix C(m1, n1);

        C.isConversionOf(A);
        addMatrix(C, rowIndex, colIndex);
     }
    else
     {
        throw SimpleException(IString("Warning In RMatrix::addMatrix (const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex), Type Not Implemented !!!"));
     }    

   return (*this);
 }


RMatrix&  RMatrix::addScalMatrix (const LongReal& z, const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex)
 {
    const AMatrixType theMatrixA = A.type();

    LongInt m  = numberOfRows();
    LongInt n  = numberOfColumns();
    LongInt m1 = A.numberOfRows();
    LongInt n1 = A.numberOfColumns();
    
    if(theMatrixA==isRkRMatrix)
     {

        const RkRMatrix& R = (const RkRMatrix&)A;
        LongInt k = R.rank();

        LongInt mR = R.numberOfRows();
        LongInt nR = R.numberOfColumns(); 
      
        // dgemm ( &notConjTrans, &conjTrans, &MIN(m,mR), &MIN(n,nR), &k, (LongReal*) &z,
	//             (LongReal*) &(R.attr_A(rowIndex, 0)), &mR, 
        //         (LongReal*) &(R.attr_B(colIndex, 0)), &nR, (LongReal*) &(Complex::complexUnit),
	//             (LongReal*) &(*this)(0, 0), &m);
        TensorCalculus::Blas<double>::gemm (notConjTrans, conjTrans, MIN(m,mR), MIN(n,nR), k, z,
	            &(R.attr_A(rowIndex, 0)), mR, &(R.attr_B(colIndex, 0)), nR, 1.0,
	            &(*this)(0, 0), m);

     }
    else if(theMatrixA==isRMatrix)
     {
        const RMatrix& C = (const RMatrix&)A;        
        int ic = 1;
        for(LongInt i=0; i<n; i++)
         {
	   // daxpy (&m,  (LongReal*) &z, (LongReal*) &C(rowIndex, colIndex+i), 
           //         &ic, (LongReal*) &(*this)(0, i), &ic);
	   TensorCalculus::Blas<double>::axpy (m,  z, &C(rowIndex, colIndex+i), ic, &(*this)(0, i), ic);
         }
     }
    else if(theMatrixA==isHMatrixInterface)
     {    
        RMatrix C(m1, n1);

        C.isConversionOf(A);
        addScalMatrix(z, C, rowIndex, colIndex);
     }
    else
     {
        throw SimpleException(IString("Warning In RMatrix::addScalMatrix (const LongReal& z, const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex), Type Not Implemented !!!"));
     }    

   return (*this);
 }



RMatrix& RMatrix::operator *= (const LongReal& alpha)
 {
    integer n  = attr_numberOfRows*attr_numberOfColumns;
    integer ic = 1;

    // dscal (&n, (LongReal*) &alpha, (LongReal*) (&this->_pelm[0]), &ic);
    TensorCalculus::Blas<double>::scal (n, alpha, (&this->_pelm[0]), ic);
   return (*this);
 }


RMatrix& RMatrix::operator /= (const LongReal& alpha)
 {
    (*this) *= 1.0/alpha; 
   return (*this);
 }


RMatrix& RMatrix::update(const LongReal& a, RMatrix& A, RMatrix& B, const LongReal& b,
                         const LongInt& rowOffset, const LongInt& colOffset)
 {
    integer l = A.numberOfRows();
    integer n = B.numberOfColumns();
    integer m = A.numberOfColumns();
    integer lda = l;
    integer ldb = m;
             
    if(A.isTranspose() || A.isTranspose())
     {
        l   = A.numberOfColumns();
        m   = A.numberOfRows();
        lda = m;
        ldb = m;
     }

    if(B.isTranspose() || B.isTranspose())
     {
        m   = B.numberOfColumns();
        n   = B.numberOfRows();
        ldb = n;
     }

    integer ldc = l;

    (*this).setDimension(l, n);

    // dgemm (A.correctChar(), B.correctChar(), &l, &n, &m, 
    //        (LongReal*) &a, (LongReal*) &A._pelm[rowOffset], &lda, 
    //        (LongReal*) &B._pelm[colOffset], &ldb, (LongReal*) &b, 
    //        (LongReal*) (&this->_pelm[0]), &ldc);
    TensorCalculus::Blas<double>::gemm (*A.correctChar(), *B.correctChar(), l, n, m, 
           a, &A._pelm[rowOffset], lda, 
           &B._pelm[colOffset], ldb, b, 
           (&this->_pelm[0]), ldc);

      clearTranspose();
      A.clearTranspose();
      B.clearTranspose();

   return (*this);
 }


RMatrix& RMatrix::operator *= (RMatrix& B)
 {
    RMatrix A((*this));

    (*this).update(COMPLEX_UNIT.real(), A, B, COMPLEX_NULL.real());

   return (*this);
 }


RMatrix& RMatrix::operator *= (const RkRMatrix& R)
 {
    (*this) *= (RMatrix &)R.attr_A;
    (*this) *= (((RMatrix &)R.attr_B).setTranspose());
   return (*this);
 }

RMatrix& RMatrix::leftMultiplied( RMatrix& A)
 {
    RMatrix B((*this));

    (*this).update(COMPLEX_UNIT.real(), A, B, COMPLEX_NULL.real());

   return (*this); 
 }


RMatrix&  RMatrix::isConversionOf (const AMatrix& A)
 {    
    const AMatrixType theMatrixA = A.type();
    
    if(theMatrixA==isRMatrix)
     {
        (*this) = (RMatrix&) A;
     }
    else if(theMatrixA==isRkRMatrix)
     {
        (*this) = ((RkRMatrix&) A).convert2RMatrix();
     }
    else 
     {
        const LongInt m = A.numberOfRows();
        const LongInt n = A.numberOfColumns();
        
        RVector temp(n);

        setDimension(m, n);

        for(LongInt j=0; j<n; j++)
         {
            RVector helpVp_w(m, &(*this)(0,j));
            temp(j) = 1.0;
            A(0, helpVp_w, 0, temp);
            temp(j) = 0.0;
         }           
     }

   return (*this);
 }


RMatrix& RMatrix::update(const AMatrix& A, const AMatrix& B)
 {
    const AMatrixType theMatrixA = A.type();
    const AMatrixType theMatrixB = B.type();

    if(theMatrixA==isRMatrix && theMatrixB==isRMatrix)
     {
        update(COMPLEX_UNIT.real(), (RMatrix &)A , (RMatrix &) B, COMPLEX_UNIT.real());
     }
    else if (theMatrixA!=isRMatrix && theMatrixB==isRMatrix)
     {
        const LongInt m = A.numberOfRows();
        const LongInt n = B.numberOfColumns();

        RMatrix& C = (RMatrix&)B;

        for(LongInt j=0; j<n; j++)
         {
            RVector helpVp_v(m, &C(0,j)); 
            RVector helpVp_w(m, &(*this)(0,j));
  
            A.addEvaluateAt(helpVp_v, helpVp_w);
         }
     }
    else if (theMatrixA==isRMatrix && theMatrixB!=isRMatrix)
     {
        LongInt m = A.numberOfRows();
        LongInt k = A.numberOfColumns();
        LongInt n = B.numberOfColumns();
        int     e = 1;

        RVector v(k), w(n);
        RMatrix& C = (RMatrix&)A;
        for(int i=0; i<m; i++)            
         {
	   TensorCalculus::Blas<double>::copy (k, &C(i,0), m, &v(0), e);
           TensorCalculus::Blas<double>::copy (n, &(*this)(i,0), m, &w(0), e);
	   B.addEvaluateTransposeAt(v, w);
	   TensorCalculus::Blas<double>::copy (n, &w(0), e, &(*this)(i,0), m);
            // dcopy (&k, (LongReal*) &C(i,0), &m, (LongReal*) &v(0), &e);
            // dcopy (&n, (LongReal*) &(*this)(i,0), &m, (LongReal*) &w(0), &e);
            // B.addEvaluateTransposeAt(v, w);
            // dcopy (&n, (LongReal*) &w(0), &e, (LongReal*) &(*this)(i,0), &m);
         }
     }
    else if (theMatrixA!=isRMatrix && theMatrixB!=isRMatrix)
     {
        const LongInt m = A.numberOfRows();
        const LongInt n = A.numberOfColumns();
        const LongInt o = B.numberOfColumns();
        
        RMatrix S(m,n);
        S.isConversionOf(A);

        RMatrix T(n,o);
        T.isConversionOf(B);
    
        update(COMPLEX_UNIT.real(), S , T, COMPLEX_UNIT.real());
     }

   return (*this);
 }


RMatrix RMatrix :: operator + (const RMatrix& A)
 {	
    RMatrix result((*this));
 
    result += A;

   return result;
 }


RMatrix RMatrix :: operator - (const RMatrix& A)
 {	
    RMatrix result((*this));
 
    result -= A;

   return result;
 }


RMatrix RMatrix :: operator * ( RMatrix& A)
 {
    RMatrix result((*this));
 
    result *= A;

   return result;
 }


RMatrix  operator / (const RMatrix& A, const LongReal alpha)
 {
    RMatrix result(A);
 
    result /= alpha;

   return result;
 }


RMatrix  operator * (const LongReal alpha, const RMatrix& A)
 {	
    RMatrix result(A);
 
    result *= alpha;

   return result;
 }


RVector RMatrix :: operator * (const RVector& v)
 {
    integer incx = 1;
    integer m    = attr_numberOfRows;
    integer n    = attr_numberOfColumns;
    integer lda  = attr_numberOfRows;
    
	RVector result (attr_numberOfRows);

	// dgemv ( correctChar(), &m, &n, (LongReal*) &(Complex::complexUnit), (LongReal*) (&_pelm[0]), 
        //     &lda, (LongReal*) &v(0), &incx, (LongReal*) &(Complex::complexNull), 
        //     (LongReal*)&result(0), &incx
        //   );
	TensorCalculus::Blas<double>::gemv ( *correctChar(), m, n, 1.0, (&_pelm[0]), 
            lda, &v(0), incx, 0.0, &result(0), incx);

        clearTranspose();
   return result;
 }


bool RMatrix::operator () (const LongInt& i, RVector& w, const LongInt& j, const RVector& v)
 {
    bool value = true;

    integer incx = 1;
    integer m    = numberOfRows();
    integer n    = numberOfColumns();
    integer lda  = m;

    // dgemv (  correctChar(), &m, &n, (LongReal*) &(Complex::complexUnit), (LongReal*) (&_pelm[0]), 
    //         &lda, (LongReal*) &v(j), &incx, (LongReal*) &(Complex::complexUnit), 
    //         (LongReal*) &w(i), &incx
    //       );
    TensorCalculus::Blas<double>::gemv(*correctChar(), m, n, 1.0, &_pelm[0], 
            lda, &v(j), incx, 1.0, &w(i), incx);

    clearTranspose();
   return value;
 }


bool RMatrix::partialEvaluateTransposeAt (const LongInt& i, RVector& w, const LongInt& j, const RVector& v) const
{
    bool value = true;

    integer incx = 1;
    integer n    = numberOfRows();
    integer m    = numberOfColumns();
    integer lda  = n;
    

    // dgemv (  &trans, &n, &m, (LongReal*) &(Complex::complexUnit), (LongReal*) (&_pelm[0]), 
    //          &lda, (LongReal*) &v(j), &incx, (LongReal*) &(Complex::complexUnit), 
    //         (LongReal*)&w(i), &incx
    //       );
    TensorCalculus::Blas<double>::gemv(trans, n, m, 1.0, (&_pelm[0]), lda, &v(j), incx, 1.0, &w(i), incx);

   return value;
 }
 

bool RMatrix::operator == (const RMatrix& A) const
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



RMatrix&  RMatrix::invert()
 {
    integer m   = numberOfRows();
    integer n   = numberOfColumns();
    integer lda = MAX(1,m);
    integer dim = MIN(m,n);
    integer info = 0; 

    integer* ipiv = new integer[dim];

    // dgetrf (&m, &n, (LongReal*) &_pelm[0], &lda, ipiv, &info);
    info = TensorCalculus::Lapack<double>::getrf(m, n, &_pelm[0], lda, ipiv);

    integer lwork = m*n;
    RVector work(lwork);

    // dgetri (&m, (LongReal*) &_pelm[0], &lda, ipiv, (LongReal*) &work(0), &lwork, &info);
    info = TensorCalculus::Lapack<double>::getri(m, &_pelm[0], lda, ipiv, &work(0), lwork);

    delete [] ipiv;
    
   return (*this);
 }


RMatrix& RMatrix::isInverseOf(const AMatrix& s, AMatrix& w)
 {
    const AMatrixType theMatrixS = s.type();
    const AMatrixType theMatrixW = w.type();
    
    if(theMatrixS==isRMatrix && theMatrixW==isRMatrix)
     {
        const RMatrix& S = (const RMatrix&)s;
        const RMatrix& W = (const RMatrix&)w;        

        (*this) = S; 

        integer m   = numberOfRows();
        integer n   = numberOfColumns();
        integer lda = MAX(1,m);
        integer dim = MIN(m,n);
        integer info = 0; 
        integer lwork = m*n;

        integer* ipiv = new integer[dim];

        // dgetrf (&m, &n, (LongReal*) &_pelm[0], &lda, ipiv, &info);
        // dgetri (&m, (LongReal*) &_pelm[0], &lda, ipiv, (LongReal*) &W(0,0), &lwork, &info);
        info = TensorCalculus::Lapack<double>::getrf (m, n, &_pelm[0], lda, ipiv);
        info = TensorCalculus::Lapack<double>::getri (m, &_pelm[0], lda, ipiv, &W(0,0), lwork);

        delete [] ipiv;
     }
    else
     {
        throw SimpleException(IString("Warning In RMatrix::isInverseOf(const AMatrix& s, const AMatrix& w), Bad Type!!!"));     
     }

   return (*this);
 }


RMatrix& RMatrix::isReorganized(const RMatrix& A)
 {    
    const LongInt n = A.numberOfColumns();

    if(true) //!todo Pr�fe ob n eine Quadratzahl ist!
     {
       const LongInt m = (LongInt) (sqrt((LongReal)n)+0.5);

       setDimension(n, n);
    
       LongInt i=0, j=0, iA=0, jA=0;

       for(LongInt i1=0; i1<m; i1++)
        {
           for(LongInt i2=0; i2<m; i2++)
            {
               for(LongInt j1=0; j1<m; j1++)
                {
                   for(LongInt j2=0; j2<m; j2++)
                    {
                       i  = i1+j1*m;
                       j  = i2+j2*m;

                       iA = i1+i2*m;
                       jA = j1+j2*m;
     
                       (*this)(i,j) = A(iA,jA);
                    }
                }
            }
        }
     }
    else
     {
        throw SimpleException(IString("Warning In RMatrix::isReorganized(const RMatrix& A), rows!=columns!!!"));
     }

   return (*this);
 }


char* RMatrix::correctChar() const
 {
    if(isTranspose())
     {
       return &trans;
     }
    else
     {
       return &notConjTrans;
     }
 }


RMatrix& RMatrix::attachColumns (const RMatrix& A)
 {
    const LongInt m  = numberOfRows();
    const LongInt n1 = numberOfColumns();
    const LongInt n2 = A.numberOfColumns();
    const LongInt n  = n1+n2;

    RMatrix temp((*this));
  
    setDimension(m, n);

    integer dim1 = m*n1;
    integer nx   = 1;

    // dcopy (&dim1, (LongReal*)(&temp._pelm[0]), &nx, (LongReal*)(&_pelm[0]), &nx);
    TensorCalculus::Blas<double>::copy (dim1, (&temp._pelm[0]), nx, (&_pelm[0]), nx);

    integer dim2 = m*n2;

    // dcopy (&dim2, (LongReal*)(&A._pelm[0]), &nx, (LongReal*)(&_pelm[dim1]), &nx);
    TensorCalculus::Blas<double>::copy (dim2, (&A._pelm[0]), nx, (&_pelm[dim1]), nx);

   return (*this);
 }


RMatrix&  RMatrix::attachRows (const RMatrix& A)
 {
    //! \todo zcopy_ verwenden
    const LongInt n  = numberOfColumns();
    const LongInt m1 = numberOfRows();
    const LongInt m2 = A.numberOfRows();
    const LongInt m  = m1+m2; 

    RMatrix temp((*this));

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


RMatrix& RMatrix::dyadicProduct (const LongReal& alpha, const RVector& x, const RVector& y)
 {
    LongInt m   = x.dimension();
    LongInt n   = y.dimension();
    LongInt inc = 1;

    //setDimension(m, n);

    // DGER(&m, &n, (LongReal*)&alpha, &x(0), &inc, &y(0), &inc, _pelm, &m);
    TensorCalculus::Blas<double>::ger(m, n, alpha, &x(0), inc, &y(0), inc, _pelm, m);

   return (*this);
 }


ostream&  operator << (ostream & s,const RMatrix &A)
 {
    LongInt m = A.numberOfRows();
    LongInt n = A.numberOfColumns();

	s << m << '\t'<< n <<'\n';

	for (LongInt i=0; i<m; i++)
     {
	    for(LongInt k=0; k<n; k++)
         {
		    s << A(i,k) << '\t';
		 }
		s << '\n';
	 }	

   return s; 
 }


istream&  operator >> (istream & s, RMatrix &A)
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


RMatrix& RMatrix::setId(const LongReal& alpha)
 {
    const LongInt m = MIN(numberOfColumns(), numberOfRows());

    setNullMatrix();

    for(LongInt i=0; i<m; i++)
     {
        (*this)(i,i) = alpha;
     }

   return (*this);
 }


RMatrix& RMatrix::setIdEps(const LongReal& alpha, const LongReal& eps)
 {
    setId(alpha);
    (*this)(0,0) = eps;

   return (*this);
 }


RMatrix& RMatrix::setDiag(const RMatrix& A)
 {
    const LongInt m = A.numberOfRows();
    const LongInt n = A.numberOfColumns();
    const LongInt d = MIN(m,n);
        
    (*this)
     .setDimension(d, d)
    ;
    
    for(LongInt i=0; i<d; i++)
     {
        (*this)(i,i) = A(i,i);
     }    
 
   return (*this);
 }


RMatrix& RMatrix::setRand(const LongReal& eps)
 {
    const LongInt m = numberOfColumns();
    const LongInt n = numberOfRows();
  
    srand( (unsigned)time( NULL ) );

    for(LongInt i=0; i<m; i++)
     {
        for(LongInt j=0; j<n; j++)
         {
            (*this)(i,j) = (LongReal)rand()*eps; 
         }
     }

   return (*this);
 }


RMatrix& RMatrix::setSymRand(const LongReal& eps)
 {
    const LongInt m = numberOfColumns();
    const LongInt n = numberOfRows();
  
    srand( (unsigned)time( NULL ) );

    for(LongInt i=0; i<m; i++)
     {
        for(LongInt j=0; j<=i; j++)
         {
            const LongReal value = (LongReal)rand()*eps;

            (*this)(i,j) = value;
            (*this)(j,i) = value;
         }
     }

   return (*this);
 }


LongReal RMatrix::cond2()const
 {
    RMatrix A((*this));

    LongReal value = norm2(A);

    A.invert();
  
    value *= norm2(A);

   return value;
 }


RMatrix& RMatrix :: setDimension(LongInt n, LongInt m)
 {
    if(n!=attr_numberOfRows || m!=attr_numberOfColumns)
     {
        delete []_pelm;	

	    attr_numberOfRows    = n;
	    attr_numberOfColumns = m;

        _pelm =  (LongRealPointer) new LongReal [attr_numberOfRows*attr_numberOfColumns];

        const LongInt dim = m*n;

        LongRealPointer p = &_pelm[0];

        for(LongInt i=0; i<dim; i++)
         {
            *p = 0.0;
             p++;
         }
     }
    else
     {
        setNullMatrix();
     }

   return (*this);
 }


RMatrix& RMatrix::transposed()
 {
    LongInt i, j;

    const LongInt m = numberOfRows();
    const LongInt n = numberOfColumns();

    RMatrix value((*this));

    setDimension(n,m);

	for (i=0; i<m ; i++)
     {
	    for(j=0; j<n ; j++)
         {
            (*this)(j,i) = value(i,j);
         }
     }
    
   return (*this);
 }


RMatrix&  RMatrix::conjugated ()
 {

   return (*this);
 }


RMatrix& RMatrix :: setColumn(LongInt k, const RVector &v)
 {
    LongInt i;
	for (i=0 ; i < attr_numberOfRows ; i++)
		(*this)(i,k)=v(i);
   return (*this);
 }


RVector RMatrix :: getColumn(LongInt k)
 {
    LongInt i;
	RVector temp(attr_numberOfRows);
	for (i=0 ; i < attr_numberOfRows ; i++)
		temp(i)=(*this)(i,k);
   return temp;
 }


RMatrix& RMatrix :: setRow(LongInt i, const RVector &v)
 {
	for (LongInt k=0 ; k < attr_numberOfColumns ; k++)
		(*this)(i,k)=v(k);
   return (*this);
 }


RVector RMatrix :: getRow(LongInt i)
 {
	RVector temp(attr_numberOfColumns);
	for (LongInt k=0 ; k < attr_numberOfColumns ; k++)
		temp(k)=(*this)(i,k);
   return temp;
 }


AMatrix& RMatrix::setNullMatrix()
 {
    const LongInt m = numberOfColumns();
    const LongInt n = numberOfRows();  

    for(LongInt i=0; i<m; i++)
     {
        for(LongInt j=0; j<n; j++)
         {
            (*this)(i,j) = 0.0; 
         }
     }
    
   return (*this);
 }


LongReal RMatrix::frobeniusNorm () const
 {
    LongReal temp = 0.0;

    integer dim  = numberOfRows()*numberOfColumns();
    integer ic = 1;
 
    // temp = dnrm2 (&dim, (LongReal*) &(*this)(0,0), &ic);
    temp = TensorCalculus::Blas<double>::nrm2 (dim, &(*this)(0,0), ic);

   return temp;
 }
