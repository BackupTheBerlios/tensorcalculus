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

// AMatrix.cpp: Implementierung der Klasse AMatrix.
//
//////////////////////////////////////////////////////////////////////

#include "AMatrix.hpp"
#include "CMatrix.hpp"
#include "RkCMatrix.hpp"
#include "HMatrixInterface.hpp"
#include <time.h>

char AMatrix::conjTrans    = 'c';
char AMatrix::trans        = 't';
char AMatrix::notConjTrans = 'n';

AMatrix::AMatrix(const LongInt& m, const LongInt& n)
 {
    (*this)
     .setNumberOfRows(m)
     .setNumberOfColumns(n)
     .disableTranspose()
     .setType(isAMatrix)
    ;
 }


AMatrix::~AMatrix()
 {

 }


AMatrix& AMatrix::operator = (const AMatrix& A)
 {
    const AMatrixType theType = type();

    if(theType==isHMatrixInterface)
     {
        ((HMatrixInterface*)(this))->isConversionOf(A);
     }
    else if(theType==isRMatrix)
     {
        ((RMatrix*)(this))->isConversionOf(A);
     }
    else if(theType==isRkRMatrix)
     {
        ((RkRMatrix*)(this))->isConversionOf(A);
     }
    else if(theType==isCMatrix)
     {
        ((CMatrix*)(this))->isConversionOf(A);
     }
    else if(theType==isRkCMatrix)
     {
        ((RkCMatrix*)(this))->isConversionOf(A);
     }
    else
     {
        throw SimpleException(IString("Warning In AMatrix::operator = (const AMatrix& A), Bad AMatrixType !!!"));
     }

   return (*this);
 }

AMatrix& AMatrix::operator += (const AMatrix* A)
 {

    addMatrix(*A, 0, 0);
   return (*this);
 }


AMatrix& AMatrix::addMatrix (const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex)
 {
    const AMatrixType theType = type();
 
    if(theType==isHMatrixInterface)
     {
        ((HMatrixInterface*)(this))->addMatrix(A, rowIndex, colIndex);
     }
    else if(theType==isRMatrix)
     {
        ((RMatrix*)(this))->addMatrix(A, rowIndex, colIndex);
     }
    else if(theType==isRkRMatrix)
     {
        ((RkRMatrix*)(this))->addMatrix(A, rowIndex, colIndex);
     }
    else if(theType==isCMatrix)
     {
        ((CMatrix*)(this))->addMatrix(A, rowIndex, colIndex);
     }
    else if(theType==isRkCMatrix)
     {
        ((RkCMatrix*)(this))->addMatrix(A, rowIndex, colIndex);
     }
    else
     {
        throw SimpleException(IString("Warning In AMatrix::addMatrix (const AMatrix& A, const LongInt& rowIndex, const LongInt  colIndex), Bad AMatrixType !!!"));
     }

   return (*this);
 }


AMatrix& AMatrix::operator *= (const AMatrix* A)
 {
//! \todo impl.
   return (*this);
 }


AMatrix& AMatrix::addScalMatrix(const LongReal& z, const AMatrix& B, const LongInt& rowIndex, const LongInt& colIndex)
 {
    const AMatrixType theType = type();
 
    if(theType==isHMatrixInterface)
     {
        ((HMatrixInterface*)(this))->addScalMatrix(z, B, rowIndex, colIndex);
     }
    else if(theType==isRMatrix)
     {
        ((RMatrix*)(this))->addScalMatrix(z, B, rowIndex, colIndex);
     }
    else if(theType==isRkRMatrix)
     {
        ((RkRMatrix*)(this))->addScalMatrix(z, B, rowIndex, colIndex);
     }
    else if(theType==isCMatrix)
     {
        ((CMatrix*)(this))->addScalMatrix(z, B, rowIndex, colIndex);
     }
    else if(theType==isRkCMatrix)
     {
        ((RkCMatrix*)(this))->addScalMatrix(z, B, rowIndex, colIndex);
     }
    else
     {
        throw SimpleException(IString("Warning In AMatrix::addScalMatrix(const Complex& z, const AMatrix& B, const LongInt& rowIndex, const LongInt& colIndex), Bad AMatrixType !!!"));
     }
    
   return (*this);
 }


AMatrix& AMatrix::addScalMatrix(const std::complex<double>& z, const AMatrix& B, const LongInt& rowIndex, const LongInt& colIndex)
 {
    const AMatrixType theType = type();
 
    if(theType==isHMatrixInterface)
     {
        ((HMatrixInterface*)(this))->addScalMatrix(z, B, rowIndex, colIndex);
     }
    else if(theType==isCMatrix)
     {
        ((CMatrix*)(this))->addScalMatrix(z, B, rowIndex, colIndex);
     }
    else if(theType==isRkCMatrix)
     {
        ((RkCMatrix*)(this))->addScalMatrix(z, B, rowIndex, colIndex);
     }
    else
     {
        throw SimpleException(IString("Warning In AMatrix::addScalMatrix(const Complex& z, const AMatrix& B, const LongInt& rowIndex, const LongInt& colIndex), Bad AMatrixType !!!"));
     }
    
   return (*this);
 }


AMatrix& AMatrix::operator *= (const std::complex<double>& alpha)
 {
    const AMatrixType theType = type();
 
    if(theType==isHMatrixInterface)
     {
        *(HMatrixInterface*)(this) *= alpha;
     }
    else if(theType==isCMatrix)
     {
        *(CMatrix*)((this)) *= alpha;
     }
    else if(theType==isRkCMatrix)
     {
        *(RkCMatrix*)(this) *= alpha;             
     }
    else
     {
        throw SimpleException(IString("Warning In AMatrix::operator *= (const Complex&), Type Not Implemented !!!"));
     }

   return (*this);
 }


AMatrix& AMatrix::operator *= (const LongReal& alpha)
 {
    const AMatrixType theType = type();
 
    if(theType==isHMatrixInterface)
     {
        *(HMatrixInterface*)(this) *= alpha;
     }
    else if(theType==isRMatrix)
     {
        *(RMatrix*)((this)) *= alpha;
     }
    else if(theType==isRkRMatrix)
     {
        *(RkRMatrix*)(this) *= alpha;             
     }
    else if(theType==isCMatrix)
     {
        *(CMatrix*)((this)) *= alpha;
     }
    else if(theType==isRkCMatrix)
     {
        *(RkCMatrix*)(this) *= alpha;             
     }
    else
     {        
        throw SimpleException(IString("Warning In AMatrix::operator *= (const LongReal&), Type Not Implemented !!!"));
     }

   return (*this);
 }


AMatrix& AMatrix::addIdentityScaledBy(const std::complex<double>& z)
 {
    const AMatrixType theType = type();
 
    if(theType==isHMatrixInterface)
     {
        ((HMatrixInterface*)(this))->addIdentityScaledBy(z);
     }
    else if(theType==isCMatrix)
     {
        ((CMatrix*)(this))->addIdentityScaledBy(z);
     }
    else if(theType==isRkCMatrix)
     {
        ((RkCMatrix*)(this))->addIdentityScaledBy(z);
     }
    else
     {
        throw SimpleException(IString("Warning In AMatrix::operator *= (const Complex&), Type Not Implemented !!!"));
     }

   return (*this);
 }


AMatrix& AMatrix::addIdentityScaledBy(const LongReal& z)
 {
    const AMatrixType theType = type();
 
    if(theType==isHMatrixInterface)
     {
        ((HMatrixInterface*)(this))->addIdentityScaledBy(z);
     }
    else if(theType==isRMatrix)
     {
        ((RMatrix*)(this))->addIdentityScaledBy(z);
     }
    else if(theType==isRkRMatrix)
     {
        ((RkRMatrix*)(this))->addIdentityScaledBy(z);
     }
    else if(theType==isCMatrix)
     {
        ((CMatrix*)(this))->addIdentityScaledBy(z);
     }
    else if(theType==isRkCMatrix)
     {
        ((RkCMatrix*)(this))->addIdentityScaledBy(z);
     }
    else
     {
        throw SimpleException(IString("Warning In AMatrix::operator *= (const LongReal&), Type Not Implemented !!!"));
     }

   return (*this);
 }


AMatrix&  AMatrix::isConversionOf (const AMatrix& A)
 {
    const AMatrixType theType = type();

    if(theType==isHMatrixInterface)
     {
        ((HMatrixInterface*)(this))->isConversionOf(A);
     }
    else if(theType==isRMatrix)
     {
        ((RMatrix*)(this))->isConversionOf(A);
     }
    else if(theType==isRkRMatrix)
     {
        ((RkRMatrix*)(this))->isConversionOf(A);
     }
    else if(theType==isCMatrix)
     {
        ((CMatrix*)(this))->isConversionOf(A);
     }
    else if(theType==isRkCMatrix)
     {
        ((RkCMatrix*)(this))->isConversionOf(A);
     }
    else
     {
        throw SimpleException(IString("Warning In AMatrix::isConversionOf (const AMatrix& A), Bad AMatrixType !!!"));
     }

   return (*this);
 }


AMatrix& AMatrix::update(const AMatrix* A, const AMatrix* B)
 {
    const AMatrixType theType = type();

    if(theType==isHMatrixInterface)
     {
        ((HMatrixInterface*) this)->update(*A,*B);
     }
    else if(theType==isRMatrix)
     {
        ((RMatrix*) this)->update(*A,*B);
     }
    else if(theType==isRkRMatrix)
     {
        ((RkRMatrix*) this)->update(*A,*B);       
     }
    else if(theType==isCMatrix)
     {
        ((CMatrix*) this)->update(*A,*B);
     }
    else if(theType==isRkCMatrix)
     {        
        ((RkCMatrix*) this)->update(*A,*B);       
     }
    else
     {
        throw SimpleException(IString("Warning In AMatrix::update(const AMatrix* A, const AMatrix* B), Type Not Implemented !!!"));
     }

   return (*this);
 }


AMatrix& AMatrix::isMultiplicationOf(const AMatrix* A, const AMatrix* B)
 {
    const AMatrixType theType = type();

    if(theType==isHMatrixInterface)
     {
        ((HMatrixInterface*) this)->isMultiplicationOf(*A,*B);
     }
    else if(theType==isRMatrix)
     {
        //! \todo setNullMatrix();
        setNullMatrix();
        ((RMatrix*) this)->update(*A,*B);
     }
    else if(theType==isRkRMatrix)
     {
        ((RkRMatrix*) this)->isMultiplicationOf(*A,*B);       
     }
    else if(theType==isCMatrix)
     {
        //! \todo setNullMatrix();
        setNullMatrix();
        ((CMatrix*) this)->update(*A,*B);
     }
    else if(theType==isRkCMatrix)
     {
        ((RkCMatrix*) this)->isMultiplicationOf(*A,*B);       
     }
    else
     {
        throw SimpleException(IString("Warning In AMatrix::isMultiplicationOf(const AMatrix* A, const AMatrix* B), Type Not Implemented !!!"));
     }

   return (*this);
 }


AMatrix& AMatrix::isInverseOf(const AMatrix& S, AMatrix& W)
 {
    const AMatrixType theType = type();

    if(theType==isHMatrixInterface)
     {
        ((HMatrixInterface*) this)->isInverseOf(S, W);
     }
    else if(theType==isRMatrix)
     {
        ((RMatrix*) this)->isInverseOf(S, W);
     }
    else if(theType==isCMatrix)
     {
        ((CMatrix*) this)->isInverseOf(S, W);
     }
    else
     {
        throw SimpleException(IString("Warning In AMatrix::isInverseOf(const AMatrix& S, const AMatrix& W), Type Not Implemented or Imposilble !!!"));
     }
    
   return (*this);
 }


AMatrix& AMatrix::convert2ResolventeAt (const std::complex<double>& z, AMatrix& A)
 {
    A *= -1.0;
    A.addIdentityScaledBy(z);

    this->setNullMatrix();
//! \todo noch machen
//    this->isInverseOf(A, W);

   return (*this);
 }


AMatrix& AMatrix::fastConvert2ResolventeAt (const std::complex<double>& z, AMatrix& A, AMatrix& W)
 {
    A *= -1.0;
    A.addIdentityScaledBy(z);
    this->setNullMatrix();
    this->isInverseOf(A, W);

   return (*this);
 }


CVector  AMatrix::operator * (const CVector& v) const
 {
    const LongInt dim = this->numberOfRows();
    CVector value(dim);

    evaluateAt(v, value);

   return value;
 }


bool AMatrix::addEvaluateAt(const RVector& x, RVector& valueAt) const
 {
    this->operator () (0, valueAt, 0, x);

   return true;
 }


bool AMatrix::evaluateAt(const RVector& x, RVector& valueAt) const
 {
    valueAt.setNull();
    this->operator () (0, valueAt, 0, x);

   return true;
 }


bool AMatrix::addEvaluateAt(const CVector& x, CVector& valueAt) const
 {
    this->operator () (0, valueAt, 0, x);

   return true;
 }


bool AMatrix::evaluateAt(const CVector& x, CVector& valueAt) const
 {
    valueAt.setNull();
    this->operator () (0, valueAt, 0, x);

   return true;
 }


bool AMatrix::operator()(const LongInt& i, CVector& w, const LongInt& j, const CVector& v) const
 {
    const AMatrixType theType = type();
 
    if(theType==isHMatrixInterface)
     {
        (*(HMatrixInterface*)(this))(i,w,j,v);
     }
    else if(theType==isCMatrix)
     {
        (*(CMatrix*)((this)))(i,w,j,v);
     }
    else if(theType==isRkCMatrix)
     {
        (*(RkCMatrix*)(this))(i,w,j,v);             
     }
    else
     {
        throw SimpleException(IString("Warning In AMatrix::operator operator()(const LongInt& i, CVector& w, const LongInt& j, const CVector& v), Type Not Implemented !!!"));
     }

   return true;
 }


bool AMatrix::operator()(const LongInt& i, RVector& w, const LongInt& j, const RVector& v) const
 {
    const AMatrixType theType = type();
 
    if(theType==isHMatrixInterface)
     {
        (*(HMatrixInterface*)(this))(i,w,j,v);
     }
    else if(theType==isRMatrix)
     {
        (*(RMatrix*)((this)))(i,w,j,v);
     }
    else if(theType==isRkRMatrix)
     {
        (*(RkRMatrix*)(this))(i,w,j,v);             
     }
    else
     {
        throw SimpleException(IString("Warning In AMatrix::operator operator()(const LongInt& i, CVector& w, const LongInt& j, const CVector& v), Type Not Implemented !!!"));
     }

   return true;
 }


bool AMatrix::addEvaluateTransposeAt(const RVector& x, RVector& valueAt) const
 {
    partialEvaluateTransposeAt(0, valueAt, 0, x);
   return true;
 }


bool AMatrix::evaluateTransposeAt(const RVector& x, RVector& valueAt) const
 {
    valueAt.setNull();
    partialEvaluateTransposeAt(0, valueAt, 0, x);
   return true;
 }


bool AMatrix::partialEvaluateTransposeAt (const LongInt& i, RVector& w, const LongInt& j, const RVector& v) const
 {
    const AMatrixType theType = type();
 
    if(theType==isHMatrixInterface)
     {
        ((HMatrixInterface*)(this))->partialEvaluateTransposeAt(i,w,j,v);
     }
    else if(theType==isRMatrix)
     {
        ((RMatrix*)((this)))->partialEvaluateTransposeAt(i,w,j,v);
     }
    else if(theType==isRkRMatrix)
     {
        ((RkRMatrix*)(this))->partialEvaluateTransposeAt(i,w,j,v);             
     }
    else
     {
        throw SimpleException(IString("Warning In AMatrix::partialEvaluateTransposeAt (const LongInt& i, RVector& w, const LongInt& j, const RVector& v) const, Type Not Implemented !!!"));
     }
    
   return true;
 }


bool AMatrix::addEvaluateTransposeAt(const CVector& x, CVector& valueAt) const
 {
    partialEvaluateTransposeAt(0, valueAt, 0, x);
   return true;
 }


bool AMatrix::evaluateTransposeAt(const CVector& x, CVector& valueAt) const
 {
    valueAt.setNull();
    partialEvaluateTransposeAt(0, valueAt, 0, x);
   return true;
 }


bool AMatrix::partialEvaluateTransposeAt (const LongInt& i, CVector& w, const LongInt& j, const CVector& v) const
 {
    const AMatrixType theType = type();
 
    if(theType==isHMatrixInterface)
     {
        ((HMatrixInterface*)(this))->partialEvaluateTransposeAt(i,w,j,v);
     }
    else if(theType==isCMatrix)
     {
        ((CMatrix*)((this)))->partialEvaluateTransposeAt(i,w,j,v);
     }
    else if(theType==isRkCMatrix)
     {
        ((RkCMatrix*)(this))->partialEvaluateTransposeAt(i,w,j,v);             
     }
    else
     {
        throw SimpleException(IString("Warning In AMatrix::partialEvaluateTransposeAt (const LongInt& i, CVector& w, const LongInt& j, const CVector& v) const, Type Not Implemented !!!"));
     }
    
   return true;
 }


bool AMatrix::addEvaluateConjugateTransposeAt(const CVector& x, CVector& valueAt) const
 {
    partialEvaluateConjugateTransposeAt(0, valueAt, 0, x);
   return true;
 }


bool AMatrix::evaluateConjugateTransposeAt(const CVector& x, CVector& valueAt) const
 {
    valueAt.setNull();
    partialEvaluateConjugateTransposeAt(0, valueAt, 0, x);
   return true;
 }


bool AMatrix::partialEvaluateConjugateTransposeAt (const LongInt& i, CVector& w, const LongInt& j, const CVector& v) const
 {
    const AMatrixType theType = type();
 
    if(theType==isHMatrixInterface)
     {
        ((HMatrixInterface*)(this))->partialEvaluateConjugateTransposeAt(i,w,j,v);
     }
    else if(theType==isCMatrix)
     {
        ((CMatrix*)((this)))->partialEvaluateConjugateTransposeAt(i,w,j,v);
     }
    else if(theType==isRkCMatrix)
     {
        ((RkCMatrix*)(this))->partialEvaluateConjugateTransposeAt(i,w,j,v);             
     }
    else
     {
        throw SimpleException(IString("Warning In AMatrix::partialEvaluateTransposeAt (const LongInt& i, CVector& w, const LongInt& j, const CVector& v) const, Type Not Implemented !!!"));
     }
    
   return true;
 }


LongReal norm2ProdMinusId (const AMatrix& Si, const AMatrix& S, const Int& nr)
 {    
    const AMatrixType theMatrixSi = Si.type();  

    bool isComplexArithmetic = false;

    if(theMatrixSi==isCMatrix || theMatrixSi==isRkCMatrix)
     {
        isComplexArithmetic = true;
     }
    else if(theMatrixSi==isRMatrix || theMatrixSi==isRkRMatrix)
     {
        isComplexArithmetic = false;
     }
    else if(theMatrixSi==isHMatrixInterface)
     {
        const HMatrixInterface& Hi = (const HMatrixInterface& ) Si;

        isComplexArithmetic = Hi.isComplexArithmetic();
     }
    else
     {
        throw SimpleException(IString("Warning In friend AMatrix LongReal norm2ProdMinusId (const AMatrix& Si, const AMatrix& S, const Int& nr)"));
     }

    const LongInt rows = Si.numberOfRows();
    const LongInt cols = Si.numberOfColumns();

    LongReal norm  = 1.0;
    LongReal value = 1.0;

    if(isComplexArithmetic==true)
     {
        CVector v(rows), w(cols), v2(rows), v3(rows);

        srand( (unsigned)time( NULL ) );

        for(LongInt i=0; i<rows; i++)
         {
            v(i) = std::complex<double>(-LongReal(rand()), LongReal(rand())*1.0e-3 );
         }

        value = 1.0/L2(v);
        v *= value;

        for(LongInt j=0; j<nr; j++)
         {
            Si.evaluateConjugateTransposeAt(v, w);
            S.evaluateConjugateTransposeAt(w, v3);

            v2  = v;
            v2 -= v3;
            v3 -= v;

            S.evaluateAt(v3, w);
            Si.evaluateAt(w, v3);

            v2 += v3;
 
            norm = sqrt(std::abs(innerProduct(v2,v)));

            value = L2(v2);
            if(EPS_NULL<value)
             {
                value = 1.0/value;

                v2 *= value;
                v   = v2;
             }
            else
             {
                j=nr;
             }
         }
     }
    else
     {
        RVector v(rows), w(cols), v2(rows), v3(rows);

        srand( (unsigned)time( NULL ) );

        for(LongInt i=0; i<rows; i++)
         {
            v(i) = (LongReal)(-rand());
         }

        value = 1.0/L2(v);
        v *= value;

        for(LongInt j=0; j<nr; j++)
         {
            Si.evaluateTransposeAt(v, w);
            S.evaluateTransposeAt(w, v3);

            v2  = v;
            v2 -= v3;
            v3 -= v;

            S.evaluateAt(v3, w);
            Si.evaluateAt(w, v3);

            v2 += v3;
 
            norm = sqrt(fabs(innerProduct(v2,v)));

            value = L2(v2);
            if(EPS_NULL<value)
             {
                value = 1.0/value;

                v2 *= value;
                v   = v2;
             }
            else
             {
                j=nr;
             }
         }
     }
   return norm;
 }


LongReal norm2Diff (const AMatrix& A, const AMatrix& B, const Int& nr)
 {
    const LongInt rows = A.numberOfRows();
    const LongInt cols = A.numberOfColumns();

    const AMatrixType theMatrixA = A.type();  

    bool isComplexArithmetic = false;

    if(theMatrixA==isCMatrix || theMatrixA==isRkCMatrix)
     {
        isComplexArithmetic = true;
     }
    else if(theMatrixA==isRMatrix || theMatrixA==isRkRMatrix)
     {
        isComplexArithmetic = false;
     }
    else if(theMatrixA==isHMatrixInterface)
     {
        const HMatrixInterface& H = (const HMatrixInterface& ) A;

        isComplexArithmetic = H.isComplexArithmetic();
     }
    else
     {
        throw SimpleException(IString("Warning In friend AMatrix LongReal norm2ProdMinusId (const AMatrix& Si, const AMatrix& S, const Int& nr)"));
     }

    LongReal norm  = 1.0;
    LongReal value = 1.0;

    if(isComplexArithmetic==true)
     {
        CVector v(cols), w(rows), v2(cols);

        srand( (unsigned)time( NULL ) );

        for(LongInt i=0; i<cols; i++)
         {
            v(i) = std::complex<double>(-LongReal(rand()), LongReal(rand())*1.0e-9 );
         }

        value = 1.0/L2(v);
        v *= value;

        for(LongInt j=0; j<nr; j++)
         {
            A.evaluateAt(v, w);
            w *= -1.0;
            B.addEvaluateAt(v, w);

            A.evaluateConjugateTransposeAt(w, v2);
            v2 *= -1.0;
            B.addEvaluateConjugateTransposeAt(w, v2);

            norm = sqrt(std::abs(innerProduct(v2,v)));

            value = L2(v2);
            if(EPS_NULL<value)
             {
                value = 1.0/value;

                v2 *= value;
                v   = v2;
             }
            else
             {
                j=nr;
             }
         }
     }
    else
     {
        RVector v(cols), w(rows), v2(cols);

        srand( (unsigned)time( NULL ) );

        for(LongInt i=0; i<cols; i++)
         {
            v(i) = -LongReal(rand());
         }

        value = 1.0/L2(v);
        v *= value;

        for(LongInt j=0; j<nr; j++)
         {
            A.evaluateAt(v, w);
            w *= -1.0;
            B.addEvaluateAt(v, w);

            A.evaluateTransposeAt(w, v2);
            v2 *= -1.0;
            B.addEvaluateTransposeAt(w, v2);

            norm = fabs(innerProduct(v2,v));

            value = L2(v2);
            if(EPS_NULL<value)
             {
                value = 1.0/value;

                v2 *= value;
                v   = v2;
             }
            else
             {
                j=nr;
             }
         }
     }

   return sqrt(norm);
 }


LongReal norm2(const AMatrix& A, const Int& nr)
 {
    const LongInt rows = A.numberOfRows();
    const LongInt cols = A.numberOfColumns();

    const AMatrixType theMatrixA = A.type();  

    bool isComplexArithmetic = false;

    if(theMatrixA==isCMatrix || theMatrixA==isRkCMatrix)
     {
        isComplexArithmetic = true;
     }
    else if(theMatrixA==isRMatrix || theMatrixA==isRkRMatrix)
     {
        isComplexArithmetic = false;
     }
    else if(theMatrixA==isHMatrixInterface)
     {
        const HMatrixInterface& H = (const HMatrixInterface& ) A;

        isComplexArithmetic = H.isComplexArithmetic();
     }
    else
     {
        throw SimpleException(IString("Warning In friend AMatrix LongReal norm2ProdMinusId (const AMatrix& Si, const AMatrix& S, const Int& nr)"));
     }

    LongReal norm  = 1.0;
    LongReal value = 1.0;

    if(isComplexArithmetic==true)
     {
        CVector v(cols), w(rows), v2(cols);

        srand( (unsigned)time( NULL ) );

        for(LongInt i=0; i<cols; i++)
         {
            v(i) = std::complex<double>(-LongReal(rand()), LongReal(rand())*1.0e-9 );
         }

        value = 1.0/L2(v);
        v *= value;

        for(LongInt j=0; j<nr; j++)
         {
            A.evaluateAt(v, w);
            A.evaluateConjugateTransposeAt(w, v2);

            norm = std::abs(innerProduct(v2,v));

            value = L2(v2);
            if(EPS_NULL<value)
             {
                value = 1.0/value;

                v2 *= value;
                v   = v2;
             }
            else
             {
                j=nr;
             }
         }
     }
    else
     {
        RVector v(cols), w(rows), v2(cols);

        srand( (unsigned)time( NULL ) );
         
        v.setRand();

        value = 1.0/L2(v);
        v *= value;

        for(LongInt j=0; j<nr; j++)
         {
            A.evaluateAt(v, w);
            A.evaluateTransposeAt(w, v2);

            norm = innerProduct(v2,v);

            value = L2(v2);

            if(EPS_NULL<value)
             {
                value = 1.0/value;

                v2 *= value;
                v   = v2;
             }
            else
             {
                j=nr;
             }
         }
     }

   return sqrt(norm);
 }

/*
double 
norm2nr_supermatrix(psupermatrix s, int nr){
  double norm2 = 0.0;
  double tmp = 0.0;
  double *v, *w, *v2;
  int i,rows,cols;
  assert(s != 0x0);
  assert(s->u == 0x0);
  rows = s->rows;
  cols = s->cols;
  v = (double*) malloc(cols*sizeof(double));
  w = (double*) malloc(rows*sizeof(double));
  v2 = (double*) malloc(cols*sizeof(double));
  assert(v!=0x0 && w!=0x0 && v2!=0x0);
  for(i=0; i<cols; i++){
    v[i] = rand();
  }
  tmp = 1.0/dnrm2_(&cols,v,eins_);
  dscal_(&cols,&tmp,v,eins_);
  for(i=0; i<nr; i++){
    eval_supermatrix(s,v,w);
    evaltrans_supermatrix(s,w,v2);
    norm2 = dnrm2_(&cols,v2,eins_);
    if(tmp>1e-32){
      tmp=1.0/norm2;
      dscal_(&cols,&tmp,v2,eins_);
      dcopy_(&cols,v2,eins_,v,eins_);
    }else{
      i=nr;
    }
  }
  free(v);
  free(v2);
  free(w);
  return sqrt(norm2);
}
*/

/*

  tmp = 1.0/dnrm2_(&cols,v,eins_);

  dscal_(&cols,&tmp,v,eins_);

  for(i=0; i<nr; i++){
    eval_supermatrix(s,v,w);
    evaltrans_supermatrix(s,w,v2);
    norm2 = dnrm2_(&cols,v2,eins_);
    if(tmp>1e-32){
      tmp=1.0/norm2;
      dscal_(&cols,&tmp,v2,eins_);
      dcopy_(&cols,v2,eins_,v,eins_);
    }else{
      i=nr;
    }
  }

*/
LongReal AMatrix::frobeniusNorm () const
 {
    LongReal norm = 0.0;

    const AMatrixType theType = type();
 
    if(theType==isHMatrixInterface)
     {
        norm = ((HMatrixInterface*)(this))->frobeniusNorm ();
     }
    else if(theType==isRMatrix)
     {
        norm = ((RMatrix*)((this)))->frobeniusNorm ();
     }
    else if(theType==isRkRMatrix)
     {
        norm = ((RkRMatrix*)(this))->frobeniusNorm ();
     }
    else if(theType==isCMatrix)
     {
        norm = ((CMatrix*)((this)))->frobeniusNorm ();
     }
    else if(theType==isRkCMatrix)
     {
        norm = ((RkCMatrix*)(this))->frobeniusNorm ();
     }
    else
     {
        throw SimpleException(IString("Warning In AMatrix::frobeniusNorm () const, Type Not Implemented !!!"));
     }    
     
   return norm;
 }
