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

/*!
 ******************************************************************************
 * \class    AMatrix
 *
 * \brief     
 *
 *
 * \author   Mike Espig
 *
 * \ingroup  MathBase
 *
 *
 ******************************************************************************
 *
 *
 ******************************************************************************
 */

#ifndef __AMatrix__
#define __AMatrix__

#include "Macros.h"
#include "SimpleException.hpp"
#include "Complex.hpp"
#include "CVector.hpp"

enum AMatrixType
 {	
    isAMatrix,
    isARMatrix,
	   isACMatrix,
    isRMatrix,
    isRkRMatrix,
	   isCMatrix,
	   isRkCMatrix,
	   isHMatrixInterface
 };


class AMatrix
 {
    friend class CMatrix;
    friend class RMatrix;
    friend class RkCMatrix;
    friend class RkRMatrix;
    friend class HMatrixInterface;
    friend class HMatrixFem;

    friend class SumOfRMatrixKroneckerProducts;
    friend class KTS2;

    DECLARE (AMatrix)

    FLAG (Transpose)

    CRITICAL_ATTRIBUTE (LongInt,     numberOfRows,    setNumberOfRows)
    CRITICAL_ATTRIBUTE (LongInt,     numberOfColumns, setNumberOfColumns)
    CRITICAL_ATTRIBUTE (AMatrixType, type,            setType)
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:

     static char  conjTrans;
     static char  trans;
     static char  notConjTrans;
   
     AMatrix (const LongInt& m, const LongInt& n);
     virtual ~AMatrix ();
  

  //@}

  /*!
  ***************************************************************************************
  * \name                             �berladene Operatoren
  ***************************************************************************************/

  //@{

  public:

     AMatrix& operator  = (const AMatrix& A);
     AMatrix& operator += (const AMatrix* A);
     AMatrix& operator *= (const AMatrix* A);
     AMatrix& operator *= (const std::complex<double>& alpha);
     AMatrix& operator *= (const LongReal& alpha);
     bool     operator () (const LongInt& i, RVector& w, const LongInt& j, const RVector& v) const;
     bool     operator () (const LongInt& i, CVector& w, const LongInt& j, const CVector& v) const;
     CVector  operator *  (const CVector& v) const ;

  //@}

  /*!
  ***************************************************************************************
  * \name                        �ffentliche Klassendienste
  ***************************************************************************************/

  //@{

  public:

     virtual AMatrix& setNullMatrix              () = 0;
             AMatrix& update                     (const AMatrix* A, const AMatrix* B);
             AMatrix& isMultiplicationOf         (const AMatrix* A, const AMatrix* B);
             AMatrix& isInverseOf                (const AMatrix& S, AMatrix& W);
             AMatrix& isConversionOf             (const AMatrix& A);
             
             AMatrix& addMatrix                  (const AMatrix& A, const LongInt& rowIndex=0,
                                                  const LongInt& colIndex=0);                                           
             AMatrix& addScalMatrix              (const LongReal& z, const AMatrix& B, const LongInt& rowIndex=0, 
                                                  const LongInt& colIndex=0);
             AMatrix& addScalMatrix              (const std::complex<double>& z, const AMatrix& B, const LongInt& rowIndex=0,
                                                  const LongInt& colIndex=0);
                                                  
             AMatrix& addIdentityScaledBy        (const LongReal& z);
             AMatrix& addIdentityScaledBy        (const std::complex<double>& z);
             AMatrix& convert2ResolventeAt       (const std::complex<double>& z, AMatrix& A);
             AMatrix& fastConvert2ResolventeAt   (const std::complex<double>& z, AMatrix& A, AMatrix& W);
             
             LongReal frobeniusNorm              () const;

             bool     addEvaluateAt              (const RVector& x, RVector& valueAt) const;  
             bool     evaluateAt                 (const RVector& x, RVector& valueAt) const;  
             bool     addEvaluateAt              (const CVector& x, CVector& valueAt) const;  
             bool     evaluateAt                 (const CVector& x, CVector& valueAt) const;  

             bool     addEvaluateTransposeAt     (const RVector& x, RVector& valueAt) const;
             bool     evaluateTransposeAt        (const RVector& x, RVector& valueAt) const;
             bool     partialEvaluateTransposeAt (const LongInt& i, RVector& w, const LongInt& j, const RVector& v) const;           
             bool     addEvaluateTransposeAt     (const CVector& x, CVector& valueAt) const;
             bool     evaluateTransposeAt        (const CVector& x, CVector& valueAt) const;
             bool     partialEvaluateTransposeAt (const LongInt& i, CVector& w, const LongInt& j, const CVector& v) const;

             bool     addEvaluateConjugateTransposeAt     (const CVector& x, CVector& valueAt) const;
             bool     evaluateConjugateTransposeAt        (const CVector& x, CVector& valueAt) const;
             bool     partialEvaluateConjugateTransposeAt (const LongInt& i, CVector& w, const LongInt& j, const CVector& v) const;

     friend  LongReal norm2              (const AMatrix& A, const Int& nr=30);
     friend  LongReal norm2Diff          (const AMatrix& A, const AMatrix& B, const Int& nr=30);     
     friend  LongReal norm2ProdMinusId   (const AMatrix& A, const AMatrix& B, const Int& nr=30);

  //@}
  

  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{

  private:


  //@}
 };


typedef  AMatrix* AMatrixPointer;

#endif // not defined __AMatrix__
