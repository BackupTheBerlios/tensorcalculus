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
 * \defgroup HMatrix
 * \brief Diese Bibliothek enth�lt grundlegende Dienstklassen f�r H-Matrizen
 *   
 ******************************************************************************
*/


/*!
 ******************************************************************************
 * \class    HMatrixInterface
 *
 * \brief    HMatrixInterface 
 *
 *
 * \author   Mike Espig
 *
 * \ingroup  HMatrix
 *
 *
 ******************************************************************************
 *
 *
 ******************************************************************************
 */

#ifndef __HMatrixInterface__
#define __HMatrixInterface__

#include "Macros.h"
#include <iostream>
#include "ACMatrix.hpp"
#include "RkCMatrix.hpp"
#include "CMatrix.hpp"
#include "RkRMatrix.hpp"
#include "RMatrix.hpp"

using namespace std;

class  HMatrixInterface : public AMatrix
 {
    DECLARE (HMatrixInterface) 

    friend class RkCMatrix;
    friend class RkRMatrix;
    friend class HMatrixFem;
    
    FLAG (ComplexArithmetic)

    CRITICAL_ATTRIBUTE(LongInt, numberOfRowBlocks,    setNumberOfRowBlocks)
    CRITICAL_ATTRIBUTE(LongInt, numberOfColumnBlocks, setNumberOfColumnBlocks)
   

  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  protected:
              HMatrixInterface (const LongInt& m, const LongInt&n,
                                const LongInt& block_m, const LongInt& block_n);

  public:
 
     virtual ~HMatrixInterface ();   


  //@}

  /*!
  ***************************************************************************************
  * \name                             �berladene Operatoren
  ***************************************************************************************/

  //@{

  public:

     HMatrixInterface& operator += (const AMatrix& A);
     HMatrixInterface& operator *= (const LongReal& alpha);
     HMatrixInterface& operator *= (const std::complex<double>& alpha);
     HMatrixInterface& operator /= (const std::complex<double>& alpha);
     
     bool              operator () (const LongInt& i, RVector& w, const LongInt& j, const RVector& v)  const;
     bool              operator () (const LongInt& i, CVector& w, const LongInt& j, const CVector& v)  const;

     AMatrixPointer&   operator () (const LongInt i, const LongInt j) const 
                        {
                           const LongInt index = i+j*attr_numberOfRowBlocks; 
                          return (attr_block[index]); 
                        }

     friend ostream& operator << (ostream & s, const HMatrixInterface& H);
     //friend istream& operator >> (istream & s, CMatrix &A);

  //@}

  /*!
  ***************************************************************************************
  * \name                        �ffentliche Klassendienste
  ***************************************************************************************/

  //@{

  public:

     HMatrixInterface& update              (const AMatrix& A,  const AMatrix& B);
     HMatrixInterface& addScalMatrix       (const LongReal& z, const AMatrix& B, const LongInt& rowIndex=0, const LongInt& colIndex=0);
     HMatrixInterface& addScalMatrix       (const std::complex<double>& z, const AMatrix& B, const LongInt& rowIndex=0, const LongInt& colIndex=0);
     HMatrixInterface& isMultiplicationOf  (const AMatrix& A, const AMatrix& B);
     HMatrixInterface& isInverseOf         (const AMatrix& S, AMatrix& W);
     HMatrixInterface& isConversionOf      (const AMatrix& A);
     HMatrixInterface& addIdentityScaledBy (const std::complex<double>& z);
     HMatrixInterface& addIdentityScaledBy (const LongReal& z);
     HMatrixInterface& addMatrix           (const AMatrix& A, const LongInt& rowIndex=0, const LongInt& colIndex=0);
     LongInt           rowIndex            (const LongInt& i, const LongInt& j) const;
     LongInt           columnIndex         (const LongInt& i, const LongInt& j) const;
     AMatrix&          setNullMatrix       ();
     
     LongReal frobeniusNorm                () const;

     bool evaluateTransposeAt                 (const CVector& x, CVector& valueAt) const;
     bool partialEvaluateTransposeAt          (const LongInt& i, CVector& w, const LongInt& j, const CVector& v) const;

     bool evaluateTransposeAt                 (const RVector& x, RVector& valueAt) const;
     bool partialEvaluateTransposeAt          (const LongInt& i, RVector& w, const LongInt& j, const RVector& v) const;

     bool evaluateConjugateTransposeAt        (const CVector& x, CVector& valueAt) const;
     bool partialEvaluateConjugateTransposeAt (const LongInt& i, CVector& w, const LongInt& j, const CVector& v) const;

  //@}


  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{

  protected:

     AMatrixPointer* attr_block;

  private:




  //@}
 };


typedef  HMatrixInterface* HMatrixInterfacePointer;

#endif // not defined __HMatrixInterface__
