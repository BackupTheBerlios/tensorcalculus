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
 * \class    RkCMatrix
 *
 * \brief    RkCMatrix aus IC(n,n)
 *
 *
 * \author   Mike Espig
 *
 * \ingroup  Mathbase
 *
 *
 ******************************************************************************
 *
 *
 ******************************************************************************
 */

#ifndef __RkCMatrix__
#define __RkCMatrix__

#include "Macros.h"
#include <iostream>
#include "CMatrix.hpp"
#include "ComplexSvd.hpp"
#include "ComplexQR.hpp"
#include "SimpleException.hpp"
#include "ACMatrix.hpp"
// #include "HMatrixInterface.hpp"

using namespace std;

class HMatrixInterface;

class  RkCMatrix : public ACMatrix
 {
    DECLARE (RkCMatrix)
 
    friend class HMatrixModelProblem1D;
    friend class HMHarmonicModelProblem1D;
    friend class HMatrixHamilton1D;
    friend class CMatrix;
    friend class HMatrixMTest;
    friend class HMatrixBTest;

    CRITICAL_ATTRIBUTE (LongInt, rank,        setRank)

  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
    //! Konstruktor
    /*! Dies ist der Standardkonstruktor
    \param n die Anzahl der Spalten und Zeilen*/
    RkCMatrix (const LongInt n, const LongInt m, const LongInt rank);

    RkCMatrix (const RkCMatrix& aRkMatrix);

    RkCMatrix (const CMatrix& Z, const LongInt rank);

    //! Destruktor
    virtual ~RkCMatrix ();

  //@}

  /*!
  ***************************************************************************************
  * \name                             �berladene Operatoren
  ***************************************************************************************/

  //@{

  public:

    CVector    operator  * (const CVector& x);
    RkCMatrix& operator  = (const RkCMatrix& aRkMatrix);
    RkCMatrix& operator += (const RkCMatrix& aRkMatrix);
    RkCMatrix& operator += (const CMatrix&   aMatrix);
    RkCMatrix& operator *= (RkCMatrix& aRkMatrix);
    RkCMatrix& operator *= (CMatrix& aCMatrix);
    RkCMatrix& operator *= (const std::complex<double>& alpha);
    RkCMatrix& operator /= (const std::complex<double>& alpha);
    bool       operator () (const LongInt& i, CVector& w, const LongInt& j, const CVector& v); 

    friend ostream&  operator << (ostream& s, const RkCMatrix& R);
    friend istream&  operator >> (istream& s, RkCMatrix& R);
    friend RkCMatrix operator + (const RkCMatrix& aRk, const RkCMatrix& aRkMatrix);
    friend RkCMatrix operator * (const RkCMatrix& aRk, RkCMatrix& aRkMatrix);
    friend RkCMatrix operator * (const RkCMatrix& aRk, CMatrix&   aCMatrix);
    friend RkCMatrix operator * (const CMatrix&   aCMatrix, RkCMatrix& aRk);
    friend RkCMatrix operator * (const RkCMatrix& aRk, const std::complex<double>& alpha);
    friend RkCMatrix operator / (const RkCMatrix& aRk, const std::complex<double>& alpha);

  //@}

  /*!
  ***************************************************************************************
  * \name                        �ffentliche Klassendienste
  ***************************************************************************************/

  //@{
    
  public:
     
     virtual AMatrix&  setNullMatrix ();
  
     CMatrix    A               () const  {return attr_A;}
     CMatrix    B               () const  {return attr_B;}    
     RkCMatrix& setA            (const CMatrix& aCMatrix);
     RkCMatrix& setB            (const CMatrix& aCMatrix);
     RkCMatrix& setAB           (const CMatrix& aCMatrix1, const CMatrix& aCMatrix2);
     RkCMatrix& setDimension    (const LongInt& m, const LongInt& n, const LongInt& k);
     RkCMatrix& turncation2     (const LongInt newRank);
     RkCMatrix& isConversionOf  (const AMatrix& A);
     CMatrix    convert2CMatrix ();
     LongReal   frobeniusNorm   ()const;

     bool       partialEvaluateTransposeAt          (const LongInt& i, CVector& w, const LongInt& j, const CVector& v);
     bool       partialEvaluateConjugateTransposeAt (const LongInt& i, CVector& w, const LongInt& j, const CVector& v);
     RkCMatrix& update                              (const AMatrix& A, const AMatrix& B);
     RkCMatrix& isMultiplicationOf                  (const AMatrix& A, const AMatrix& B);
     RkCMatrix& addMatrix                           (const AMatrix& A, const LongInt& rowIndex=0, const LongInt& colIndex=0);
     RkCMatrix& addIdentityScaledBy                 (const std::complex<double>& z);
     RkCMatrix& addScalMatrix                       (const std::complex<double>& z, const AMatrix& B, const LongInt& rowIndex=0, const LongInt& colIndex=0);
  //@}


  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{

  private:

    bool isPreparedConversionOf (const HMatrixInterface& H, const LongInt& newRank, const LongInt& maxRank);

    CMatrix attr_A;
    CMatrix attr_B;

    //! 

  //@}
 };

#endif // not defined __RkCMatrix__
