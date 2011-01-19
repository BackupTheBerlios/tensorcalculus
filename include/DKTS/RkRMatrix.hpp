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
 * \class    RkRMatrix
 *
 * \brief    RkRMatrix aus IR(m, n)
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

#ifndef __RkRMatrix__
#define __RkRMatrix__

#include "Macros.h"
#include <iostream>
#include "RMatrix.hpp"
#include "LongRealSvd.hpp"
#include "LongRealQR.hpp"
#include "SimpleException.hpp"
#include "ARMatrix.hpp"
#include "HMatrixInterface.hpp"

using namespace std;

class  RkRMatrix : public ARMatrix
 {
    friend class RMatrix;
    friend class HMatrixModelProblem1D;

    friend class KTS2;

    DECLARE (RkRMatrix)
    
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
    RkRMatrix (const LongInt n, const LongInt m, const LongInt rank);
    RkRMatrix (const LongInt n, const LongInt m, const LongInt rank, const bool& flag);

    RkRMatrix (const RkRMatrix& aRkMatrix);

    RkRMatrix (const RMatrix& Z, const LongInt rank);

    //! Destruktor
    virtual ~RkRMatrix ();

  //@}

  /*!
  ***************************************************************************************
  * \name                             �berladene Operatoren
  ***************************************************************************************/

  //@{

  public:

    RVector    operator  * (const RVector& x);
    RkRMatrix& operator  = (const RkRMatrix& aRkMatrix);
    RkRMatrix& operator += (const RkRMatrix& aRkMatrix);
    RkRMatrix& operator += (const RMatrix&   aMatrix);
    RkRMatrix& operator *= (RkRMatrix& aRkMatrix);
    RkRMatrix& operator *= (RMatrix& aRMatrix);
    RkRMatrix& operator *= (const LongReal& alpha); 
    RkRMatrix& operator /= (const LongReal& alpha);
    bool       operator () (const LongInt& i, RVector& w, const LongInt& j, const RVector& v); 

    friend ostream&  operator << (ostream& s, const RkRMatrix& R);
    friend istream&  operator >> (istream& s, RkRMatrix& R);
    friend RkRMatrix operator + (const RkRMatrix& aRk, const RkRMatrix& aRkMatrix);
    friend RkRMatrix operator * (const RkRMatrix& aRk, RkRMatrix& aRkMatrix);
    friend RkRMatrix operator * (const RkRMatrix& aRk, RMatrix&   aRMatrix);
    friend RkRMatrix operator * (const RMatrix&   aRMatrix, RkRMatrix& aRk);
    friend RkRMatrix operator * (const RkRMatrix& aRk, const LongReal& alpha); 
    friend RkRMatrix operator / (const RkRMatrix& aRk, const LongReal& alpha); 

  //@}

  /*!
  ***************************************************************************************
  * \name                        �ffentliche Klassendienste
  ***************************************************************************************/

  //@{
    
  public:
     
     virtual AMatrix&  setNullMatrix ();
  
     RMatrix    A               () const  {return attr_A;}
     RMatrix    B               () const  {return attr_B;}    
     RkRMatrix& setA            (const RMatrix& aRMatrix);
     RkRMatrix& setB            (const RMatrix& aRMatrix);
     RkRMatrix& setAB           (const RMatrix& aRMatrix1, const RMatrix& aRMatrix2);
     RkRMatrix& setDimension    (const LongInt& m, const LongInt& n, const LongInt& k);
     RkRMatrix& turncation2     (const LongInt newRank);
     RkRMatrix& isConversionOf  (const AMatrix& A);
     RMatrix    convert2RMatrix ();
     LongReal   frobeniusNorm   ()const;

     bool       partialEvaluateTransposeAt (const LongInt& i, RVector& w, const LongInt& j, const RVector& v);
     bool       partialevaluateTransposeAt (const LongInt& i, RVector& w, const LongInt& j, const RVector& v);
     RkRMatrix& update                     (const AMatrix& A, const AMatrix& B);
     RkRMatrix& isMultiplicationOf         (const AMatrix& A, const AMatrix& B);
     RkRMatrix& addMatrix                  (const AMatrix& A, const LongInt& rowIndex=0, const LongInt& colIndex=0);
     RkRMatrix& addIdentityScaledBy        (const LongReal& z);
     RkRMatrix& addScalMatrix              (const LongReal& z, const AMatrix& B, const LongInt& rowIndex=0, const LongInt& colIndex=0);
  //@}


  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{

  private:

    bool isPreparedConversionOf (const HMatrixInterface& H, const LongInt& newRank, const LongInt& maxRank);

    RMatrix attr_A;
    RMatrix attr_B;

    //! 

  //@}
 };

#endif // not defined __RkRMatrix__
