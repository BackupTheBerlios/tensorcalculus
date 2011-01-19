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
 * \class    RMatrix
 *
 * \brief    allgemeine RMatrix aus IR(m, n)
 *
 *
 * \author   Mike Espig
 *
 * \ingroup  Mathbase
 *
 ******************************************************************************
 *
 *
 ******************************************************************************
 */

#ifndef __RMatrix__
#define __RMatrix__


#include "RVector.hpp"
#include "RMatrix.hpp"
#include "ARMatrix.hpp"


class RkRMatrix;

class  RMatrix : public ARMatrix
 {

  friend class RVector;
  friend class HMatrixFem;

  DECLARE (RMatrix)

  FLAG (ConjugateTranspose)


  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
    //! Konstruktor
    /*! Dies ist der Standardkonstruktor
    \param n die Anzahl der Spalten
    \param m die Anzahl der Zeilen                                 */
    RMatrix (LongInt n=0, LongInt m=0);

    //! Konstruktor
    /*! Dies ist der Copy-Konstruktor                              */
    RMatrix (const RMatrix &A);

    //! Destruktor
    virtual ~RMatrix ();

  //@}

  /*!
  ***************************************************************************************
  * \name                             �berladene Operatoren
  ***************************************************************************************/

  //@{

  public:
    bool       operator == (const RMatrix& A) const;
    RMatrix&   operator  = (const RMatrix& A);
    RMatrix&   operator *= (RMatrix& A);
    RMatrix&   operator *= (const RkRMatrix& R);    
    RMatrix&   operator += (const RMatrix& A);
    RMatrix&   operator -= (const RMatrix& A);
    RMatrix&   operator *= (const LongReal& alpha);
    RMatrix&   operator /= (const LongReal& alpha);
    RMatrix    operator  + (const RMatrix& A);
    RMatrix    operator  - (const RMatrix& A);
    RMatrix    operator  * (RMatrix& A);
    RVector    operator  * (const RVector& v);
    bool       operator () (const LongInt& i, RVector& w, const LongInt& j, const RVector& v);
   
    RMatrix&   leftMultiplied(RMatrix& A);

    inline LongReal&   operator () (const LongInt i, const LongInt j) const 
                            {
                               const LongInt index = i+j*attr_numberOfRows; 
                              return (_pelm[index]); 
                            }

  //@}

  /*!
  ***************************************************************************************
  * \name                        �ffentliche Klassendienste
  ***************************************************************************************/

  //@{

  public:

    //! tansponiert diese RMatrix
    RMatrix&  transposed();

    //! tansponiert und konjugiert diese RMatrix
    RMatrix&  conjugateTransposed ();

    //! konjugiert diese RMatrix
    RMatrix&  conjugated ();

    //! setzt alle RMatrix-Komponenten auf 0
    virtual AMatrix&  setNullMatrix ();

    //! setzt die Gr��e der RMatrix neu
    /*!
    \param n die Anzahl der Spalten
    \param m die Anzahl der Zeilen                                 */
    RMatrix&	 setDimension   (LongInt n, LongInt m);

    RMatrix&	 setTranspose   (bool b = true);

    //! liefert einen SpaltenRVector
    /*
    \param k der 0-basierte Index des SpaltenRVectors */
    RVector   getColumn      (LongInt k);

    //! setzt einen SpaltenRVector
    /*
    \param k der 0-basierte Index des SpaltenRVectors
    \param v der RVector, auf den der SpaltenRVector gesetzt werden soll */
    RMatrix&  setColumn      (LongInt k, const RVector &v);

    //! liefert einen ZeilenRVector
    /*
    \param i der 0-basierte Index des ZeilenRVectors */
    RVector   getRow         (LongInt i);

    //! setzt einen ZeilenRVector
    /*
    \param i der 0-basierte Index des ZeilenRVectors
    \param v der RVector, auf den der ZeilenRVector gesetzt werden soll */
    RMatrix&  setRow         (LongInt i, const RVector &v);

    RMatrix&  setDiag   (const RMatrix& A);
    RMatrix&  setRand   (const LongReal& eps=2.5e-6);    
    RMatrix&  setSymRand(const LongReal& eps=2.5e-6);
    RMatrix&  setId     (const LongReal& alpha=1.0);
    RMatrix&  setIdEps  (const LongReal& alpha=1.0, const LongReal& eps=0.5);

    RMatrix&  update    (const AMatrix& A, const AMatrix& B);
    RMatrix&  update    (const LongReal& a, const RMatrix& B);
    RMatrix&  update    (const LongReal& a, RMatrix& A, RMatrix& B, 
                         const LongReal& b, const LongInt& rowOffset=0, const LongInt& colOffset=0);

    RMatrix&  addMatrix (const AMatrix& A, const LongInt& rowIndex=0, const LongInt& colIndex=0);
    
    RMatrix&  addScalMatrix              (const LongReal& z, const AMatrix& B, const LongInt& rowIndex=0, const LongInt& colIndex=0);
    RMatrix&  addIdentityScaledBy        (const LongReal& z);
    bool      partialEvaluateTransposeAt (const LongInt& i, RVector& w, const LongInt& j, const RVector& v) const;


    RMatrix&  isReorganized   (const RMatrix& A);
    RMatrix&  isConversionOf  (const AMatrix& A);
    RMatrix&  attachColumns   (const RMatrix& A);
    RMatrix&  attachRows      (const RMatrix& A);
    RMatrix&  invert          ();
    RMatrix&  isInverseOf     (const AMatrix& s, AMatrix& w);
    LongReal  frobeniusNorm   () const;
    LongReal  cond2           () const;
    char*     correctChar     () const;

    RMatrix&  dyadicProduct (const LongReal& alpha, const RVector& x, const RVector& y);

  //@}

  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{

  private:

    friend ostream& operator << (ostream & s,const RMatrix &A);
    friend istream& operator >> (istream & s, RMatrix &A);
    friend RMatrix  operator *  (const LongReal alpha, const RMatrix& A);
    friend RMatrix  operator /  (const RMatrix& A, const LongReal alpha);
    friend LongReal  BiLform    (const RVector &v1,const RMatrix& A,const RVector &v2);


    //! Zeiger auf RMatrixelemente
    LongRealPointer _pelm;

  //@}
 };


typedef  RMatrix* RMatrixPointer;

#endif // not defined __RMatrix__
