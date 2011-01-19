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
 * \class    CMatrix
 *
 * \brief    allgemeine CMatrix aus IC(n,m)
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

#ifndef __CMatrix__
#define __CMatrix__


#include "CVector.hpp"
#include "RMatrix.hpp"
#include "ACMatrix.hpp"




class RkCMatrix;

class  CMatrix : public ACMatrix
 {

  friend class CVector;
  friend class HMatrixFem;

  DECLARE (CMatrix)

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
    CMatrix (LongInt n, LongInt m);

    //! Konstruktor
    /*! Dies ist der Copy-Konstruktor                              */
    CMatrix (const CMatrix &A);

    CMatrix (const RMatrix &A);
    //! Destruktor
    virtual ~CMatrix ();

  //@}

  /*!
  ***************************************************************************************
  * \name                             �berladene Operatoren
  ***************************************************************************************/

  //@{

  public:
    bool       operator == (const CMatrix& A) const;
    CMatrix&   operator  = (const CMatrix& A);
    CMatrix&   operator *= (CMatrix& A);
    CMatrix&   operator *= (const RkCMatrix& R);    
    CMatrix&   operator += (const CMatrix& A);
    CMatrix&   operator -= (const CMatrix& A);
    CMatrix&   operator -= (const RMatrix& A);
    CMatrix&   operator *= (const std::complex<double>& alpha);
    CMatrix&   operator /= (const std::complex<double>& alpha);
    CMatrix    operator  + (const CMatrix& A);
    CMatrix    operator  - (const CMatrix& A);
    CMatrix    operator  * (CMatrix& A);
    CVector    operator  * (const CVector& v);
    bool       operator () (const LongInt& i, CVector& w, const LongInt& j, const CVector& v);
   
    CMatrix&   leftMultiplied(CMatrix& A);

    inline std::complex<double>&   operator () (const LongInt i, const LongInt j) const
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

    //! tansponiert und konjugiert diese CMatrix
    CMatrix&  conjugateTransposed ();

    //! konjugiert diese CMatrix
    CMatrix&  conjugated ();

    //! setzt alle CMatrix-Komponenten auf 0
    virtual AMatrix&  setNullMatrix ();

    //! liefert den Realteil der complexen Matrix
    RMatrix realPart() const;

    //! liefert den Imagin�rteil der complexen Matrix
    RMatrix imagPart() const;


    //! setzt die Gr��e der CMatrix neu
    /*!
    \param n die Anzahl der Spalten
    \param m die Anzahl der Zeilen                                 */
    CMatrix&	 setDimension   (LongInt n, LongInt m);

    //! liefert einen SpaltenCVector
    /*
    \param k der 0-basierte Index des SpaltenCVectors */
    CVector   getColumn      (LongInt k);

    //! setzt einen SpaltenCVector
    /*
    \param k der 0-basierte Index des SpaltenCVectors
    \param v der CVector, auf den der SpaltenCVector gesetzt werden soll */
    CMatrix&  setColumn      (LongInt k, const CVector &v);

    //! liefert einen ZeilenCVector
    /*
    \param i der 0-basierte Index des ZeilenCVectors */
    CVector   getRow         (LongInt i);

    //! setzt einen ZeilenCVector
    /*
    \param i der 0-basierte Index des ZeilenCVectors
    \param v der CVector, auf den der ZeilenCVector gesetzt werden soll */
    CMatrix&  setRow         (LongInt i, const CVector &v);


    CMatrix&  update    (const AMatrix& A, const AMatrix& B);
    CMatrix&  update    (const std::complex<double>& a, const CMatrix& B);
    CMatrix&  update    (const std::complex<double>& a, CMatrix& A, CMatrix& B,
                         const std::complex<double>& b, const LongInt& rowOffset=0, const LongInt& colOffset=0);

    CMatrix&  addMatrix (const AMatrix& A, const LongInt& rowIndex=0, const LongInt& colIndex=0);
    
    CMatrix&  addScalMatrix                       (const std::complex<double>& z, const AMatrix& B, const LongInt& rowIndex=0, const LongInt& colIndex=0);
    CMatrix&  addIdentityScaledBy                 (const std::complex<double>& z);
    bool      partialEvaluateTransposeAt          (const LongInt& i, CVector& w, const LongInt& j, const CVector& v) const;
    bool      partialEvaluateConjugateTransposeAt (const LongInt& i, CVector& w, const LongInt& j, const CVector& v) const;

    CMatrix&  isConversionOf (const AMatrix& A);

    CMatrix&  attachColumns (const CMatrix& A);
    CMatrix&  attachRows    (const CMatrix& A);
    CMatrix&  invert        ();
    CMatrix&  isInverseOf   (const AMatrix& s, AMatrix& w);
    LongReal  frobeniusNorm () const;
    char*     correctChar   () const;

  //@}

  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{

  private:

    friend ostream& operator << (ostream & s,const CMatrix &A);
    friend istream& operator >> (istream & s, CMatrix &A);
    friend CMatrix  operator *  (const std::complex<double> alpha, const CMatrix& A);
    friend CMatrix  operator /  (const CMatrix& A, const std::complex<double> alpha);
    friend std::complex<double>  BiLform     (const CVector &v1,const CMatrix& A,const CVector &v2);


    //! Zeiger auf CMatrixelemente
    std::complex<double>* _pelm;

  //@}
 };

#endif // not defined __CMatrix__
