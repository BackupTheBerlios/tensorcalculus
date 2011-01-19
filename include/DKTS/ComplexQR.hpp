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
 * \class    ComplexQR
 *
 * \brief    ComplexQR aus IC(n,n)
 *
 *
 * \author   Mike Espig
 *
 * \ingroup  LapackInterface
 *
 *
 ******************************************************************************
 *
 *
 ******************************************************************************
 */

#ifndef __ComplexQR__
#define __ComplexQR__

#include "Macros.h"
#include <iostream>
#include "CMatrix.hpp"
#include "LapackInterface.hpp"

using namespace std;

class ComplexQR : public LapackInterface
 {
    DECLARE (ComplexQR)

    friend class RkCMatrix;

    CRITICAL_ATTRIBUTE (CMatrix, Q,     setQ)
    CRITICAL_ATTRIBUTE (CMatrix, R,     setR)
    CRITICAL_ATTRIBUTE (LongInt, order, setOrder)

  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
    //! Konstruktor
    /*! Dies ist der Standardkonstruktor
    \param n die Anzahl der Spalten und Zeilen*/
    ComplexQR (const CMatrix& A, const LongInt& k);


    //! Destruktor
    virtual ~ComplexQR ();

  //@}

  /*!
  ***************************************************************************************
  * \name                             �berladene Operatoren
  ***************************************************************************************/

  //@{

  public:

   friend ostream& operator << (ostream&, const ComplexQR & qr);
   friend istream& operator >> (istream&, ComplexQR &qr);

  //@}

  /*!
  ***************************************************************************************
  * \name                        �ffentliche Klassendienste
  ***************************************************************************************/

  //@{

  public:

  //@}


  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{


  private:

    ComplexQR&  qr          (const CMatrix& A);
    ComplexQR&  generateQby (const CMatrix& A);
    ComplexQR&  generateRby (const CMatrix& A);

    //! 

  //@}
 };

#endif // not defined __ComplexQR__
