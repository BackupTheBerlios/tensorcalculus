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
 * \class    LongRealQR
 *
 * \brief    LongRealQR 
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

#ifndef __LongRealQR__
#define __LongRealQR__

#include "Macros.h"
#include <iostream>
#include "RMatrix.hpp"
#include "LapackInterface.hpp"


using namespace std;

class LongRealQR : public LapackInterface
 {
    DECLARE (LongRealQR)

    friend class RkRMatrix;
    friend class KTS2;

    CRITICAL_ATTRIBUTE (RMatrix, Q,     setQ)
    CRITICAL_ATTRIBUTE (RMatrix, R,     setR)
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
    LongRealQR (const RMatrix& A, const LongInt& k);

    LongRealQR (const RMatrix& A);

    //! Destruktor
    virtual ~LongRealQR ();

  //@}

  /*!
  ***************************************************************************************
  * \name                             �berladene Operatoren
  ***************************************************************************************/

  //@{

  public:

   friend ostream& operator << (ostream&, const LongRealQR & qr);
   friend istream& operator >> (istream&, LongRealQR &qr);

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

    LongRealQR&  qr          (const RMatrix& A);
    LongRealQR&  generateQby (const RMatrix& A);
    LongRealQR&  generateRby (const RMatrix& A);

    //! 

  //@}
 };

#endif // not defined __LongRealQR__
