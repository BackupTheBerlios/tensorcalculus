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
 * \class    LongRealSvd
 *
 * \brief    LongRealSvd
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

#ifndef __LongRealSvd__
#define __LongRealSvd__

#include "Macros.h"
#include <iostream>
#include "RMatrix.hpp"
#include "LapackInterface.hpp"

using namespace std;

class LongRealSvd : public LapackInterface
 {
    DECLARE (LongRealSvd)

    friend class RkRMatrix;
    friend class KTS2;

    static const char   job[];

    CRITICAL_ATTRIBUTE (RMatrix, U,     setU)
    CRITICAL_ATTRIBUTE (RMatrix, V,     setV)
    CRITICAL_ATTRIBUTE (RVector, sigma, setSigma)

  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
    //! Konstruktor
    /*! Dies ist der Standardkonstruktor
    \param n die Anzahl der Spalten und Zeilen*/
    LongRealSvd (const RMatrix& A, const LongInt k);

    LongRealSvd (const RMatrix& A);


    //! Destruktor
    virtual ~LongRealSvd ();

  //@}

  /*!
  ***************************************************************************************
  * \name                             �berladene Operatoren
  ***************************************************************************************/

  //@{

  public:

  //@}

  /*!
  ***************************************************************************************
  * \name                        �ffentliche Klassendienste
  ***************************************************************************************/

  //@{

  public:

    LongRealSvd&  multiplyUSigma ();


  //@}


  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{


  private:

    LongRealSvd&  svd         (const RMatrix& A);
    LongRealSvd&  generateUby (const RMatrix& U);
    LongRealSvd&  generateVby (const RMatrix& U);

    //! 

  //@}
 };

#endif // not defined __LongRealSvd__
