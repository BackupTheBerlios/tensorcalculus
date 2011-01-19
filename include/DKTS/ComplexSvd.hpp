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
 * \class    ComplexSvd
 *
 * \brief    ComplexSvd
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

#ifndef __ComplexSvd__
#define __ComplexSvd__

#include "Macros.h"
#include <iostream>
#include "CMatrix.hpp"
#include "LapackInterface.hpp"

using namespace std;

class ComplexSvd : public LapackInterface
 {
    DECLARE (ComplexSvd)

    friend class RkCMatrix;

    static const char   job[];

    CRITICAL_ATTRIBUTE (CMatrix, U,     setU)
    CRITICAL_ATTRIBUTE (CMatrix, V,     setV)
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
    ComplexSvd (const CMatrix& A, const LongInt k);


    //! Destruktor
    virtual ~ComplexSvd ();

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

    ComplexSvd&  multiplyUSigma ();


  //@}


  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{


  private:

    ComplexSvd&  svd         (const CMatrix& A);
    ComplexSvd&  generateUby (const CMatrix& U);
    ComplexSvd&  generateVby (const CMatrix& U);

    //! 

  //@}
 };

#endif // not defined __ComplexSvd__
