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
 * \class    ACMatrix
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

#ifndef __ACMatrix__
#define __ACMatrix__

#include "Macros.h"
#include <iostream>
#include "AMatrix.hpp"
#include "CVector.hpp"

using namespace std;

class  ACMatrix : public AMatrix
 {

    DECLARE (ACMatrix)

  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
    //! Konstruktor
    /*! Dies ist der Standardkonstruktor
    \param m die Anzahl der Spalten 
    \param n die Anzahl der Zeilen 
    */
    ACMatrix (const LongInt& m, const LongInt& n);

    //! Destruktor
    virtual ~ACMatrix ();

  //@}

  /*!
  ***************************************************************************************
  * \name                             �berladene Operatoren
  ***************************************************************************************/

  //@{

  public:
    virtual CVector  operator *  (const CVector& v)=0;    

  //@}

  /*!
  ***************************************************************************************
  * \name                        �ffentliche Klassendienste
  ***************************************************************************************/

  //@{

  public:

     virtual AMatrix& setNullMatrix (){return (*this);}

  //@}


  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{


  private:

    //! 

  //@}
 };

typedef ACMatrix* ACMatrixPointer;
typedef ACMatrix& ACMatrixReference;

#endif // not defined __ACMatrix__
