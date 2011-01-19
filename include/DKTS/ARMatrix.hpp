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
 * \class    ARMatrix
 *
 * \brief    ARMatrix aus IR(n,n)
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

#ifndef __ARMatrix__
#define __ARMatrix__

#include "Macros.h"
#include <iostream>
#include "RVector.hpp"
#include "AMatrix.hpp"


class  ARMatrix : public AMatrix
 {
    DECLARE (ARMatrix)
   
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
    //! Konstruktor
    /*! Dies ist der Standardkonstruktor
    */
    ARMatrix (const LongInt& m, const LongInt& n);

    //! Destruktor
    virtual ~ARMatrix ();

  //@}

  /*!
  ***************************************************************************************
  * \name                             �berladene Operatoren
  ***************************************************************************************/

  //@{

  public:
   // virtual RVector    operator *  (const RVector& v) const =0;    
   // virtual LongReal&  operator () (const LongInt i, const LongInt k) const =0;
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

   // friend ARMatrix& operator * (LongReal alpha, const ARMatrix& A);

  private:
    //! 

  //@}
 };

#endif // not defined __ARMatrix__
