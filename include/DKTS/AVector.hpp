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
 * \class    AVector
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

#ifndef __AVector__
#define __AVector__

#include "Macros.h"
#include "SimpleException.hpp"

enum VectorType
 {	
    isAVector,
    isRVector,
	   isCVector
 };

class AVector
 {
    friend class CVector;
    friend class RVector;
    friend class HMatrixLTest;
    
    DECLARE (AVector)

    CRITICAL_ATTRIBUTE (LongInt,    dimension, setDimension)
    CRITICAL_ATTRIBUTE (VectorType, type,      setType)

    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
   
     AVector (const LongInt& m);
    ~AVector ();
  

  //@}

  /*!
  ***************************************************************************************
  * \name                             �berladene Operatoren
  ***************************************************************************************/

  //@{

  public:

     AVector& operator += (const AVector* A);

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


  //@}
 };


typedef  AVector* AVectorPointer;

#endif // not defined __AVector__
