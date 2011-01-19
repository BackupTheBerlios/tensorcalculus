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
 * \class    KTVector
 *
 * \brief    KTVector 
 *
 *
 * \author   Mike Espig
 *
 * \ingroup  KTSD
 *
 *
 ******************************************************************************
 *
 *
 ******************************************************************************
 */

#ifndef __KTVector__
#define __KTVector__

#include "Macros.h"
#include <iostream>


using namespace std;

class  KTVector
 {
    DECLARE (KTVector)
   
    friend class DKTVector;
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:

     KTVector (LongIntPointer n);
     
     virtual ~KTVector ();


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
  
     LongInt& n() const
              {
                return *attr_n;
              }

  //@}


  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{


  private:
  
     LongIntPointer attr_n;
  
    //! 

  //@}
 };


typedef  KTVector* KTVectorPointer;

#endif // not defined __KTVector__
