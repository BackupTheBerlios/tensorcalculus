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
 * \defgroup LapackInterface
 * \brief Diese Bibliothek enth�lt grundlegende Dienstklassen f�r 
 *   Lapack- Algorithmen
 ******************************************************************************
*/

/*!
 ******************************************************************************
 * \class    LapackInterface
 *
 * \brief    LapackInterface 
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

#ifndef __LapackInterface__
#define __LapackInterface__

#include "Macros.h"
#include <iostream>
// #include <mkl.h>

using namespace std;

//! \todo Eigentlich sollte man das LapackInterface neu modelieren, ist im jetztigen Zusstand "quick & dirty".


class  LapackInterface
 {
    DECLARE (LapackInterface)
   
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
   
  

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


typedef  LapackInterface* LapackInterfacePointer;

#endif // not defined __LapackInterface__
