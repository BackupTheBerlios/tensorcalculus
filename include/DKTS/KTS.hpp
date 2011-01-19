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
 * \defgroup KTS
 * \brief Diese Bibliothek enth�lt grundlegende Dienstklassen f�r Summen von Kronecker-Tensor-Produkten
 *   
 ******************************************************************************
*/
/*!
 ******************************************************************************
 * \class    KTS
 *
 * \brief    KTS 
 *
 *
 * \author   Mike Espig
 *
 * \ingroup  KTS
 *
 *
 ******************************************************************************
 *
 *
 ******************************************************************************
 */

#ifndef __KTS__
#define __KTS__

#include "Macros.h"
#include <iostream>

using namespace std;


class  KTS
 {
    DECLARE (KTS)
    
    friend class DKTS;
   
    ATTRIBUTE (LongInt, d, setD)
    ATTRIBUTE (LongInt, k, setK)    
    ATTRIBUTE (LongInt, n, setN)
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:

     KTS(const LongInt& d=1, const LongInt& k=1, const LongInt& n=1);
    
     virtual ~KTS();
   
  

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


typedef  KTS* KTSPointer;

#endif // not defined __KTS__
