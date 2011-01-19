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
 * \class    Random
 *
 * \brief    Random 
 *
 *
 * \author   Mike Espig
 *
 * \ingroup  
 *
 *
 ******************************************************************************
 *
 *
 ******************************************************************************
 */

#ifndef __Random__
#define __Random__

#include "Macros.h"
#include <iostream>


using namespace std;

#include <stdlib.h>
#include <time.h>


class  Random
 {
    DECLARE (Random)      
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
  
     Random ();
     
     virtual ~Random();

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
  
      LongInt  randomLongInt  (const LongInt&  max, const LongInt& min=0) const;
      LongReal randomLongReal (const LongReal& max, const LongReal& min=0.0) const;

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


typedef  Random* RandomPointer;

#endif // not defined __Random__
