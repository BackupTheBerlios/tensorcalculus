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
 * \class    DKTSTruncation2Eps
 *
 * \brief    DKTSTruncation2Eps 
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

#ifndef __DKTSTruncation2Eps__
#define __DKTSTruncation2Eps__

#include "Macros.h"
#include <iostream>

#include "DKTSTruncation.hpp"

using namespace std;

class  DKTSTruncation2Eps : public DKTSTruncation
 {
    DECLARE (DKTSTruncation2Eps)
   
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
   
     DKTSTruncation2Eps (const LongInt& d=1, const LongInt& R=1, const LongInt& n=1);
     DKTSTruncation2Eps (const IString& fileName, const LongInt& rMin);
     DKTSTruncation2Eps (const IString& fileName, const IString& coFileName, const IString& ortFileName, const LongInt& rMin);
     
     DKTSTruncation2Eps (const DKTS& A);
     
     virtual ~DKTSTruncation2Eps();

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

     DKTSTruncation2Eps& truncateF2                 (const LongReal& eps, const LongReal& epsN, const LongReal& preC);
     DKTSTruncation2Eps& truncate2                  (const LongReal& eps, const DKTS& A, DKTS& X, const bool& plotInfo=true, const bool& useX=false);
     DKTSTruncation2Eps& truncateSuccessiveAls2     (const LongReal& eps, const DKTS& A, DKTS& X, const bool& plotInfo=true);
     LongReal            truncate2EstimatedAccuracy (const DKTS& A, DKTS& X, const bool& plotInfo=true, const LongInt& st=1);          
          
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


typedef  DKTSTruncation2Eps* DKTSTruncation2EpsPointer;

#endif // not defined __DKTSTruncation2Eps__
