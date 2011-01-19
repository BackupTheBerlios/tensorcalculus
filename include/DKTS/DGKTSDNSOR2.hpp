/*
 * Copyright (C) 2011 Mike Espig
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
 * \class    DGKTSDNSOR2
 *
 * \brief    DGKTSDNSOR2 
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

#ifndef __DGKTSDNSOR2__
#define __DGKTSDNSOR2__


#include "Macros.h"
#include <iostream>
#include <iomanip>

#include "DKTS.hpp"
#include "DGKTSDDataBaseNSOR.hpp"
#include "DKTSDIterationInfo.hpp"
#include "DKTSDDataBlock.hpp"


using namespace std;

class  DGKTSDNSOR2 : public DGKTSDDataBaseNSOR
 {
    DECLARE (DGKTSDNSOR2) 
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
   
     DGKTSDNSOR2 ();

     virtual ~DGKTSDNSOR2();
  

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

     DKTSDIterationInfo truncate2Eps    (const LongReal& eps, const DKTS& A, DKTS& X, const bool& plotInfo=true, const bool& useX=false);
     DKTSDIterationInfo startIteration  (DKTS& x, const DKTS& a);
     DKTSDIterationInfo decompose       (DKTS& x, const DKTS& a, const LongReal& normA, DKTSDDataBlock& dataBlock);
  
     LongReal f                         () const;     
          
     DGKTSDNSOR2& computeAllInnerProducts  (const DKTS& x, const DKTS& a);
     DGKTSDNSOR2& setDirection             (const DKTS& d);
     
     DGKTSDNSOR2& computeGradient          (const DKTS& x, const DKTS& a);
     
     
     DGKTSDNSOR2& computeAu (const LongReal& t);         
     DGKTSDNSOR2& computeAo (const LongReal& t);
     DGKTSDNSOR2& computeXu (const LongReal& t);         
     DGKTSDNSOR2& computeXo (const LongReal& t);

     DGKTSDNSOR2& computeMaxPivot (LongInt& j1, LongInt& mu1);
     
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


typedef  DGKTSDNSOR2* DGKTSDNSOR2Pointer;

#endif // not defined __DGKTSDNSOR2__
