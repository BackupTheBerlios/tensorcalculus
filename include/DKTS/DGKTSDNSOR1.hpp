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
 * \class    DGKTSDNSOR1
 *
 * \brief    DGKTSDNSOR1 
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

#ifndef __DGKTSDNSOR1__
#define __DGKTSDNSOR1__


#include "Macros.h"
#include <iostream>
#include <iomanip>

#include "DKTS.hpp"
#include "DGKTSDDataBaseNSOR.hpp"
#include "DKTSDIterationInfo.hpp"
#include "DKTSDDataBlock.hpp"


using namespace std;

class  DGKTSDNSOR1 : public DGKTSDDataBaseNSOR
 {
    DECLARE (DGKTSDNSOR1) 
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
   
     DGKTSDNSOR1 ();

     virtual ~DGKTSDNSOR1();
  

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
          
     DGKTSDNSOR1& computeAllInnerProducts  (const DKTS& x, const DKTS& a);
     DGKTSDNSOR1& setDirection             (const DKTS& d);
     
     DGKTSDNSOR1& computeGradient          (const DKTS& x, const DKTS& a);
     
     
     DGKTSDNSOR1& computeAu (const LongReal& t);         
     DGKTSDNSOR1& computeAo (const LongReal& t);
     DGKTSDNSOR1& computeXu (const LongReal& t);         
     DGKTSDNSOR1& computeXo (const LongReal& t);

     DGKTSDNSOR1& computeMaxPivot (LongInt& j1, LongInt& mu1);
     
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


typedef  DGKTSDNSOR1* DGKTSDNSOR1Pointer;

#endif // not defined __DGKTSDNSOR1__
