/*
 * Copyright (C) Mike Espig, Henry Auer
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
 * \class    DGKTSDCG
 *
 * \brief    DGKTSDCG 
 *
 *
 * \author   Henry Auer; Mike Espig
 *
 * \ingroup  KTSD
 *
 *
 ******************************************************************************
 *
 *
 ******************************************************************************
 */

#ifndef __DGKTSDCG__
#define __DGKTSDCG__


#include "Macros.h"
#include <iostream>
#include <iomanip>

#include "DKTS.hpp"
#include "DGKTSDDataBaseCG.hpp"
#include "DKTSDIterationInfo.hpp"
#include "DKTSDDataBlock.hpp"


using namespace std;

class  DGKTSDCG : public DGKTSDDataBaseCG
 {
    DECLARE (DGKTSDCG) 
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
   
     DGKTSDCG ();

     virtual ~DGKTSDCG();
  

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
  
     DGKTSDCG& preComputeValuesForStepSizeRule  (const DKTS& x, const DKTS& a);
     
     LongReal computeStepSize           (const DKTS& x, const DKTS& a);
     LongReal valueFromCurryRule        (const DKTS& x, const DKTS& a);

     LongReal derivativeOfPhi           (const LongReal& t, const DKTS& x, const DKTS& a);
     LongReal f                         () const;
     
     
     
     DKTSDIterationInfo truncate2Eps    (const LongReal& eps, const DKTS& A, DKTS& X, const bool& plotInfo=true, const bool& useX=false);
     DKTSDIterationInfo startIteration  (DKTS& x, const DKTS& a);
     DKTSDIterationInfo decompose       (DKTS& x, const DKTS& a, const LongReal& normA, DKTSDDataBlock& dataBlock);
     
     DGKTSDCG& computeAllInnerProducts  (const DKTS& x, const DKTS& a);
     DGKTSDCG& setDirection             (const DKTS& d);
     
     DGKTSDCG& computeDirection         (const DKTS& x, const DKTS& a);
     DGKTSDCG& computeGradient          (const DKTS& x, const DKTS& a);
     DGKTSDCG& computeAinv              ();     
     
     bool      evaluateA                (const LongReal& alpha, const DKTS& v, DKTS& w) const;
     bool      evaluateAinv             (const LongReal& alpha, const DKTS& v, DKTS& w) const;
     
     DGKTSDCG& computeAu (const LongReal& t);         
     DGKTSDCG& computeAo (const LongReal& t);
     DGKTSDCG& computeXu (const LongReal& t);         
     DGKTSDCG& computeXo (const LongReal& t);
          
     
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


typedef  DGKTSDCG* DGKTSDCGPointer;

#endif // not defined __DGKTSDCG__
