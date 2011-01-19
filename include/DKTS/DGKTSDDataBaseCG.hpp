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
 * \class    DGKTSDDataBaseCG
 *
 * \brief    DGKTSDDataBaseCG 
 *
 *
 * \author   Henry Auer, Mike Espig
 *
 * \ingroup  KTSD
 *
 *
 ******************************************************************************
 *
 *
 ******************************************************************************
 */

#ifndef __DGKTSDDataBaseCG__
#define __DGKTSDDataBaseCG__


#include "Macros.h"
#include <iostream>
#include <iomanip>

#include "DKTS.hpp"

#include "DKTSDSimpleInnerProducts.hpp"

using namespace std;

class  DGKTSDDataBaseCG
 {
    DECLARE (DGKTSDDataBaseCG) 
 
    friend class DGKTSDCG;
 
    /*!
     The Order of the Tensor-Sum*/
    ATTRIBUTE (LongInt, d, setD)

    /*!
     The small Tensor-Rank of the Tensor-Sum*/
    ATTRIBUTE (LongInt, r, setR)

    
    /*!
    The Tensor-Rank of a*/
    ATTRIBUTE (LongInt, k, setK)
    
    ATTRIBUTE (LongInt, m, setM)
    

    ATTRIBUTE (LongReal, normA2, setNormA2)
    
    ATTRIBUTE (LongReal, epsilon,         setEpsilon)
    ATTRIBUTE (LongReal, error,           setError)    
    ATTRIBUTE (LongReal, accuracy,        setAccuracy)
    
    ATTRIBUTE (LongReal, alpha,           setAlpha)
    
    ATTRIBUTE (LongReal, C, setC)
    ATTRIBUTE (LongReal, D, setDD)
    ATTRIBUTE (LongReal, epsilonStepSize, setEpsilonStepSize)    
    ATTRIBUTE (LongReal, scalingFactor,   setScalingFactor)
    //edit
    ATTRIBUTE (LongReal, scalingFactorOld, setScalingFactorOld) 
    
    ATTRIBUTE (LongInt,  maxIterations,   setMaxIterations)
    ATTRIBUTE (LongInt,  maxIterationsCR, setMaxIterationsCR)
        
    ATTRIBUTE (bool,     printCout,     setPrintCout)
    ATTRIBUTE (bool,     printLog,      setPrintLog)
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
   
     DGKTSDDataBaseCG ();

     virtual ~DGKTSDDataBaseCG();
  

  //@}

  /*!
  ***************************************************************************************
  * \name                             �berladene Operatoren
  ***************************************************************************************/

  //@{

  public:
  
     LongReal& Au(const LongInt& i, const LongInt& j1, const LongInt& mu1) const
      {
         // speicheranordnung: mu1(j1(i))
         const LongInt index = i + (j1 + mu1*attr_r)*attr_k;
                            
        return (attr_valuesAu[index]);                
      }

     LongReal& Ao(const LongInt& i, const LongInt& j1, const LongInt& mu1) const
      {
         const LongInt index = i + (j1 + mu1*attr_r)*attr_k;
                            
        return (attr_valuesAo[index]);                
      }
     LongReal& Xu(const LongInt& i, const LongInt& j1, const LongInt& mu1) const
      {
         const LongInt index = i + (j1 + mu1*attr_r)*attr_r;
                            
        return (attr_valuesXu[index]);                
      }

     LongReal& Xo(const LongInt& i, const LongInt& j1, const LongInt& mu1) const
      {
         const LongInt index = i + (j1 + mu1*attr_r)*attr_r;
                            
        return (attr_valuesXo[index]);                
      }
     
  //@}

  /*!
  ***************************************************************************************
  * \name                        �ffentliche Klassendienste
  ***************************************************************************************/

  //@{

  public:

     LongReal& Ainv(const LongInt& j1, const LongInt& j2, const LongInt& mu) const
      {
         const LongInt index = (MAX(j1,j2) + mu*attr_r)*attr_r + MIN(j1,j2);
 
        return attr_inversOfA[index];
      }


     bool resize (const DKTS& x, const DKTS& a);
     
     
  //@}


  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{


    DKTS attr_gradient;
    DKTS attr_gradientOld;
    DKTS attr_direction;
    
    DKTS attr_work;
    
    DKTS attr_gradientS;        
    
    DKTSDSimpleInnerProducts attr_xa;
    DKTSDSimpleInnerProducts attr_da;
    
    DKTSDSimpleInnerProducts attr_xx;
    DKTSDSimpleInnerProducts attr_xd;
    DKTSDSimpleInnerProducts attr_dd;
    
    RVector attr_v1;
    RVector attr_v2;
   

  private:

    DGKTSDDataBaseCG&    allocateWorkSpace (const LongInt& d, const LongInt& k, const LongInt& r, const LongInt& m);
    DGKTSDDataBaseCG&    deleteWorkSpace   ();    

           
    DGKTSDDataBaseCG& setDefaultParameter  (const LongReal& normA);

    LongRealPointer attr_inversOfA;
    
    LongRealPointer attr_valuesAu;    
    LongRealPointer attr_valuesAo;
    LongRealPointer attr_valuesXu;    
    LongRealPointer attr_valuesXo;
    
    //!

  //@}
 };


typedef  DGKTSDDataBaseCG* DGKTSDDataBaseCGPointer;

#endif // not defined __DGKTSDDataBaseCG__
