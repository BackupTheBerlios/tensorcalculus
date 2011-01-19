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
 * \class    DGKTSDDataBaseNSOR
 *
 * \brief    DGKTSDDataBaseNSOR 
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

#ifndef __DGKTSDDataBaseNSOR__
#define __DGKTSDDataBaseNSOR__


#include "Macros.h"
#include <iostream>
#include <iomanip>

#include "DKTS.hpp"

#include "DKTSDSimpleInnerProducts.hpp"

using namespace std;

class  DGKTSDDataBaseNSOR
 {
    DECLARE (DGKTSDDataBaseNSOR) 
 
    friend class DGKTSDNSOR1;
    friend class DGKTSDNSOR2;
    friend class DGKTSDNSOR3;
 
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
    
    ATTRIBUTE (LongInt,  maxIterations,   setMaxIterations)
        
    ATTRIBUTE (bool,     printCout,     setPrintCout)
    ATTRIBUTE (bool,     printLog,      setPrintLog)
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
   
     DGKTSDDataBaseNSOR ();

     virtual ~DGKTSDDataBaseNSOR();
  

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


     bool resize (const DKTS& x, const DKTS& a);
     
     
  //@}


  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{


    DKTS attr_gradient;
    
    DKTSDSimpleInnerProducts attr_xa;
    DKTSDSimpleInnerProducts attr_xx;
    
   

  private:

    DGKTSDDataBaseNSOR&    allocateWorkSpace (const LongInt& d, const LongInt& k, const LongInt& r, const LongInt& m);
    DGKTSDDataBaseNSOR&    deleteWorkSpace   ();    

           
    DGKTSDDataBaseNSOR& setDefaultParameter  (const LongReal& normA);

    LongRealPointer attr_inversOfA;
    
    LongRealPointer attr_valuesAu;    
    LongRealPointer attr_valuesAo;
    LongRealPointer attr_valuesXu;    
    LongRealPointer attr_valuesXo;
    
    //!

  //@}
 };


typedef  DGKTSDDataBaseNSOR* DGKTSDDataBaseNSORPointer;

#endif // not defined __DGKTSDDataBaseNSOR__
