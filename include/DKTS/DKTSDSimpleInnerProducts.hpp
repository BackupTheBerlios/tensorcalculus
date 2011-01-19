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
 * \class    DKTSDSimpleInnerProducts
 *
 * \brief    DKTSDSimpleInnerProducts 
 *
 *
 * \author   Mike Espig, Henry Auer
 *
 * \ingroup  
 *
 *
 ******************************************************************************
 *
 *
 ******************************************************************************
 */

#ifndef __DKTSDSimpleInnerProducts__
#define __DKTSDSimpleInnerProducts__

#include "Macros.h"
#include <iostream>

#include "DKTS.hpp"

using namespace std;

class  DKTSDSimpleInnerProducts
 {
    DECLARE (DKTSDSimpleInnerProducts)    
    
    ATTRIBUTE (LongInt, d,  setD)
    ATTRIBUTE (LongInt, r1, setR1)
    ATTRIBUTE (LongInt, r2, setR2)
    
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
  
     DKTSDSimpleInnerProducts (const LongInt& d=0, const LongInt& r1=0, const LongInt& r2=0);
     DKTSDSimpleInnerProducts (const DKTS& x, const DKTS& y);
     
     virtual ~DKTSDSimpleInnerProducts();

  //@}

  /*!
  ***************************************************************************************
  * \name                             �berladene Operatoren
  ***************************************************************************************/

  //@{

  public:
  
     DKTSDSimpleInnerProducts& operator *= (const LongReal& alpha);
     
     LongReal& innerProductAt(const LongInt& i1, const LongInt& i2, const LongInt& mu) const
                {
                   const LongInt index = attr_d*(attr_r2*i1 + i2) + mu;
                                      
                  return (attr_valuesS[index]);                
                }
               
     LongReal& operator () (const LongInt& i1, const LongInt& mu) const
                {                
                  return (tensorInnerProductAt(i1, 0, mu));
                }

                
     LongReal& tensorInnerProductAt(const LongInt& i1, const LongInt& i2, const LongInt& mu) const
                {
                   const LongInt index = i2 + (i1 + mu*attr_r1)*attr_r2;
                                      
                  return (attr_valuesG[index]);                
                }

     LongReal& aU(const LongInt& i1, const LongInt& i2, const LongInt& mu) const
                {
                   const LongInt index = i2 + (i1 + mu*attr_r1)*attr_r2;
                                      
                  return (attr_valuesAu[index]);                
                }

     LongReal& aO(const LongInt& i1, const LongInt& i2, const LongInt& mu) const
                {
                   const LongInt index = i2 + (i1 + mu*attr_r1)*attr_r2;
                                      
                  return (attr_valuesAo[index]);                
                }

                
  //@}

  /*!
  ***************************************************************************************
  * \name                        �ffentliche Klassendienste
  ***************************************************************************************/

  //@{

  public:

     DKTSDSimpleInnerProducts& resize       (const DKTS& x, const DKTS& y);
     DKTSDSimpleInnerProducts& resize       (const LongInt& d, const LongInt& r1, const LongInt& r2);
     DKTSDSimpleInnerProducts& setAllEntrys ();

     DKTSDSimpleInnerProducts& computeAllInnerProducts    (const DKTS& x, const DKTS& y);
     DKTSDSimpleInnerProducts& computeInnerProducts       (const DKTS& x, const DKTS& y);
     DKTSDSimpleInnerProducts& computeSimpleInnerProducts (); 
     
     //Henry
     DKTSDSimpleInnerProducts& computeaU ();         
     DKTSDSimpleInnerProducts& computeaO ();
     DKTSDSimpleInnerProducts& computeSimpleLeaveOut ();

  //@}


  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{


  private:

     DKTSDSimpleInnerProducts& allocateDataSpace (const LongInt& d, const LongInt& r1, const LongInt& r2);
     DKTSDSimpleInnerProducts& deleteDataSpace   ();

     LongRealPointer attr_valuesS;     //innerProductAt
     LongRealPointer attr_valuesG;     //tensorInnerProductAt

     LongRealPointer attr_valuesAu;    
     LongRealPointer attr_valuesAo;
     

  //@}
 };


typedef  DKTSDSimpleInnerProducts* DKTSDSimpleInnerProductsPointer;

#endif // not defined __DKTSDSimpleInnerProducts__
