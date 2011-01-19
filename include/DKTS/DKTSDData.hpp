/*
 * Copyright (C) Marcel Schuster
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
 * \class    DKTSDData
 *
 * \brief    DKTSDData 
 *
 *
 * \author   Marcel Schuster
 *
 * \ingroup  
 *
 *
 ******************************************************************************
 *
 *
 ******************************************************************************
 */

#ifndef __DKTSDData__
#define __DKTSDData__

#include "Macros.h"
#include <iostream>

#include "SimpleIterationData.hpp"

using namespace std;

class  DKTSDData : public SimpleIterationData
 {
    DECLARE (DKTSDData)
				      
				/*! dienen zum setzen der einzelnen Werte */
    ATTRIBUTE (LongReal, normOfGradient,   setNormOfGradient)
    ATTRIBUTE (LongReal, delta,            setDelta)
    ATTRIBUTE (LongReal, stepSize,	        setStepSize)
				ATTRIBUTE (LongReal, relFunctionValue, setRelFunctionValue)
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
  
     DKTSDData(const LongInt&  st =0, const LongReal& er =0,
					          const LongReal& nOG=0, const LongReal& del=0, const LongReal& sS=0, const LongReal& rFV=0);
     DKTSDData(const DKTSDData& dD);
     
     virtual ~DKTSDData();

  //@}

  /*!
  ***************************************************************************************
  * \name                             �ueberladene Operatoren
  ***************************************************************************************/

  //@{

  public:

     virtual bool    operator == (const ProtocolData& D);
     DKTSDData&      operator =  (const DKTSDData& data);
     friend ostream& operator << (ostream& s, const DKTSDData& data);
     friend istream& operator >> (istream& s, DKTSDData& data);

  //@}

  /*!
  ***************************************************************************************
  * \name                        oeffentliche Klassendienste
  ***************************************************************************************/

  //@{

  public:

     /*! Ausgabe der Strings auf dem uebergebenen ostream */
     virtual bool plot                  (ostream& log,              const ProtocolProperties& pP) const;
     
     /*! rein virtuelle Methoden                          */
     virtual bool    generatedAll       (const IString& blockIndex, const ProtocolProperties& pP) const {return false;}
     virtual bool    generatedPlotFiles (const IString& blockIndex, const ProtocolProperties& pP) const {return false;}
     virtual IString generatedTeXFiles  (const IString& blockIndex, const ProtocolProperties& pP) const {return false;}

  //@}


  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{


  private:
  
     DKTSDData setDefaultParameter();
  
    //! 

  //@}
 };


typedef  DKTSDData* DKTSDDataPointer;

#endif // not defined __DKTSDData__
