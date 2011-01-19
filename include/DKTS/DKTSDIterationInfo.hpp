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
 * \class    DKTSDIterationInfo
 *
 * \brief    DKTSDIterationInfo 
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

#ifndef __DKTSDIterationInfo__
#define __DKTSDIterationInfo__

#include "Macros.h"
#include <iostream>

#include "SimpleIterationData.hpp"

using namespace std;

class  DKTSDIterationInfo : public SimpleIterationData
 {
    DECLARE (DKTSDIterationInfo)
    
    ATTRIBUTE (LongReal, startError,          setStartError)
    ATTRIBUTE (LongReal, numberOfNewtonSteps, setNumberOfNewtonSteps)
    ATTRIBUTE (LongReal, calculationTime,     setCalculationTime)
    ATTRIBUTE (LongReal, relativeDifferenz,   setRelativeDifferenz)
				ATTRIBUTE (LongReal, gradient,            setGradient)
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
  
     DKTSDIterationInfo (const LongInt&  st=0,   const LongReal& er=0.0, const LongReal& ser =0.0, 
                         const LongReal& ra=0.0, const LongReal& gr=0.0, const LongReal& nons=0.0, 
                         const LongReal& ct=0.0);
     DKTSDIterationInfo (const DKTSDIterationInfo& data);
     
     virtual ~DKTSDIterationInfo();

  //@}

  /*!
  ***************************************************************************************
  * \name                             �ueberladene Operatoren
  ***************************************************************************************/

  //@{

  public:
     
     virtual bool                operator == (const ProtocolData& data);
             DKTSDIterationInfo& operator =  (const DKTSDIterationInfo& data);
     friend  ostream&            operator << (ostream& s, const DKTSDIterationInfo& data);
					friend  istream&            operator >> (istream& s, DKTSDIterationInfo& data);

  //@}

  /*!
  ***************************************************************************************
  * \name                        oeffentliche Klassendienste
  ***************************************************************************************/

  //@{

  public:
  
     /*! Ausgabe der Strings auf dem uebergebenen ostream             */
     virtual bool    plot               (ostream& log,              const ProtocolProperties& pP) const;
     
     /*! rein virtuelle Methoden                                      */
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
     
     DKTSDIterationInfo setDefaultParameter();
  
  
    //! 

  //@}
 };


typedef  DKTSDIterationInfo* DKTSDIterationInfoPointer;

#endif // not defined __DKTSDIterationInfo__
