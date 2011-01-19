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
 * \class    SimpleIterationData
 *
 * \brief    SimpleIterationData
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

#ifndef __SimpleIterationData__
#define __SimpleIterationData__

#include "Macros.h"
#include <iostream>

#include "ProtocolData.hpp"

using namespace std;


class  SimpleIterationData : public ProtocolData
 {
    DECLARE (SimpleIterationData)      
    
				/*! dienen zum setzen der einzelnen Werte */
    ATTRIBUTE (LongInt,  step,  setStep)
    ATTRIBUTE (LongReal, error, setError) 
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
  
     SimpleIterationData(const LongInt& st=0, const LongReal& er=0);
     SimpleIterationData(const SimpleIterationData& pD);
     
     virtual ~SimpleIterationData();

  //@}

  /*!
  ***************************************************************************************
  * \name                             �ueberladene Operatoren
  ***************************************************************************************/

  //@{

  public:

     virtual bool         operator == (const ProtocolData& D);
     SimpleIterationData& operator =  (const SimpleIterationData& sD);
					friend  ostream&     operator << (ostream& s, const SimpleIterationData& sD);
					friend  istream&     operator >> (istream& s, SimpleIterationData& data);

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
     
     SimpleIterationData setDefaultParameter();
  
    //! 

  //@}
 };


typedef  SimpleIterationData* SimpleIterationDataPointer;

#endif // not defined __SimpleIterationData__
