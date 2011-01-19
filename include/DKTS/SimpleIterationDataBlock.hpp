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
 * \class    SimpleIterationDataBlock
 *
 * \brief    SimpleIterationDataBlock 
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

#ifndef __SimpleIterationDataBlock__
#define __SimpleIterationDataBlock__

#include "Macros.h"
#include <iostream>

#include "IterationDataBlock.hpp"
#include "SimpleIterationData.hpp"


using namespace std;

class  SimpleIterationDataBlock : public IterationDataBlock
 {
    DECLARE (SimpleIterationDataBlock)      
    
    CONTAINS  (SimpleIterationData, Entry)
    
				
				
				/*! Strings fuer plot */
    ATTRIBUTE (IString, stepString,  setStepString)
    ATTRIBUTE (IString, errorString, setErrorString)

				/*! Strings fuer TeX-plot */
    ATTRIBUTE (IString, stepTeXString,  setStepTeXString)
    ATTRIBUTE (IString, errorTeXString, setErrorTeXString)
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
  
     SimpleIterationDataBlock();
     SimpleIterationDataBlock(const SimpleIterationDataBlock& data);
     
     virtual ~SimpleIterationDataBlock();

  //@}

  /*!
  ***************************************************************************************
  * \name                             �ueberladene Operatoren
  ***************************************************************************************/

  //@{

  public:

     SimpleIterationDataBlock& operator =  (const SimpleIterationDataBlock& data);
     friend ostream&           operator << (ostream& s, const SimpleIterationDataBlock& sD);

  //@}

  /*!
  ***************************************************************************************
  * \name                        oeffentliche Klassendienste
  ***************************************************************************************/

  //@{

  public:
     
     /*! Ausgabe der Strings auf dem uebergebenen ostream             */
     virtual bool    plot               (ostream& log,              const ProtocolProperties& pP) const;
     
     /*! Ausgabe der Strings in eine .dat Datei (leere Methode)       */
     virtual bool    generatedPlotFiles (const IString& blockIndex, const ProtocolProperties& pP) const;
     
     /*! Ausgabe der Strings in eine .tex Datei                       */
     virtual IString generatedTeXFiles  (const IString& blockIndex, const ProtocolProperties& pP) const;
					
     /*! Zusammenfassung von generatedPlotFiles und generatedTeXFiles */
     virtual bool    generatedAll       (const IString& blockIndex, const ProtocolProperties& pP) const;
     
  //@}


  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{


  private:
  
     /*! Setzt die default Werte                                      */
     SimpleIterationDataBlock setDefaultParameter();
     
     /*! Hilfsmethode fuer plot                                       */
     bool plotEntry       (ostream& log,                  const ProtocolProperties& pP)                                        const;
     
     /*! Hilfsmethode zum zentrieren der einzelnen Eintraege in Abhaengigkeit der Laenge der Zahlen und des Tabellenkopfes */  
     bool centeringEntrys (const LongInt& lengthOfNumber, const LongInt& lengthString, LongInt& disString, LongInt& disNumber) const;
     
     /*! Hilfsmethode fuer generatedPlotFiles                         */
     bool plot2File       (const IString& blockIndex,     const ProtocolProperties& pP)                                        const;
     
     /*! Erstellt Dateinamen blockIndex-name-dateiEndung              */
     IString fileName     (const IString& name,           const IString& blockIndex, const IString& dateiEndung)               const;
     
     /*! Erstellt Dateinamen datum-zeit-blockIndex-name-dateiEndung   */
     IString timeFileName (const IString& name,           const IString& blockIndex, const IString& dateiEndung)               const;

    //! 

  //@}
 };


typedef  SimpleIterationDataBlock* SimpleIterationDataBlockPointer;

#endif // not defined __SimpleIterationDataBlock__
