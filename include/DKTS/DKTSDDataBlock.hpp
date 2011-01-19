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
 * \class    DKTSDDataBlock
 *
 * \brief    DKTSDDataBlock 
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

#ifndef __DKTSDDataBlock__
#define __DKTSDDataBlock__

#include "Macros.h"
#include <iostream>

#include "IterationDataBlock.hpp"
#include "DKTSDData.hpp"


using namespace std;

class  DKTSDDataBlock : public IterationDataBlock
 {
    DECLARE (DKTSDDataBlock)      
    
    CONTAINS (DKTSDData, Entry)
    
				/*! Strings fuer plot */
    ATTRIBUTE (IString, stepString,             setStepString)
    ATTRIBUTE (IString, normOfGradientString,   setNormOfGradientString)
    ATTRIBUTE (IString, deltaString,            setDeltaString)
    ATTRIBUTE (IString, stepSizeString,         setStepSizeString)
    ATTRIBUTE (IString, errorString,            setErrorString)
				ATTRIBUTE (IString, relFunctionValueString, setRelFunctionValueString)
    
				/*! Strings fuer TeX-plot */
    ATTRIBUTE (IString, stepTeXString,             setStepTeXString)
    ATTRIBUTE (IString, normOfGradientTeXString,   setNormOfGradientTeXString)
    ATTRIBUTE (IString, deltaTeXString,            setDeltaTeXString)
    ATTRIBUTE (IString, stepSizeTeXString,         setStepSizeTeXString)
    ATTRIBUTE (IString, errorTeXString,            setErrorTeXString)
				ATTRIBUTE (IString, relFunctionValueTeXString, setRelFunctionValueTeXString)
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
  
     DKTSDDataBlock();
     DKTSDDataBlock(const DKTSDDataBlock& data);

     virtual ~DKTSDDataBlock();

  //@}

  /*!
  ***************************************************************************************
  * \name                             �ueberladene Operatoren
  ***************************************************************************************/

  //@{

  public:

     DKTSDDataBlock& operator =  (const DKTSDDataBlock& data);
					friend ostream& operator << (ostream& s, const DKTSDDataBlock& dD);

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
     DKTSDDataBlock setDefaultParameter();
     
     /*! Hilfsmethode fuer plot                                       */
     bool plotEntry       (ostream& log,                  const ProtocolProperties& pP)                                        const;
     
     /*! Hilfsmethode zum zentrieren der einzelnen Eintraege in Abhaengigkeit der Laenge der Zahlen und des Tabellenkopfes  */
     bool centeringEntrys (const LongInt& lengthOfNumber, const LongInt& lengthString, LongInt& disString, LongInt& disNumber) const;
     
     /*! Hilfsmethode fuer generatedPlotFiles                         */
     bool plot2File       (const IString& blockIndex,     const ProtocolProperties& pP)                                        const;
     
     /*! Erstellt Dateinamen blockIndex-name-dateiEndung              */
     IString fileName     (const IString& name,           const IString& blockIndex,   const IString& dateiEndung)             const;
     
     /*! Erstellt Dateinamen datum-zeit-blockIndex-name-dateiEndung   */
     IString timeFileName (const IString& name,           const IString& blockIndex,   const IString& dateiEndung)             const;
     
    //! 

  //@}
 };


typedef  DKTSDDataBlock* DKTSDDataBlockPointer;

#endif // not defined __DKTSDDataBlock__
