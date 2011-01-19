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
 * \class    DKTSDIterationInfoBlock
 *
 * \brief    DKTSDIterationInfoBlock 
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

#ifndef __DKTSDIterationInfoBlock__
#define __DKTSDIterationInfoBlock__

#include "Macros.h"
#include <iostream>

#include "IterationDataBlock.hpp"
#include "DKTSDIterationInfo.hpp"

using namespace std;

class  DKTSDIterationInfoBlock : public IterationDataBlock
 {
    DECLARE (DKTSDIterationInfoBlock)
    
    CONTAINS (DKTSDIterationInfo, Entry)
    
				/*! Strings fuer plot */
    ATTRIBUTE (IString, stepString,                setStepString)
    ATTRIBUTE (IString, startErrorString,          setStartErrorString)
    ATTRIBUTE (IString, endErrorString,            setEndErrorString)
    ATTRIBUTE (IString, numberOfNewtonStepsString, setNumberOfNewtonStepsString)
    ATTRIBUTE (IString, calculationTimeString,     setCalculationTimeString)      
    ATTRIBUTE (IString, relativeDifferenzString,   setRelativeDifferenzString)
				ATTRIBUTE (IString, gradientString,            setGradientString)
    
				/*! Strings fuer TeX-plot */
    ATTRIBUTE (IString, stepTeXString,                setStepTeXString)
    ATTRIBUTE (IString, startErrorTeXString,          setStartErrorTeXString)
    ATTRIBUTE (IString, endErrorTeXString,            setEndErrorTeXString)
    ATTRIBUTE (IString, numberOfNewtonStepsTeXString, setNumberOfNewtonStepsTeXString)
    ATTRIBUTE (IString, calculationTimeTeXString,     setCalculationTimeTeXString)      
    ATTRIBUTE (IString, relativeDifferenzTeXString,   setRelativeDifferenzTeXString)
    ATTRIBUTE (IString, gradientTeXString,            setGradientTeXString)
				
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
  
     DKTSDIterationInfoBlock ();
     DKTSDIterationInfoBlock (const DKTSDIterationInfoBlock& data);
     
     virtual ~DKTSDIterationInfoBlock();

  //@}

  /*!
  ***************************************************************************************
  * \name                             �ueberladene Operatoren
  ***************************************************************************************/

  //@{

  public:
     
     DKTSDIterationInfoBlock& operator =  (const DKTSDIterationInfoBlock& data);
					friend ostream&          operator << (ostream& s, DKTSDIterationInfoBlock& data);


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
     
     //! Setzt die default Werte
     DKTSDIterationInfoBlock setDefaultParameter();
     
     //! Hilfsmethode fuer plot
     bool plotEntry       (ostream& log,                  const ProtocolProperties& pP)                                        const;
     
     //! Hilfsmethode zum zentrieren der einzelnen Eintraege in Abhaengigkeit der Laenge der Zahlen und des Tabellenkopfes  
     bool centeringEntrys (const LongInt& lengthOfNumber, const LongInt& lengthString, LongInt& disString, LongInt& disNumber) const;
     
     //! Hilfsmethode fuer generatedPlotFiles
     bool plot2File       (const IString& blockIndex,     const ProtocolProperties& pP)                                        const;
     
     //! Erstellt Dateinamen blockIndex-name-dateiEndung
     IString fileName     (const IString& name,           const IString& blockIndex,   const IString& dateiEndung)             const;
     
     //! Erstellt Dateinamen datum-zeit-blockIndex-name-dateiEndung
     IString timeFileName (const IString& name,           const IString& blockIndex,   const IString& dateiEndung)             const;
  
    //! 

  //@}
 };


typedef  DKTSDIterationInfoBlock* DKTSDIterationInfoBlockPointer;

#endif // not defined __DKTSDIterationInfoBlock__
