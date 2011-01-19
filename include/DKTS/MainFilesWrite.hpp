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
 * \class    MainFilesWrite
 *
 * \brief    MainFilesWrite 
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

#ifndef __MainFilesWrite__
#define __MainFilesWrite__

#include "Macros.h"
#include <iostream>

#include <fstream>
#include "IString.hpp"
#include "ProtocolProperties.hpp"

using namespace std;

class  MainFilesWrite
 {
    DECLARE (MainFilesWrite)
				
				/*! zum setzen des Ausgabe Terminals                              */
				ATTRIBUTE (IString, terminal, setTerminal)
				
				/*! zum setzen des OutputDateiNamens                              */
				ATTRIBUTE (IString, output,   setOutput)
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
  
     MainFilesWrite ();
     
     virtual ~MainFilesWrite();

  //@}

  /*!
  ***************************************************************************************
  * \name                             �ueberladene Operatoren
  ***************************************************************************************/

  //@{

  public:

  //@}

  /*!
  ***************************************************************************************
  * \name                        oeffentliche Klassendienste
  ***************************************************************************************/

  //@{

  public:
    
     //! Erstellt den Kopf der gnu-Plot Projektdatei
     bool mainGnuPlotFileHead (ofstream& log, const IString& terminal=IString("windows"), const IString& outputFileName=IString("")) const;
     
     //! Erstellt das Ende der gnu-Plotdatei
     bool mainGnuPlotFileEnd  (ofstream& log) const;
     
     //! Erzeugt einen  main TeXFile Kopf bis zum input
     bool mainTeXFileHead (ofstream& log, const ProtocolProperties& pP) const;
     
     //! Erzeugt das main TeXFile Ende
     bool mainTeXFileEnd  (ofstream& log)     const;
     
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


typedef  MainFilesWrite* MainFilesWritePointer;

#endif // not defined __MainFilesWrite__
