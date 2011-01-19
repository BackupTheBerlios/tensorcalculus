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
 * \class    ProtocolFormat
 *
 * \brief    ProtocolFormat 
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

#ifndef __ProtocolFormat__
#define __ProtocolFormat__

#include "Macros.h"
#include <iostream>

#include "IString.hpp"

using namespace std;

class  ProtocolFormat
 {
    DECLARE (ProtocolFormat)   
    
    /*! dient zum setzen der Genauigkeit                         */
    ATTRIBUTE (LongInt, precision, setPrecision)
    
    /*! dient zum setzen des Spaltenabstandes                    */
    ATTRIBUTE (LongInt, tabSize,   setTabSize)
    
    /*! plazieren der Tabellen in LaTeX                          */
    ATTRIBUTE (IString, tablePosition, setTablePosition)
				
			 /*! Gibt die max. Anzahl der Zeilen in einer TeX Tabelle an  */
    ATTRIBUTE(LongInt, maxTableSize, setMaxTableSize)         
     
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
  
     ProtocolFormat();
					ProtocolFormat(const LongInt& prec, const LongInt& tabSize, const IString& position, const LongInt& tableSize);
     ProtocolFormat(const ProtocolFormat& pF);
     
     virtual ~ProtocolFormat();

  //@}

  /*!
  ***************************************************************************************
  * \name                             �ueberladene Operatoren
  ***************************************************************************************/

  //@{

  public:

     ProtocolFormat& operator =  (const ProtocolFormat& pF);
     friend ostream& operator << (ostream& s, const ProtocolFormat& pF);


  //@}

  /*!
  ***************************************************************************************
  * \name                        oeffentliche Klassendienste
  ***************************************************************************************/

  //@{

  public:

  //@}


  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{


  private:
  
     ProtocolFormat setDefaultParameter();
  
    //! 

  //@}
 };


typedef  ProtocolFormat* ProtocolFormatPointer;

#endif // not defined __ProtocolFormat__
