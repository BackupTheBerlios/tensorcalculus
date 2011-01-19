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
 * \class    ProtocolProperties
 *
 * \brief    ProtocolProperties 
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

#ifndef __ProtocolProperties__
#define __ProtocolProperties__

#include "Macros.h"
#include <iostream>

#include "IString.hpp"
#include "Timer.hpp"

#include "ProtocolFormat.hpp"

using namespace std;

class  ProtocolProperties
 {
    DECLARE (ProtocolProperties)      
    
    /*! dient zum setzen von datum und Zeit                 */
    ATTRIBUTE(IString, time, setTime)
    
    /*! dient zum setzen von Genauigkeit und Spaltenabstand */
    PATTRIBUTE(ProtocolFormat, format, setFormat)
    
    /*! Angaben fuer die ProtocolTitelseite                 */
    ATTRIBUTE(IString, institutNamePart1, setInstitutNamePart1)
    ATTRIBUTE(IString, institutNamePart2, setInstitutNamePart2)
    ATTRIBUTE(IString, topicString,       setTopicString)
    ATTRIBUTE(IString, authorString,      setAuthorString)
    ATTRIBUTE(IString, cityString,        setCityString)
				
				/*! dient zum setzen des Namen der main - TeX Datei     */
				ATTRIBUTE(IString, TeXFileName,       setTeXFileName) 
				
				/*! dient zum setzen des Projektnamens                  */
				ATTRIBUTE(IString, projektName,       setProjektName)
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
  
     ProtocolProperties (const LongInt& prec=4, const LongInt& tabSize=3, const IString& position="ht", const LongInt& tableSize=35);
     ProtocolProperties (const ProtocolFormat& pF);
     ProtocolProperties (const ProtocolProperties& pP);
     
     virtual ~ProtocolProperties();

  //@}

  /*!
  ***************************************************************************************
  * \name                             �ueberladene Operatoren
  ***************************************************************************************/

  //@{

  public:
     
     ProtocolProperties& operator =  (const ProtocolProperties& pP);
     friend ostream&     operator << (ostream& s, const ProtocolProperties& pP);

     
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
     
     ProtocolProperties setDefaultParameter();
  
    //! 

  //@}
 };


typedef  ProtocolProperties* ProtocolPropertiesPointer;

#endif // not defined __ProtocolProperties__
