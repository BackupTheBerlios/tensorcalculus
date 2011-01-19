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
 * \class    ProtocolData
 *
 * \brief    ProtocolData 
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

#ifndef __ProtocolData__
#define __ProtocolData__

#include "Macros.h"
#include <iostream>

#include <iomanip>			//! setPrecision, scientific, uppercase, setw
#include <fstream>			//! ofstream

#include "IString.hpp"
#include "Timer.hpp"

#include "ProtocolFormat.hpp"
#include "ProtocolProperties.hpp"
//#include "FolderManipulator.hpp"
#include "MainFilesWrite.hpp"

using namespace std;

class  ProtocolData
 {
    DECLARE (ProtocolData) 
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:

     ProtocolData();
     
     virtual ~ProtocolData();


  //@}

  /*!
  ***************************************************************************************
  * \name                             �ueberladene Operatoren
  ***************************************************************************************/

  //@{

  public:
     
     virtual bool operator == (const ProtocolData& D) const {return false;}
             bool operator != (const ProtocolData& D) const;
     
  //@}

  /*!
  ***************************************************************************************
  * \name                        oeffentliche Klassendienste
  ***************************************************************************************/

  //@{

  public:
     
					/*! rein virtuelle Methoden */
     virtual bool    plot               (ostream& log,              const ProtocolProperties& pP) const {return false;}
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
  
  
    //! 

  //@}
 };


typedef  ProtocolData* ProtocolDataPointer;

#endif // not defined __ProtocolData__
