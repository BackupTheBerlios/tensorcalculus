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
 * \class    Protocol
 *
 * \brief    Protocol 
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

#ifndef __Protocol__
#define __Protocol__

#include "Macros.h"

#include <iostream>
#include <iomanip>
#include <typeinfo>

#include "SimpleException.hpp"

#include "ProtocolData.hpp"

#include "ProtocolProperties.hpp"
#include "ProtocolFormat.hpp"
#include "MainFilesWrite.hpp"


#include "DKTSDIterationInfo.hpp"
#include "SimpleIterationData.hpp"
#include "DKTSDData.hpp"

#include "DKTSDIterationInfoBlock.hpp"
#include "DescriptionData.hpp"
#include "SimpleIterationDataBlock.hpp"
#include "DKTSDDataBlock.hpp"

using namespace std;

class  Protocol
 {
    DECLARE (Protocol)      
    
    CONTAINS  (ProtocolDataPointer, Entry)
    
    PATTRIBUTE (ProtocolProperties,  properties, setProperties)
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public: 
  
    Protocol ();
    Protocol (const ProtocolProperties& pP);
    Protocol (const Protocol& p);
     
    virtual ~Protocol ();

  //@}

  /*!
  ***************************************************************************************
  * \name                             �ueberladene Operatoren
  ***************************************************************************************/

  //@{

  public:

     Protocol&       operator =  (const Protocol& p);
     friend ostream& operator << (ostream& s, const Protocol& p);

  //@}

  /*!
  ***************************************************************************************
  * \name                        oeffentliche Klassendienste
  ***************************************************************************************/

  //@{

  public:

     /*! plot2 gibt das Protokol auf dem Bildschirm aus                         */
     bool plot2         (ostream& log)                                              const;
     
     /*! plot2PlotFile schreibt die Daten des Protokols in .dat Datein          */
     bool plot2PlotFile (const IString& mainFoldername)                             const;
     
     /*! plot2TeXFiles erzeugt TeX Files aus den einzelnen Datenbloecken        */
     bool plot2TeXFile  (const IString& mainFoldername, const IString& ProjektName) const;
     
     /*! plotAll ist eine Zusammenfassung von plot2PlotFile und plot2TeXFile    */
     bool plotAll       (const IString& mainFoldername, const IString& name)        const;

     
     /*! add haengt ein Protocoldata vom typ ... ans Protocol an                */
     Protocol& add (const SimpleIterationData& data);
     Protocol& add (const DKTSDData& data);
     Protocol& add (const DKTSDIterationInfo& data);
     
     Protocol& add (const DescriptionData& data);
     Protocol& add (const SimpleIterationDataBlock& data);
     Protocol& add (const DKTSDDataBlock& data);
     Protocol& add (const DKTSDIterationInfoBlock& data);

  //@}


  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{


  private:
  
     Protocol  setDefaultParameter ();
     Protocol& include             (const ProtocolDataPointer element);     

     IString   createFolderTree    (const IString& mainFoldername, const ProtocolProperties& pP, const IString& name) const;
     IString   setZeros            (const LongInt& number)                                                            const;
     
     bool      plotAll2File        (const IString& mainFoldername, const IString& name)      		                       const;
     
   //! 
   
  //@}
 };


typedef  Protocol* ProtocolPointer;

#endif // not defined __Protocol__
