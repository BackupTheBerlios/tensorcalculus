/*
 * Copyright (C) Mike Espig
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
 * \class    DescriptionData
 *
 * \brief    DescriptionData 
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

#ifndef __DescriptionData__
#define __DescriptionData__

#include "Macros.h"
#include <iostream>

#include "BlockData.hpp"

using namespace std;

class  DescriptionData : public BlockData
 {
    DECLARE (DescriptionData)      
    
				/*! erzeugt einen IString-Block */
    CONTAINS (IString, String)
    

  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
  
     DescriptionData();
     
     virtual ~DescriptionData();

  //@}

  /*!
  ***************************************************************************************
  * \name                             �ueberladene Operatoren
  ***************************************************************************************/

  //@{

  public:
   
			  friend ostream& operator << (ostream& s, const DescriptionData& data);
					
  //@}

  /*!
  ***************************************************************************************
  * \name                        oeffentliche Klassendienste
  ***************************************************************************************/

  //@{

  public:
    
     /*! Ausgabe der Strings auf dem uebergebenen ostream             */
     virtual bool    plot               (ostream& log,              const ProtocolProperties& pP) const;
     
     /*! Zusammenfassung von generatedPlotFiles und generatedTeXFiles */
     virtual bool    generatedAll       (const IString& blockIndex, const ProtocolProperties& pP) const;

     /*! rein virtuelle Methode                                       */
     virtual bool    generatedPlotFiles (const IString& blockIndex, const ProtocolProperties& pP) const {return false;}
     
     /*! Ausgabe der Strings in eine .tex Datei                       */
     virtual IString generatedTeXFiles  (const IString& blockIndex, const ProtocolProperties& pP) const;
     
     /*! fuegt vor und nach dem String s das Symbol $ ein             */
             IString addMathString      (const IString& s);

  //@}


  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{


  private:
     
     /*! Erstellt Dateinamen datum-zeit-blockIndex-name-dateiEndung   */
     IString timeFileName (const IString& name, const IString& blockIndex, const IString& dateiEndung) const;
     
     /*! Erstellt Dateinamen blockIndex-name-dateiEndung              */
     IString fileName     (const IString& name, const IString& blockIndex, const IString& dateiEndung) const;

  
    //! 

  //@}
 };


typedef  DescriptionData* DescriptionDataPointer;

#endif // not defined __DescriptionData__
