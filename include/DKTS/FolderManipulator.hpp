/*
 * Copyright (C) Marcel Schuster
 *               2011 Philipp Wähnert
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
 * \class    FolderManipulator
 *
 * \brief    FolderManipulator 
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

#ifndef __FolderManipulator__
#define __FolderManipulator__

#include "Macros.h"
#include <iostream>

class  IString;

class  FolderManipulator
 {
    DECLARE (FolderManipulator)      
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:

     FolderManipulator () { }

     virtual ~FolderManipulator() { }

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
     
     //! Anlegen eines Ordners
     bool createFolder (const IString& name) const;
     
     //! Wechseln in den Ordner
     bool changeFolder (const IString& name) const;
     
     //! oeffnen des Ordners
     bool openFolder   (const IString& name) const;
     
     //! schliessen des Ordners
     bool closeFolder  (const IString& name) const;
     
     //! wechselt einen Ordner hoeher
     bool folderUp     () const;
     
     //! gibt das aktuelle Arbeitsverzeichnis zurueck
     IString currentFolder()                    const;
 
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


typedef  FolderManipulator* FolderManipulatorPointer;

#endif // not defined __FolderManipulator__
