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
 * \class    Timer
 *
 * \brief    Timer
 *
 *
 * \author   Mike Espig
 *
 * \ingroup  Base
 *
 *
 ******************************************************************************
 *
 *
 ******************************************************************************
 */

#ifndef __Timer__
#define __Timer__

#include "Macros.h"
#include "IString.hpp"

using namespace std;

#include <time.h>
//#include <windows.h>


class  Timer 
 {
    DECLARE (Timer)
   
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
    //! Konstruktor
    /*! Dies ist der Standardkonstruktor*/
    Timer(){};

    //! Destruktor
    virtual ~Timer(){};

  //@}

  /*!
  ***************************************************************************************
  * \name                             �berladene Operatoren
  ***************************************************************************************/

  //@{

  public:

  //@}

  /*!
  ***************************************************************************************
  * \name                        �ffentliche Klassendienste
  ***************************************************************************************/

  //@{

  public:

    // Das Zeitmessen beginnen
    Timer& startTiming();

    // die seit dem letzten aufruf von 'startTiming' verstrichene Zeit in
    // Sekunden (Millisekunden genau)
    LongReal elapsedTimeSec();

    //! Erstellt Datum
    IString date() const;
 
    //! Erstellt Zeitstempel in der Form: datum-Zeit
    IString timeStamp() const;


  //@}


  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{


  private:

    LongReal startTimeS;
    IString  setZeros(const LongInt& number)const;


    //! 

  //@}
 };

#endif // not defined __Timer__
