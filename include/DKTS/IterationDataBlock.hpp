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
 * \class    IterationDataBlock
 *
 * \brief    IterationDataBlock 
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

#ifndef __IterationDataBlock__
#define __IterationDataBlock__

#include "Macros.h"
#include <iostream>

#include "BlockData.hpp"

using namespace std;

class  IterationDataBlock : public BlockData
 {
    DECLARE (IterationDataBlock)      
   
    /*! dient zum setzen der Tabellenbeschrieftung im LaTeX */
				ATTRIBUTE(IString, captionString, setCaptionString)
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
     
     IterationDataBlock ();
     IterationDataBlock (const IterationDataBlock& dB);

     virtual ~IterationDataBlock();
     

  //@}

  /*!
  ***************************************************************************************
  * \name                             �ueberladene Operatoren
  ***************************************************************************************/

  //@{

  public:
     
     IterationDataBlock& operator = (const IterationDataBlock& dB);
     
  //@}

  /*!
  ***************************************************************************************
  * \name                        oeffentliche Klassendienste
  ***************************************************************************************/

  //@{

  public:

     /*! nur virtuelle Methoden */
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
     
     IterationDataBlock setDefaultParameter();

    //! 

  //@}
 };


typedef  IterationDataBlock* IterationDataBlockPointer;

#endif // not defined __IterationDataBlock__
