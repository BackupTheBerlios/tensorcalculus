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
 * \defgroup KTSD
 * \brief Diese Bibliothek enth�lt grundlegende Dienstklassen zum k�rzen von Summen von Kronecker-Tensor-Produkten
 *   
 ******************************************************************************
*/

/*!
 ******************************************************************************
 * \class    DGKTSDecomposer
 *
 * \brief    DGKTSDecomposer 
 *
 *
 * \author   Mike Espig
 *
 * \ingroup  KTSD
 *
 *
 ******************************************************************************
 *
 *
 ******************************************************************************
 */

#ifndef __DGKTSDecomposer__
#define __DGKTSDecomposer__


#include "Macros.h"
#include <iostream>
#include <iomanip>

#include "DKTS.hpp"
#include "Protocol.hpp"
#include "DGKTSDDataBase.hpp"
#include "DKTSDIterationInfo.hpp"


using namespace std;

class  DGKTSDecomposer : public DGKTSDDataBase
 {
    DECLARE (DGKTSDecomposer)
    
    friend class DGKTSDNewCompare;
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
   
     DGKTSDecomposer ();

     virtual ~DGKTSDecomposer();
  
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

     DKTSDIterationInfo decompose (DKTS& a, DKTS& xi, DKTSDDataBlock& dataBlock);
     DKTSDIterationInfo decompose (DKTS& a, DKTS& xi, const LongReal& normA, DKTSDDataBlock& dataBlock);

  //@}


  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{

  private:
 
     DGKTSDecomposer& setDefaultParameter (const LongReal& normA);

     virtual DKTSDIterationInfo startIteration (const DKTS& a, DKTS& xi, const LongReal& normA, DKTSDDataBlock& dataBlock) = 0;
     

  //@}
 };


typedef  DGKTSDecomposer* DGKTSDecomposerPointer;

#endif // not defined __DGKTSDecomposer__
