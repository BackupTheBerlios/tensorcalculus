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
 * \class    DGKTSDNewton
 *
 * \brief    DGKTSDNewton 
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

#ifndef __DGKTSDNewton__
#define __DGKTSDNewton__


#include "Macros.h"
#include <iostream>
#include <iomanip>

#include "DKTS.hpp"
//#include "Protocol.hpp"
#include "DGKTSDFullMethod.hpp"


using namespace std;

class  DGKTSDNewton : public DGKTSDFullMethod
 {
    DECLARE (DGKTSDNewton)
        
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
   
     DGKTSDNewton ();

     virtual ~DGKTSDNewton();
  

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

  //@}


  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{



  private:

    DGKTSDNewton&  solveDirection  (const DKTS& a,  const DKTS& xi);    
    DGKTSDNewton&  useKrylovMethod (const DKTS& a,  const DKTS& xi);

    bool     evaluateHf            (const DKTS& a, const DKTS& xi, 
                                    const DKTS& v, DKTS& w) const;                                  
    bool     evaluateA             (const DKTS& a, const DKTS& xi, 
                                    const DKTS& v, DKTS& w,
                                    const LongReal& alpha, const LongReal& beta) const;                                  
    bool     evaluateAinv          (const DKTS& v, DKTS& w) const;
    
    bool     cg                    (const DKTS& a, const DKTS& xi, const LongInt& maxSteps);        

    //!

  //@}
 };


typedef  DGKTSDNewton* DGKTSDNewtonPointer;

#endif // not defined __DGKTSDNewton__
