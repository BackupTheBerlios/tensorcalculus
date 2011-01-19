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
 * \class    DGKTSDFullMethod
 *
 * \brief    DGKTSDFullMethod 
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

#ifndef __DGKTSDFullMethod__
#define __DGKTSDFullMethod__


#include "Macros.h"
#include <iostream>
#include <iomanip>


#include "DGKTSDecomposer.hpp"
#include "DKTS.hpp"


using namespace std;

class  DGKTSDFullMethod : public DGKTSDecomposer
 {
    DECLARE (DGKTSDFullMethod)
     
    friend class DGKTSDNewton;  
    friend class DGKTSD;
        
    static char     const_u;
    static char     const_jobz;
    static char     const_notConjTrans;
    static char     const_trans;
    static LongInt  const_inc;
    static LongReal const_eins;
    static LongReal const_minus_eins;
    static LongReal const_null;        
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
   
     DGKTSDFullMethod ();

     virtual ~DGKTSDFullMethod();
  

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
  
    virtual DKTSDIterationInfo startIteration (const DKTS& a, DKTS& xi, const LongReal& normA, DKTSDDataBlock& dataBlock);
  
  
    LongReal optimalOmega ();

           
    virtual DGKTSDFullMethod&  solveDirection (const DKTS& a,  const DKTS& xi)=0;        

    DGKTSDFullMethod& initialAllValues   (const DKTS& a, DKTS& xi);
    DGKTSDFullMethod& initialAllValues   (const DKTS& a, DKTS& xi, const LongReal& alpha); 
    
    DGKTSDFullMethod& initialValuesOfax  (const DKTS& a, const DKTS& xi);
    DGKTSDFullMethod& initialValuesOfax  (const LongReal& alpha);
    DGKTSDFullMethod& prepareValuesOfadx (const DKTS& a, const DKTS& xi);

    DGKTSDFullMethod& initialValuesOfx   (const DKTS& xi);
    DGKTSDFullMethod& initialValuesOfa   (const DKTS& a, const DKTS& xi);

    DGKTSDFullMethod& invertA            (const LongReal& alpha, const LongReal& beta);

    DGKTSDFullMethod& solveGradient      (const DKTS& a,  const DKTS& xi);    
    bool              evaluateW          (const DKTS& v, DKTS& w) const;
    bool              evaluateW          (DKTS& w) const;

    LongReal F                            () const;
    LongReal f                            () const;
    LongReal g2                           () const;
    LongReal F                            (const LongReal& alpha) const;

    //!

  //@}
 };


typedef  DGKTSDFullMethod* DGKTSDFullMethodPointer;

#endif // not defined __DGKTSDFullMethod__
