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
 * \class    DGKTSD
 *
 * \brief    DGKTSD 
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

#ifndef __DGKTSD__
#define __DGKTSD__


#include "Macros.h"
#include <iostream>
#include <iomanip>

#include "DKTS.hpp"
#include "Protocol.hpp"
#include "DGKTSDFullMethod.hpp"


using namespace std;

class  DGKTSD : public DGKTSDFullMethod
 {
    DECLARE (DGKTSD)
    
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
   
     DGKTSD ();

     virtual ~DGKTSD();
  

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

    DGKTSD&  solveDirection       (const DKTS& a,  const DKTS& xi);

    
    DGKTSD&  useKrylovMethod      (const DKTS& a,  const DKTS& xi);


    DGKTSD&  setParameter         (const LongReal& norm, const LongReal& error, const LongInt& r);
    DGKTSD&  setParameter         (const LongReal& norm, const LongInt& r);

    bool     evaluateWHf         (const DKTS& a, const DKTS& xi, 
                                  const DKTS& v, DKTS& w) const;
                                  
    bool     evaluateHf          (const DKTS& a, const DKTS& xi, 
                                  const DKTS& v, DKTS& w) const;
                                  
    bool     evaluateA           (const DKTS& a, const DKTS& xi, 
                                  const DKTS& v, DKTS& w,
                                  const LongReal& alpha, const LongReal& beta) const;
                                  
    bool     evaluateAinv        (const DKTS& v, DKTS& w) const;
    
    bool     cr                  (const DKTS& a, const DKTS& xi);
    bool     crA                 (const DKTS& a, const DKTS& xi);



    //!

  //@}
 };


typedef  DGKTSD* DGKTSDPointer;

#endif // not defined __DGKTSD__
