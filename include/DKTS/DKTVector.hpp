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
 * \class    DKTVector
 *
 * \brief    DKTVector 
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

#ifndef __DKTVector__
#define __DKTVector__

#include "Macros.h"
#include "KTVector.hpp"
#include "SimpleException.hpp"

#include "Random.hpp"

#include <iostream>
#include <time.h>

using namespace std;

class  DKTVector : public KTVector
 {
    DECLARE (DKTVector)
    
    static LongInt inc;
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
   
     DKTVector  ();
     DKTVector  (LongIntPointer n, LongRealPointer element);
     virtual ~DKTVector ();



  //@}

  /*!
  ***************************************************************************************
  * \name                             �berladene Operatoren
  ***************************************************************************************/

  //@{

  public:
  
     inline LongReal&  operator () () const
                        {
                          return attr_element[0];
                        }
                 
     inline LongReal&  operator () (const LongInt& i) const
                        {                           
                          return attr_element[i];
                        }

     DKTVector& operator  = (const DKTVector& x);
     DKTVector& operator += (const DKTVector& x);
     DKTVector& operator -= (const DKTVector& x);     
     DKTVector& operator *= (const LongReal&  alpha);
     DKTVector& operator /= (const LongReal&  alpha);

     friend ostream& operator << (ostream& os, const DKTVector& x);

     DKTVector& setCopyOf (const LongReal& alpha, const DKTVector& x);

  //@}

  /*!
  ***************************************************************************************
  * \name                        �ffentliche Klassendienste
  ***************************************************************************************/

  //@{

  public:
  
     DKTVector& set          (LongInt& n, LongReal& element);
  
     DKTVector& add          (const LongReal& alpha, const DKTVector& x, const LongReal& beta);
     DKTVector& update       (const LongReal& alpha, const DKTVector& x);
     DKTVector& pwProduct    (const LongReal& alpha, const DKTVector& x, const DKTVector& y);
     
     DKTVector& setNull      ();
     DKTVector& setAllEntrys (const LongReal& v);
     DKTVector& setRand      (const LongReal& eps=1.0e2);
               
     LongReal   normalized   ();     
     LongInt    indexOfMax   () const;

     friend LongReal innerProduct  (const DKTVector& x, const DKTVector& y);
     friend LongReal l2            (const DKTVector& x);
     
     friend LongReal maximumNormOf  (const DKTVector& x);
     
     friend LongReal maximumValueOf (const DKTVector& x, LongInt& indexOfMaxValue);
     friend LongReal minimumValueOf (const DKTVector& x, LongInt& indexOfMinValue);
  
  
  //@}
  
  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{

  private:
    
     LongRealPointer attr_element;
  
    //! 

  //@}
 };

typedef  DKTVector* DKTVectorPointer;


#endif // not defined __DKTVector__
