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
 * \class    IntegerArray
 *
 * \brief    IntegerArray 
 *
 *
 * \author   Mike Espig
 *
 * \ingroup  
 *
 *
 ******************************************************************************
 *
 *
 ******************************************************************************
 */

#ifndef __IntegerArray__
#define __IntegerArray__

#include "Macros.h"
#include <iostream>

#include "Random.hpp"

using namespace std;

class  IntegerArray
 {
    DECLARE (IntegerArray)
    
    ATTRIBUTE (LongInt, order, setOrder)
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
  
     IntegerArray (const LongInt& d);
     IntegerArray (const IntegerArray& v);
     IntegerArray (const LongInt& d, const LongInt& nMax);
     
     virtual ~IntegerArray();

  //@}

  /*!
  ***************************************************************************************
  * \name                             �berladene Operatoren
  ***************************************************************************************/

  //@{

  public:
  
      LongInt&      operator () (const LongInt& i) const { return (attr_elements[i]);}
      IntegerArray& operator  = (const IntegerArray& a);      
      IntegerArray& operator *= (const LongInt& multiplier);
      bool          operator == (const IntegerArray& a) const;
      bool          operator != (const IntegerArray& a) const;
      
  //@}

  /*!
  ***************************************************************************************
  * \name                        �ffentliche Klassendienste
  ***************************************************************************************/

  //@{

  public:
  
     IntegerArray& resize       (const LongInt& d);
     IntegerArray& setAllEntrys (const LongInt& value=0);     
     IntegerArray& setRand      (const LongInt& vMax, const LongInt& vMin=0);
     
     IntegerArray& setPwProductOf (const LongInt& alpha, const IntegerArray& v1, const IntegerArray& v2); 
     
     LongInt       sumOfAllEntry2 (const LongInt& nu) const;

  //@}


  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{

     friend ostream& operator << (ostream&, const IntegerArray& A);
     friend LongInt l1           (const IntegerArray& v);
     friend LongInt maxEntryOf   (const IntegerArray& v);     

  private:
  
     IntegerArray& allocateDataSpace (const LongInt& d);
     IntegerArray& deleteDataSpace   ();
  
     LongIntPointer attr_elements;
  
    //! 

  //@}
 };


typedef  IntegerArray* IntegerArrayPointer;

#endif // not defined __IntegerArray__
