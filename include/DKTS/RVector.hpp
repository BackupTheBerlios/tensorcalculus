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
 * \class    RVector
 *
 * \brief    Allgemeiner, n-dimensionaler LongRealer Vektor
 *
 *
 * \author   Mike Espig
 *
 *
 *
 * \ingroup  MathBase 
 *
 ******************************************************************************
 *
 *
 */


#ifndef __RVector__
#define __RVector__

#include "Macros.h"

#include <iostream>
#include <time.h>

#include "AVector.hpp"
#include "DKTVector.hpp"


using namespace std;


class RMatrix;
class RkRMatrix;


class  RVector : public AVector
 {

    DECLARE (RVector)
    CRITICAL_ATTRIBUTE(int, refCounter, setRefCounter)

  /*
  ***************************************************************************************
  * ! @name                       Konstruktion / Destruktion
  ***************************************************************************************/
  //@{

  public:
     /** Dies ist der Standartkonstruktor **/
     RVector (LongInt dim=0); 
     RVector (const RVector& v);
     RVector (const DKTVector& v);
     RVector(const LongInt& dim, const LongReal* add);
    ~RVector ();

  //@}


  /*
  ***************************************************************************************
  * ! \name                       �berladene Operatoren
  ***************************************************************************************/
  //@{

  public:
    RVector&  operator  = (const RVector& v);
    RVector&  operator += (const RVector& v) ;
    RVector&  operator -= (const RVector& v) ;
    RVector&  operator *= (const LongReal& a) ;
    RVector&  operator /= (const LongReal& a) ;
    RVector	  operator  + (const RVector& v) const;
    RVector	  operator  - (const RVector& v) const;
    LongReal&  operator () (LongInt i) const { return (_pkomp[i]); }
    bool      operator == (const RVector& v) const;
    bool      operator != (const RVector& v) const;
    

  //@}


  /*
  ***************************************************************************************
  * ! \name                       �ffentliche Klassendienste
  ***************************************************************************************/
  //@{

  public:
   
    RVector&  newDimension (LongInt newdim);
    RVector&  resize       (LongInt newdim);
    RVector&  update	      (const LongReal& a, const RVector& d);
    RVector&  setNull      ();
    RVector&  setRand      (const LongReal& eps=2.5e-6);
    
    LongReal  smallestEntry() const;


  //@}



  /*
  ***************************************************************************************
  * ! \name                       Implementierungsbereich
  ***************************************************************************************/
  //@{

	friend ostream&	 operator <<   (ostream & s,const RVector &v);
	friend istream&	 operator >>   (istream & s, RVector &v);
	friend RVector	 operator *    (LongReal alpha, const RVector &v);
	friend RVector	 operator /    (const RVector &v, LongReal alpha);
	friend LongReal	 innerProduct  (const RVector &v1,const RVector &v2); 		//das kanonische Skalarprodukt
	friend LongReal	 L1			   (const RVector &v);							//L1 Norm
	friend LongReal	 L2			   (const RVector &v);							//euklidische Norm
	friend LongReal	 Max		   (const RVector &v);							//die Maximumnorm
	friend LongReal	 Min		   (const RVector &v);

  protected:

  private:

	LongReal	*_pkomp;	//Zeiger auf die RVectorkomponenten

  //@}
};


typedef  RVector* RVectorPointer;

#endif /* RVector_h */
