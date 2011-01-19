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
 * \class    CVector
 *
 * \brief    Allgemeiner, n-dimensionaler complexer Vektor
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


#ifndef __CVector__
#define __CVector__

#include "Macros.h"
#include <iostream>
#include "Complex.hpp"
#include "AVector.hpp"
#include "RVector.hpp"

using namespace std;



class CMatrix;
class RkCMatrix;

class  CVector : public AVector
 {

    DECLARE (CVector)
    CRITICAL_ATTRIBUTE(int, refCounter, setRefCounter)


  /*
  ***************************************************************************************
  * ! @name                       Konstruktion / Destruktion
  ***************************************************************************************/
  //@{

  public:
     /** Dies ist der Standartkonstruktor **/
     CVector (LongInt dim=0); 
     CVector (const CVector& v);
     CVector (const RVector& v);
     CVector(const LongInt& dim, const std::complex<double>* add);
    ~CVector ();

  //@}


  /*
  ***************************************************************************************
  * ! \name                       �berladene Operatoren
  ***************************************************************************************/
  //@{

  public:
    CVector&  operator  = (const CVector& v);
    CVector&  operator += (const CVector& v) ;
    CVector&  operator -= (const CVector& v) ;
    CVector&  operator *= (const std::complex<double>& a) ;
    CVector&  operator /= (const std::complex<double>& a) ;
    CVector	  operator  + (const CVector& v) const;
    CVector	  operator  - (const CVector& v) const;
    std::complex<double>&  operator () (LongInt i) const { return (_pkomp[i]); }
    bool      operator == (const CVector& v) const;
    bool      operator != (const CVector& v) const;
    

  //@}


  /*
  ***************************************************************************************
  * ! \name                       �ffentliche Klassendienste
  ***************************************************************************************/
  //@{

  public:
   
    CVector&  newDimension (LongInt newdim);
    RVector   realPart     () const;
    RVector   imagPart     () const;
    CVector&  update	   (const std::complex<double>& a, const CVector& d);
    CVector&  setNull      ();

  //@}



  /*
  ***************************************************************************************
  * ! \name                       Implementierungsbereich
  ***************************************************************************************/
  //@{

	friend ostream&	 operator <<   (ostream & s,const CVector &v);
	friend istream&	 operator >>   (istream & s, CVector &v);
	friend CVector	 operator *    (std::complex<double> alpha, const CVector &v);
	friend CVector	 operator /    (const CVector &v, std::complex<double> alpha);
	friend std::complex<double>	 innerProduct  (const CVector &v1,const CVector &v2); 		//das kanonische Skalarprodukt
	friend LongReal	 L1			   (const CVector &v);							//L1 Norm
	friend LongReal	 L2			   (const CVector &v);							//euklidische Norm
	friend LongReal	 Max		   (const CVector &v);							//die Maximumnorm

  protected:

  private:

	std::complex<double>	*_pkomp;	//Zeiger auf die CVectorkomponenten

  //@}
};


#endif /* CVector_h */
