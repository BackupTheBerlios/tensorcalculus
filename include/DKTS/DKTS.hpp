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
 * \class    DKTS
 *
 * \brief    DKTS 
 *
 *
 * \author   Mike Espig
 *
 * \ingroup  KTS
 *
 *
 ******************************************************************************
 *
 *
 ******************************************************************************
 */

#ifndef __DKTS__
#define __DKTS__

#include "Macros.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>

using namespace std;

#include "KTS.hpp"
#include "DKTVector.hpp"
#include "IntegerArray.hpp"

#include "RMatrix.hpp"
#include "RVector.hpp"
#include "DescriptionData.hpp"


#include <cstdlib>    
#include <new> 
    

class  DKTS : public KTS
 {    
    DECLARE (DKTS)

    static char     const_u;
    static char     const_jobz;
    static char     const_notConjTrans;
    static char     const_trans;
    static LongInt  const_inc;
    static LongReal const_eins;
    static LongReal const_minus_eins;
    static LongReal const_null;        
    
    
    ATTRIBUTE (LongInt, dimension, setDimension)
    ATTRIBUTE (LongInt, kn,        setKN)
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
   
     DKTS(const LongInt& d=1, const LongInt& k=1, const LongInt& n=1);
     DKTS(const DKTS& A);     
    
     virtual ~DKTS();
  

  //@}

  /*!
  ***************************************************************************************
  * \name                             �berladene Operatoren
  ***************************************************************************************/

  //@{

  public:

     inline LongReal& operator () (const LongInt& j, const LongInt& mu, const LongInt& i) const
                        {
                           const LongInt index = i + attr_n*j + attr_kn*mu;
                          return (attr_values[index]);
                        }

     inline DKTVector& operator () () const 
                        {
                          return (attr_vecOfTensor[0]);
                        }
                        
     inline DKTVector& operator () (const LongInt& mu) const
                        {                            
                          return (attr_vecOfDimension[mu]);
                        }

     inline DKTVector& operator () (const LongInt& j, const LongInt& mu) const
                        {                            
                           const LongInt index = j + attr_k*mu;                           
                          return (attr_vecOfRepresentationVector[index]);
                        }
      
     LongReal operator () (const IntegerArray& l) const;
     LongReal operator () (const IntegerArray& l, const LongInt& jMax, const LongInt& jMin=0) const;    
      
     DKTS& operator  = (const DKTS& A);     
     DKTS& operator *= (const LongReal& alpha);                             

  //@}

  /*!
  ***************************************************************************************
  * \name                        �ffentliche Klassendienste
  ***************************************************************************************/

  //@{

  public:     

     DKTS&    resize          (const LongInt& d, const LongInt& k, const LongInt& n);

     DKTS&    setSumOf        (const LongReal& alpha, const DKTS& A1, const LongReal& beta, const DKTS& A2);
     DKTS&    setSumOf        (const LongReal& a, const DKTS& A1, const LongReal& b, const DKTS& A2, const LongReal& c, const DKTS& A3);
     
     DKTS&    setProductOf    (const LongReal& alpha, const DKTS& A1, const DKTS& A2, const LongReal& beta);
     bool     evaluateAt      (const LongReal& alpha, const DKTS& v, const LongReal& beta, DKTS& w) const;
     
     bool     evaluateOrthProjAt     (const DKTS& v, DKTS& w) const;
     bool     evaluateOrthCompProjAt (const DKTS& v, DKTS& w) const;

     DKTS&    reScaled         ();
     DKTS&    setRand          (const LongReal eps=1.0e-3);
     DKTS&    setNull          ();
     //edit 2010-01-23
     DKTS&    setAllEntries    (const LongReal& a);
     //end edit
     DKTS&    scale            (const LongReal eps=0.5, const LongReal& alpha=1.0);
     DKTS&    balanced         ();
     DKTS&    normalizedAddens ();
     
     DKTS&    setCrossApproximationOf (const DKTS& A, IntegerArray& i, const LongInt& maxSteps=3);
     DKTS&    setCrossApproximationOf (const DKTS& A, const LongInt& maxSteps=3);
     DKTS&    setAlsApproximationOf   (const DKTS& A, const LongInt& maxSteps=5, const bool& useCross=true);
     DKTS&    successiveApproximation (const DKTS& A, const LongReal& eps, const LongInt& maxSteps=8);
     
     DKTS&    improveApproximation    (const DKTS& a, const LongInt& maxStepsAdd=10, const LongInt& maxStepsAls=2);
     
     DKTS&    setCopyOf       (const DKTS& A);    
     DKTS&    isPartialCopyOf (const DKTS& A);
     DKTS&    copyFromTo      (const DKTS& A, const LongInt& os, const LongInt& j0, const LongInt& j1);
     DKTS&    setSummandOf    (const DKTS& A, const LongInt& index);     
     LongReal setMainPartOf   (const DKTS& A, const LongReal& eps);     
     LongInt  addR1Addent     (const DKTS& x);
         
     DKTS&    setHadamardProductOf (const LongReal& a, const DKTS& A, const DKTS& B);
          
     DKTS&    readDataFrom    (const IString& file);
     DKTS&    readDataFrom    (const IString& file, const LongInt& k1);
     bool     writeDataTo     (const IString& file) const;     
     
     bool     writeDiagonalDataTo (const IString& file) const;
          
     bool     writeIndexSet   (ostream& s) const;
     bool     writeIndexSet   (DescriptionData& s) const;
     bool     writeAllNorms   (ostream& s) const;
     bool     writeParameter  (ostream& s) const;

     LongReal distanceTo              (const DKTS& A) const;
     LongReal frobeniusNormOfSummand  (const LongInt& i) const;
     LongReal distanceTo              (const LongInt& j1, const LongInt& j2) const;

     DKTS&    regeneratedBy           (const DKTS& x, const DKTS& Z);
     DKTS&    setOrthogonal2          (const DKTS& A, const LongReal& error=1.0e-6);
     DKTS&    setCoefficientsSystemOf (const DKTS& A, const DKTS& Z);
     DKTS&    sortIndexSet            ();
     DKTS&    setModellExample        (const LongInt& d, const LongInt& n);

     bool     analyseVectorSystem     () const;
     bool     analyseGSystem          () const;

     friend ostream& operator <<      (ostream& s, const DKTS& A);     
     friend LongReal frobeniusNorm    (const DKTS& A);
     friend LongReal innerProduct     (const DKTS& A, const DKTS& B);     

     friend LongReal maximumValueOf   (const DKTS& a, IntegerArray& indexOfMaximum);
     friend LongReal minimumValueOf   (const DKTS& a, IntegerArray& indexOfMinimum);
     
     friend LongReal upperBoundOfMaximumValue (const DKTS& a);
     friend LongReal lowerBoundOfMinimumValue (const DKTS& a);
     

  //@}

  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{

  private:

     LongInt  partition    (LongInt low, LongInt high, LongRealPointer f, LongIntPointer indexSet);
     DKTS&    quickSort    (LongInt low, LongInt high, LongRealPointer f, LongIntPointer indexSet);
     DKTS&    swapIndex    (const LongInt& i, const LongInt& j, LongRealPointer f, LongIntPointer indexSet);
     
     
     friend LongReal computeMaximumValueOf (const DKTS& a, IntegerArray& indexOfMaximum);
       
     DKTS& allocateDataSpace    (const LongInt& d, const LongInt& k, const LongInt& n);
     DKTS& deleteDataSpace      ();     
     
     LongRealPointer  attr_values;
     DKTVectorPointer attr_vecOfTensor;
     DKTVectorPointer attr_vecOfDimension;
     DKTVectorPointer attr_vecOfRepresentationVector;

    //! 

  //@}
 };

typedef  DKTS* DKTSPointer;

#endif // not defined __DKTS__
