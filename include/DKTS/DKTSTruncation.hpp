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
 * \class    DKTSTruncation
 *
 * \brief    DKTSTruncation 
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

#ifndef __DKTSTruncation__
#define __DKTSTruncation__

#include "Macros.h"
#include <iostream>

#include "DGKTSD.hpp"
#include "DGKTSDNewton.hpp"
#include "Timer.hpp"
#include "Random.hpp"


#include "IntegerArray.hpp"

#include "Protocol.hpp"

using namespace std;

class  DKTSTruncation
 {
    DECLARE (DKTSTruncation)
    
    PATTRIBUTE(Protocol, truncationLog,   setTruncationLog)
    PATTRIBUTE(Protocol, truncationLogR1, setTruncationLogR1)

    ATTRIBUTE (LongReal, preCalculationTime, setPreCalculationTime)
    ATTRIBUTE (IString,  inputFileName,      setInputFileName)
    
    ATTRIBUTE (LongReal, normA, setNormA)
    
    ATTRIBUTE (LongInt, originalRank, setOriginalRank)
    
    friend class DKTSTest;
    friend class DKTSTruncation2Eps;
    friend class ProblemCases3D;
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
     
     DKTSTruncation (const LongInt& d=1, const LongInt& r=1, const LongInt& R=1, const LongInt& n=1);
     DKTSTruncation (const IString& fileName, const LongInt& r=1);
     DKTSTruncation (const IString& fileName, const IString& coFileName, const IString& ortFileName, const LongInt& r=1);
     DKTSTruncation (const DKTS& a, const DKTS& A, const DKTS& Z, const LongInt& r);
     DKTSTruncation (const DKTS& A, const LongInt& r);
     DKTSTruncation (const DKTS& A);
     
     
     virtual ~DKTSTruncation();

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
  
     DKTSDIterationInfo truncate ();
  
     DKTSTruncation& generateRandomExample     (const LongInt& d, const LongInt& r, const LongInt& R, 
                                                const LongInt& n, const LongReal& eps=0.5);
     DKTSTruncation& setInitialGuestBiggestSum (const LongInt& r);
     
     DKTSTruncation& readInitialGuestFrom      (const IString& fileName);               
     DKTSTruncation& writeSolutionTo           (const IString& fileName);
     
     DKTSTruncation& writePreComputationDataTo (const IString& fileNameCt=IString("coefficientSystem.ten"), 
                                                const IString& fileNameOb=IString("orthogonalBasis.ten"));     
     DKTSTruncation& writeCoefficientSystemTo  (const IString& fileName=IString("coefficientSystem.ten"));
     DKTSTruncation& writeOrthogonalBasisTo    (const IString& fileName=IString("orthogonalBasis.ten"));

     DKTSTruncation& prepareInputData    (const DKTS& A);
     DKTSTruncation& prepareInputData    (const DKTS& A, const LongInt& r);
     DKTSTruncation& prepareInputData    (const DKTS& A, const DKTS& X);
     
     bool            writeParameter      (ostream& s) const;
     DescriptionData parameterDescription()           const;

     bool   writeNormsOfA (DescriptionData& dDIn) const;
     bool   writeNormsOfA (ostream& s)            const;     

     inline DKTS& a() { return attr_a; }
     inline DKTS& x() { return attr_x; }
     
     inline DKTS& Z() { return attr_Z; }
     
     inline DGKTSDecomposer& decomposer() { return attr_decomposer;}

     DKTSDIterationInfo truncate2      (const LongInt& rT, DKTS& a, DKTS& x, const LongReal& normA);
     
     DKTSDIterationInfo truncate2Eps   (const LongReal& eps, DKTS& a, DKTS& x, const bool& plotInfo=true, const bool& useX=false);
     DKTSDIterationInfo truncate2EpsSA (const LongReal& eps, DKTS& a, DKTS& x, const bool& plotInfo=true);
     DKTSDIterationInfo truncate2EA    (DKTS& a, DKTS& x, const bool& plotInfo=true, const LongInt& steps=1);

  //@}


  /*!
  ***************************************************************************************
  * \name                          Implementierungsbereich
  ***************************************************************************************/

  //@{
  
     DKTS attr_a;     
     DKTS attr_x;
     
     DKTS attr_Z;

  private:       
  
     DKTSDIterationInfo bestR1Truncation (const LongInt& r, DKTS& a, DKTS& x, const LongReal& normD);     
     DKTSDIterationInfo startTruncation  (DKTS& a, DKTS& x, const LongReal& normA);         
     
     DKTSTruncation& computeMainPart     (const LongReal& eps);  
  
     DKTSTruncation& resize              (const LongInt& d, const LongInt& r, const LongInt& R, const LongInt& n);
     DKTSTruncation& allocateDataSpace   (const LongInt& d, const LongInt& r, const LongInt& R, const LongInt& n);
     DKTSTruncation& deleteDataSpace     ();
     DKTSTruncation& resizeInitialGuest  (const LongInt& d, const LongInt& r, const LongInt& R, const LongInt& n);
  
     DKTSTruncation& sortIndexSet ();
     LongInt         partition    (LongInt low, LongInt high, LongRealPointer f);
     DKTSTruncation& quickSort    (LongInt low, LongInt high, LongRealPointer f);
     DKTSTruncation& swapIndex    (const LongInt& i, const LongInt& j);    

     LongRealPointer attr_values;
     LongIntPointer  attr_indexSet;

     DKTSTruncation& addInputTensorInformation (const IString& date, Protocol& truncationLog);
     DKTSTruncation& addBeginTrancation        (const LongInt& r, Protocol& truncationLog);
     DKTSTruncation& addBeginTrancationR1      (const LongInt& r, Protocol& truncationLog);
     DKTSTruncation& addTrancationR1           (const LongInt& run ,const LongInt& rang, Protocol& truncationLog);
     DKTSTruncation& addDataBlock              (const DKTSDDataBlock& dataBlock, const DKTSDIterationInfo& infoBlock, 
                                                const LongInt& k, Protocol& truncationLog);
     DKTSTruncation& addEndTrancation          (const DKTSDIterationInfo& infoBlock, const LongInt& r,
                                                Protocol& truncationLog);                                                                                                

     DKTSTruncation& addInfoBlock              (const DKTSDIterationInfoBlock& infoBlock, const LongReal& totalTime,
                                                const LongReal& eps, const LongReal& epsN, const LongReal& preC);
     DKTSTruncation& addInfoBlockR1            (const DKTSDIterationInfoBlock& infoBlock);

//     DGKTSD attr_decomposer;
     DGKTSDNewton attr_decomposer;
     //DGKTSDALS attr_decomposer;
    
				 DKTSTruncation& setDefaultParameter();
				
    //! 

  //@}
 };


typedef  DKTSTruncation* DKTSTruncationPointer;

#endif // not defined __DKTSTruncation__
