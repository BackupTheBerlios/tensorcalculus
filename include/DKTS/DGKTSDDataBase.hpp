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
 * \class    DGKTSDDataBase
 *
 * \brief    DGKTSDDataBase 
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

#ifndef __DGKTSDDataBase__
#define __DGKTSDDataBase__


#include "Macros.h"
#include <iostream>
#include <iomanip>

#include "DKTS.hpp"
#include "Protocol.hpp"


using namespace std;

class  DGKTSDDataBase
 {
    DECLARE (DGKTSDDataBase) 
 
    friend class DGKTSDFullMethod;
    friend class DGKTSDecomposer;
    friend class DGKTSD;
    friend class DGKTSDNewton;
 
    /*!
     The Order of the Kronecker-Tensor-Sum*/
    ATTRIBUTE (LongInt, d, setD)

    /*!
     The small Kronecker-Rank of the Kronecker-Tensor-Sum*/
    ATTRIBUTE (LongInt, k, setK)

    /*!
     The big Kronecker-Rank of the Kronecker-Tensor-Sum*/
    ATTRIBUTE (LongInt, l, setL)

    /*!
     The dimension of the image*/
    ATTRIBUTE (LongInt, m, setM)
    
    ATTRIBUTE (LongInt, preC, setPreC)

    /*!
     Parameter for the Newton-iteration*/   
    ATTRIBUTE (LongReal, epsilon,       setEpsilon)
    ATTRIBUTE (LongInt,  maxSteps,      setMaxSteps)    

    /*!
     Parameter for the cr-procedure*/
    ATTRIBUTE (LongReal, epsilonCR,     setEpsilonCR)
    ATTRIBUTE (LongInt,  maxStepsCR,    setMaxStepsCR)    
    ATTRIBUTE (LongInt,  maxStepsCRL1,  setMaxStepsCRL1)    
    ATTRIBUTE (LongInt,  maxStepsCRL2,  setMaxStepsCRL2)    

    ATTRIBUTE (LongInt,  quality, setQuality)
    
    ATTRIBUTE (bool, useNewtonMethode, setUseNewtonMethode)
    
    /*!
     Parameter for the Armijo-rule*/
    ATTRIBUTE (LongInt,  maxStepsAmijo, setMaxStepsAmijo)
    ATTRIBUTE (LongReal, sigma,         setSigma)
        
    ATTRIBUTE (LongReal, precision,     setPrecision)

    /*!
     Parameter for the target function*/

    ATTRIBUTE (LongReal, mu,     setMu);    
    ATTRIBUTE (LongReal, kappa,  setKappa);
    
    ATTRIBUTE (LongReal, gamma,  setGamma);    
    ATTRIBUTE (LongReal, delta,  setDelta);
        
    ATTRIBUTE (LongReal, kappa1, setKappa1);
    ATTRIBUTE (LongReal, kappa2, setKappa2);
    ATTRIBUTE (LongReal, kappa3, setKappa3);
    
    ATTRIBUTE (LongReal, lambda1, setLambda1);
    ATTRIBUTE (LongReal, lambda2, setLambda2);

    
    ATTRIBUTE (LongInt,  memory,        setMemory)

    ATTRIBUTE (LongReal, upperBound,    setUpperBound)

    ATTRIBUTE (LongReal, optW,          setOptW)
    
    ATTRIBUTE (bool,     printCout,     setPrintCout)
    
  /*!
  ***************************************************************************************
  * \name                        Konstruktion / Destruktion
  ***************************************************************************************/

  //@{

  public:
   
     DGKTSDDataBase ();

     virtual ~DGKTSDDataBase();
  

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

    bool               resize            (const DKTS& a, DKTS& xi);
    DGKTSDDataBase&    allocateWorkSpace (const LongInt& d, const LongInt& k, const LongInt& l, const LongInt& m);
    DGKTSDDataBase&    deleteWorkSpace   ();    

    LongReal& A   (const LongInt& mu) const;
    LongReal& A   (const LongInt& j1, const LongInt& j2, const LongInt& mu)  const;
    LongReal& W   (const LongInt& j,  const LongInt& mu) const;

    LongReal& ad  (const LongInt& j,  const LongInt& i,  const LongInt& mu)  const;
    LongReal& dd  (const LongInt& j1, const LongInt& j2, const LongInt& mu)  const;
    LongReal& xd  (const LongInt& j1, const LongInt& j2, const LongInt& mu)  const;

    
    LongReal& x   (const LongInt& j1, const LongInt& j2, const LongInt& mu)  const;
    LongReal& a   (const LongInt& j,  const LongInt& i,  const LongInt& mu)  const;
    LongReal& phi (const LongInt& j,  const LongInt& i,  const LongInt& mu1, const LongInt& mu2) const;
    LongReal& psi (const LongInt& j1, const LongInt& j2, const LongInt& mu1, const LongInt& mu2) const;
    
    
    LongInt isIndexOffsetOfx   (const LongInt& j1,  const LongInt& j2)  const;
    LongInt isIndexOffsetOfa   (const LongInt&  j,  const LongInt& i)   const;
    LongInt isIndexOffsetOfxd  (const LongInt& j1,  const LongInt& j2)  const;
    LongInt isIndexOffsetOfmuA (const LongInt& mu1, const LongInt& mu2) const;
    LongInt isIndexOffsetOfmuX (const LongInt& mu1, const LongInt& mu2) const;

           
    DGKTSDDataBase& setDefaultParameter();

    //WorkSpace
    LongRealPointer attr_a;
    LongRealPointer attr_x;
 
    LongRealPointer attr_xd;
    LongRealPointer attr_ad;
    LongRealPointer attr_dd;

    LongRealPointer attr_phi;
    LongRealPointer attr_psi;
        
    LongInt         attr_dl;
        
    // newton    
    LongRealPointer attr_workK;
    LongRealPointer attr_workL;

    LongReal        attr_error;

    LongRealPointer attr_A;
    LongRealPointer attr_W;

    LongRealPointer attr_eigWork;
    LongRealPointer attr_eigA;
    LongRealPointer attr_eigValue;

    DKTS attr_gradient;
    DKTS attr_gradientV;
    DKTS attr_direction;
    DKTS attr_work;


    // cr
    DKTS attr_w;
    DKTS attr_r;

    DKTS attr_ap;
    DKTS attr_apo;

    DKTS attr_b;
    DKTS attr_bo;

    DKTS attr_p;
    DKTS attr_po;

    //!

  //@}
 };


typedef  DGKTSDDataBase* DGKTSDDataBasePointer;

#endif // not defined __DGKTSDDataBase__
