/*
 * Copyright (C) 2011 Mike Espig
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

// DGKTSDNSOR2.cpp: Implementierung der Klasse DGKTSDNSOR2.
//
//////////////////////////////////////////////////////////////////////

#include "DGKTSDNSOR2.hpp"
#include "BlasInterface.hpp"


DGKTSDNSOR2::DGKTSDNSOR2()
:DGKTSDDataBaseNSOR()
 {
 }


DGKTSDNSOR2::~DGKTSDNSOR2()
 {
 }


DKTSDIterationInfo DGKTSDNSOR2::truncate2Eps(const LongReal& eps, const DKTS& A, DKTS& X, const bool& plotInfo, const bool& useX)
 {  
    DKTS Z, a, x;

    Z.setOrthogonal2(A, 1.0e-6);
    a.setCoefficientsSystemOf(A, Z);
    x.setCoefficientsSystemOf(X, Z);

    const LongInt d = a.d();
    const LongInt R = a.k();
    const LongInt n = a.n();

    if(plotInfo)
     {
        cout << endl;
        cout << "d = " << d << endl;
        cout << "R = " << R << endl;
        cout << "t = " << n << endl; 
        cout << "n = " << A.n() << endl;        
     }    
    
    const LongReal normA = frobeniusNorm(a);
            
    DKTSDIterationInfo infoEntry;            
    DKTSDDataBlock dataBlock;    
    
    LongReal error = 1.0;
    
    if(useX==false)
     {        
        x.setAlsApproximationOf(a, 10);
        infoEntry = decompose(x, a, normA, dataBlock);        
        
        if(plotInfo)
         {
            cout << infoEntry << endl;            
         }        
     }
    else
     {	       
        infoEntry = decompose(x, a, normA, dataBlock);            
                
        if(plotInfo)
         {
            cout << infoEntry << endl;
         }
         
        //addDataBlock(dataBlock, infoEntry, x.k(), attr_truncationLog);     
     }            
    
    DKTS x0(d, 1, n), residuum(d, R, n);
        
    LongInt r = x.k();        
    
    error = infoEntry.error();
    
    while(eps<error && r<R)
     {     
        residuum.setSumOf(1.0, a, -1.0, x);        
        x0.setAlsApproximationOf(residuum, 10);

       	DKTS xt(x);
	       x.setSumOf(1.0, xt, 1.0, x0);	       
	               
        infoEntry = decompose(x, a, normA, dataBlock);

        error = infoEntry.error();                

        if(plotInfo)
         {
            cout << infoEntry << endl;
         }

//        addDataBlock(dataBlock, infoEntry, x.k(), attr_truncationLog);

        r++;
     }
    
    X.regeneratedBy(x, Z);
    
   return infoEntry; 
 }


DKTSDIterationInfo DGKTSDNSOR2::decompose(DKTS& x, const DKTS& a, const LongReal& normA, DKTSDDataBlock& dataBlock)
 {
    DKTSDIterationInfo infoEntry;
            
    if(resize(x, a)==true)
     {
        (*this)
         .setNormA2(normA*normA)
        ;

        infoEntry = startIteration(x, a);
     }
    else
     {
        throw SimpleException(IString("Warning In DKTSDIterationInfo DGKTSDNSOR2::decompose(DKTS& a, DKTS& xi, const LongReal& normA, DKTSDDataBlock& dataBlock), Error in resize !!!"));
     }
 
   return infoEntry;
 }


DKTSDIterationInfo DGKTSDNSOR2::startIteration(DKTS& x, const DKTS& a)
 {
    Timer timer;
    timer.startTiming();
 
    DKTSDIterationInfo infoEntry;
     
    const LongInt maxSteps   = DGKTSDNSOR2::maxIterations(); 

    LongInt numberOfIterations = 1;
    
    const LongInt d      = x.d();
    const LongInt kRankr = x.k();
    const LongInt kRankR = a.k();
    
    if(kRankr==kRankR)
     {
        numberOfIterations = maxSteps+1;
        x = a;
     }

    DKTVector& grad    = attr_gradient();
    DKTVector& vecX    = x(); 

    (*this)
     .computeAllInnerProducts(x, a)
     .computeGradient(x, a)     
    ;

    LongReal errorSq    = innerProduct(grad, grad);
    LongReal errorSqOld = errorSq;

    LongReal distance    = sqrt(fabs(1.0+2.0*f()));
    LongReal distanceOld = distance;
    LongReal decr        = 1.0; 
    
    const LongReal dOld = distance;   

    attr_error = sqrt(errorSq);        
    
    Timer time;    
    LongInt mu1, j1;
    
    while((attr_epsilon<attr_error || attr_accuracy <= decr) && numberOfIterations <= maxSteps)    
     {  
        computeMaxPivot(j1, mu1);
        
        /*
        cout << "numberOfIterations = " << numberOfIterations << endl;
        cout << "j1  = " << j1 << endl;
        cout << "mu1 = " << mu1 << endl;
        */        
        
        DKTVector& x_jmu = x(j1, mu1);

        x_jmu.add(-attr_alpha/attr_xx.tensorInnerProductAt(j1, j1, mu1), attr_gradient(j1, mu1), 1.0);       
                
        (*this)
         .computeAllInnerProducts(x, a)
         .computeGradient(x, a)
        ;                  

        errorSqOld  = errorSq;
        errorSq     = innerProduct(grad, grad);
        attr_error  = sqrt(errorSq);

        distanceOld = distance;
        distance    = sqrt(fabs(1.0 + 2.0*f()));

        decr        = (distanceOld-distance)/distance;                                                      
        
/*
        cout << setprecision(4) << scientific;
        cout << numberOfIterations << ", " << attr_error << ", " << distance << endl;
*/         
        numberOfIterations++;
     }
    
    numberOfIterations--;
    
    const LongReal diff = dOld-distance;
    const LongReal sec  = timer.elapsedTimeSec();
    
    if(kRankr==kRankR)
     {
        numberOfIterations = 0;
     }
      
    infoEntry
     .setStartError(dOld)
     .setNumberOfNewtonSteps(numberOfIterations)
     .setRelativeDifferenz(diff/dOld*100)
     .setGradient(attr_error)
     .setStep(kRankr)
     .setError(distance)
    ;
    
				infoEntry.setGradient(attr_error);
				infoEntry.setCalculationTime(sec);

   return infoEntry;
 }


DGKTSDNSOR2& DGKTSDNSOR2::computeMaxPivot (LongInt& j1, LongInt& mu1)
 {
    const LongInt d = DGKTSDNSOR2::d();
    const LongInt r = DGKTSDNSOR2::r();
    
    LongReal maxValue     = -1.0;
    LongReal currentValue = 1.0;
    
    
    for(LongInt mu=0; mu<d; mu++)
     {
        for(LongInt j=0; j<r; j++)
         {
            DKTVector& grad = attr_gradient(j, mu);
            
            currentValue = innerProduct(grad, grad);
            
            if(maxValue < currentValue)
             {
                j1  = j;
                mu1 = mu;
             }
         }
     }
 
   return (*this);
 }
 

LongReal DGKTSDNSOR2::f() const
 {     
     const LongInt r = DGKTSDNSOR2::r();
     const LongInt R = DGKTSDNSOR2::k();     

     LongReal value = 0.0;
     
     for(LongInt j1=0; j1<r; j1++)
      {
         value += 0.5 * attr_xx.tensorInnerProductAt(j1, j1, 0) * attr_xx.innerProductAt(j1, j1, 0);
         
         for(LongInt i=0; i<R; i++)
          {
             value -= attr_xa.tensorInnerProductAt(j1, i, 0) * attr_xa.innerProductAt(j1, i, 0);
          }
                           
         for(LongInt j2=0; j2<j1; j2++)
          {
             value += attr_xx.tensorInnerProductAt(j1, j2, 0) * attr_xx.innerProductAt(j1, j2, 0);
          }
      }     
     
     value *= 1.0/normA2();
     
    return value;
 }



DGKTSDNSOR2& DGKTSDNSOR2::computeGradient(const DKTS& x, const DKTS& a)
 {
    const LongInt d = DGKTSDNSOR2::d();
    const LongInt m = DGKTSDNSOR2::m(); 
    
    const LongInt r = DGKTSDNSOR2::r(); // tRank of x
    const LongInt k = DGKTSDNSOR2::k(); // tRank of b    

    const LongReal normA2 = DGKTSDNSOR2::normA2();
    
    const LongReal factor1 = -1.0/normA2;
    const LongReal factor2 = -factor1;
    
    const LongReal const_eins = 1.0;
    const LongReal const_null = 0.0;
    
    const LongInt const_inc  = 1;
    char  const_notConjTrans =  'n';    
    
    for(LongInt mu1=0; mu1<d; mu1++)
     {
        for(LongInt j1=0; j1<r; j1++)
         {
            DKTVector& gradV = attr_gradient(j1, mu1);
            LongReal&  grad  = gradV(0);            
            
            TensorCalculus::Blas<double>::gemv (const_notConjTrans, m, k, factor1, &a(mu1)(0), m, &attr_xa(j1, mu1), const_inc,
                   const_null, &grad, const_inc);

            TensorCalculus::Blas<double>::gemv (const_notConjTrans, m, r, factor2, &x(mu1)(0), m, &attr_xx(j1, mu1), const_inc,
                   const_eins, &grad, const_inc);
         }     
     }

   return (*this);
 }


DGKTSDNSOR2& DGKTSDNSOR2::computeAllInnerProducts(const DKTS& x, const DKTS& a)
 {
    attr_xa.computeAllInnerProducts(x, a);    
    attr_xx.computeAllInnerProducts(x, x);
        
   return (*this);
 }


