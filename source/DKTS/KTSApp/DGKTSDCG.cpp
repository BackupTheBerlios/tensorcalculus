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

// DGKTSDCG.cpp: Implementierung der Klasse DGKTSDCG.
//
//////////////////////////////////////////////////////////////////////

#include "DGKTSDCG.hpp"
#include "LapackInterface2.hpp"
#include "BlasInterface.hpp"


DGKTSDCG::DGKTSDCG()
:DGKTSDDataBaseCG()
 {
 }


DGKTSDCG::~DGKTSDCG()
 {
 }


DKTSDIterationInfo DGKTSDCG::truncate2Eps(const LongReal& eps, const DKTS& A, DKTS& X, const bool& plotInfo, const bool& useX)
 {  
    DKTS Z, a, x;

    Z.setOrthogonal2(A, 1.0e-6);
    a.setCoefficientsSystemOf(A, Z);
    
    //a.setCopyOf(A);
    
    if(useX==true)
     {
        x.setCoefficientsSystemOf(X, Z);
     }    


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
	       x.improveApproximation(a, 20, 5);
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
        //decompose(x0, residuum, error*normA, dataBlock);    

       	DKTS xt(x);
	       x.setSumOf(1.0, xt, 1.0, x0);
	       
	       x.improveApproximation(a, 10, 5);	   
	               
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
    //X.setCopyOf(x);
    
   return infoEntry; 
 }


DKTSDIterationInfo DGKTSDCG::decompose(DKTS& x, const DKTS& a, const LongReal& normA, DKTSDDataBlock& dataBlock)
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
        throw SimpleException(IString("Warning In DKTSDIterationInfo DGKTSDCG::decompose(DKTS& a, DKTS& xi, const LongReal& normA, DKTSDDataBlock& dataBlock), Error in resize !!!"));
     }
 
   return infoEntry;
 }


DKTSDIterationInfo DGKTSDCG::startIteration(DKTS& x, const DKTS& a)
 {
    Timer timer;
    timer.startTiming();
 
    DKTSDIterationInfo infoEntry;
     
    const LongInt maxSteps   = DGKTSDCG::maxIterations(); 

    LongInt numberOfIterations = 1;
    
    const LongInt kRankr = x.k();
    const LongInt kRankR = a.k();
    
    if(kRankr==kRankR)
     {
        numberOfIterations = maxSteps+1;
        x = a;
     }

    DKTVector& grad    = attr_gradient();
    DKTVector& gradOld = attr_gradientOld();
    DKTVector& d       = attr_direction();
    DKTVector& w       = attr_work();
    DKTVector& vecX    = x(); 

    (*this)
     .computeAllInnerProducts(x, a)
     .computeGradient(x, a)     
    ;

    gradOld = grad;

/*
    // d = -grad
    d = grad;
    d *= -1.0;    
*/

    (*this)
     .computeAinv()
     .evaluateAinv(-1.0, attr_gradient, attr_direction)
    ;

    d.add(-0.5, grad, 0.5);


    LongReal errorSq    = innerProduct(grad, grad);
    LongReal errorSqOld = errorSq;

    LongReal distance    = sqrt(fabs(1.0+2.0*f()));
    LongReal distanceOld = distance;
    LongReal decr        = 1.0; 
    
    const LongReal dOld = distance;   

    attr_error = sqrt(errorSq);        
    
    while((attr_epsilon<attr_error || attr_accuracy <= decr) && numberOfIterations <= maxSteps)    
     {         
         const LongReal alpha = computeStepSize(x, a);                  
         
         vecX.add(alpha, d, 1.0);
         
         gradOld = grad;         
         
         Timer time;
         
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
         
         
         LongReal beta = MAX(0.0, (errorSq-innerProduct(grad,gradOld))/errorSqOld); //(2)
         
         (*this)
          .computeDirection(x, a)
         ;
         
         w.add(-0.5, grad, 0.5);
                          
         d.add(-1.0, w, beta);         

//         d.add(-1.0, grad, beta);
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
 
LongReal DGKTSDCG::derivativeOfPhi(const LongReal& t, const DKTS& x, const DKTS& a)
 {
    
    (*this)
     .computeAo(t)
     .computeAu(t)
     .computeXo(t)
     .computeXu(t)
    ;
        
    LongReal value = 0.0;
    
    const LongInt d = DGKTSDCG::d();
    const LongInt m = DGKTSDCG::m(); 
    
    const LongInt r = DGKTSDCG::r();
    const LongInt k = DGKTSDCG::k();    
    
    const LongReal normA2 = DGKTSDCG::normA2();
    
    const LongReal factor1 = -1.0/normA2;
    const LongReal factor2 = -factor1;
    const LongReal factor3 = t*factor2;
    
    const LongReal const_eins = 1.0;
    const LongReal const_null = 0.0;
    
    const LongInt const_inc  = 1;
    char  const_notConjTrans =  'n';    
        
    // computing f'(x+td)
    for(LongInt mu1=0; mu1<d; mu1++)
     {
        for(LongInt j1=0; j1<r; j1++)
         {
            DKTVector& gradV = attr_gradientS(j1, mu1);
            LongReal&  grad  = gradV(0);
                        
            for(LongInt i=0; i<k; i++)
             {
                attr_v1(i) = Au(i,j1,mu1) * Ao(i,j1,mu1);
             }             
            
            TensorCalculus::Blas<double>::gemv (const_notConjTrans, m, k, factor1, &a(mu1)(0), m, &attr_v1(0), const_inc,
                   const_null, &grad, const_inc);

            for(LongInt j=0; j<r; j++)
             {
                attr_v2(j) = Xo(j,j1,mu1) * Xu(j,j1,mu1);                
             }             

            TensorCalculus::Blas<double>::gemv (const_notConjTrans, m, r, factor2, &x(mu1)(0), m, &attr_v2(0), const_inc,
                   const_eins, &grad, const_inc);
            
            TensorCalculus::Blas<double>::gemv (const_notConjTrans, m, r, factor3, &attr_direction(mu1)(0), m, &attr_v2(0), const_inc,
                   const_eins, &grad, const_inc);
         }     
     }
    //cout<<attr_gradientS<<endl;
    //edit 
    value = innerProduct(attr_direction(), attr_gradientS())*scalingFactor(); 
    
    
   return value;
 }


LongReal DGKTSDCG::f() const
 {     
     const LongInt r = DGKTSDCG::r();
     const LongInt R = DGKTSDCG::k();     

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


DGKTSDCG& DGKTSDCG::setDirection(const DKTS& d)
 {
    attr_direction = d;
            
   return (*this);
 }


DGKTSDCG& DGKTSDCG::computeGradient(const DKTS& x, const DKTS& a)
 {
    const LongInt d = DGKTSDCG::d();
    const LongInt m = DGKTSDCG::m(); 
    
    const LongInt r = DGKTSDCG::r(); // tRank of x
    const LongInt k = DGKTSDCG::k(); // tRank of b    

    const LongReal normA2 = DGKTSDCG::normA2();
    
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


DGKTSDCG& DGKTSDCG::computeAllInnerProducts(const DKTS& x, const DKTS& a)
 {
    attr_xa.computeAllInnerProducts(x, a);    
    attr_xx.computeAllInnerProducts(x, x);
        
   return (*this);
 }

DGKTSDCG& DGKTSDCG::preComputeValuesForStepSizeRule(const DKTS& x, const DKTS& a)
 {
    const LongInt d = DGKTSDCG::d();
  
    LongReal temp          = fabs(innerProduct(attr_direction(), attr_gradient()));
    //cout<<"scaling^d : "<<temp<<endl;
    LongReal scalingFactor = 1.0 / pow(temp, (1.0/(LongReal)d)); //scalingFactor muss positiv sein!
    
    (*this)
     .setScalingFactor(scalingFactor)
    ;
    
    attr_da.computeInnerProducts(attr_direction, a);
    attr_xd.computeInnerProducts(x, attr_direction);
    attr_dd.computeInnerProducts(attr_direction, attr_direction);
    
    attr_da *= scalingFactor;
    attr_xd *= scalingFactor;
    attr_dd *= scalingFactor;
    attr_xx *= scalingFactor;
    attr_xa *= scalingFactor; 
    
   return (*this);
 }

 
LongReal DGKTSDCG::computeStepSize(const DKTS& x, const DKTS& a)
 { 
    (*this)
     .preComputeValuesForStepSizeRule(x, a)
    ;    
    
    LongReal value = valueFromCurryRule(x, a);
            
   return value;
 }
 
 
LongReal DGKTSDCG::valueFromCurryRule(const DKTS& x, const DKTS& a) 
 {
    LongReal value = 0.0;

    LongReal aa =  0.0;
    LongReal bb =  2.0*alpha();
    
    //edit
    //LongReal fa = innerProduct(attr_gradient(), attr_direction())/scalingFactorOld(); //edit kann man -1 setzen?
    LongReal fa = -1.0;
    LongReal fb = derivativeOfPhi(bb, x, a);
    
    // cout << "fa = " << fa << ", aa = " << aa << endl;
    // cout << "fb = " << fb << ", bb = " << bb << endl;
    
    while(fa*fb>0.0)
     {
        aa = bb;
        fa = fb;
        
        bb *= 1.0e2;
        fb = derivativeOfPhi(bb, x, a);
       // cout << "fb = " << fb << ", bb = " << bb << endl;        
     }
    
        
    LongReal R   = aa;
    LongReal fba = (fb-fa)/(bb-aa);
    LongReal cc  = bb - fb/fba;
    LongReal fc  = derivativeOfPhi(cc, x, a);
    
    const LongReal C = DGKTSDCG::C();
    const LongReal D = DGKTSDCG::D();   
    
    const LongInt  maxSteps = maxIterationsCR();
    const LongReal epsilon  = epsilonStepSize();

    LongInt numberOfIterations = 1;
    
    LongReal error = fabs(fc);
    LongReal y = 1.0;        
    while(epsilon<error && numberOfIterations<=maxSteps)
     {
        if(fc*fb<0.0)
         {
            R = bb;
         }   
         
        LongReal fcb  = (fc-fb)/(cc-bb);
        LongReal fca  = (fc-fb)/(cc-bb);
            
        LongReal Qabc = ((cc-aa)*fcb + (bb-cc)*fca)/(bb-aa);
            
        if(fabs(Qabc)<EPS_NULL)
         {
            y = 0.5*(R + cc);
         }
        else
         {
            y = cc - fc/Qabc;
         }
        if((y-cc)*(y-R)<0.0 && (fabs(y-R)<C*fabs(cc-R) || fabs(fc)<D*fabs(fb)))
         {
            // y is a good value
         }
        else
         {
            y = 0.5*(R + cc);
         }               
        
        aa = bb;
        fa = fb;
        bb = cc;
        fb = fc;
        cc = y;
        fc = derivativeOfPhi(cc, x, a);
        
/*        
        cout << "i = " << numberOfIterations << ", ";
        cout << "t_i = " << cc << ", ";
        cout << "f(t_i) = " << fabs(fc) << endl;
*/              
        error = fabs(fc);
        numberOfIterations++;        
     }
/*
    numberOfIterations--;
    cout << "i = " << numberOfIterations << ", ";
    cout << "t_i = " << cc << ", ";
    cout << "f(t_i) = " << fabs(fc) << endl;
*/
    value = cc;
    
    (*this)
     .setAlpha(value)
    ;
     
   return value;
 }
 
 
bool DGKTSDCG::evaluateA(const LongReal& alpha, const DKTS& v, DKTS& w) const
 {
    bool value = true;
    
    const LongInt d = v.d();
    const LongInt r = v.k();
    
    w().setNull();
    
    for(LongInt mu1=0; mu1<d; mu1++)
     {
        for(LongInt j1=0; j1<r; j1++)
         {
            DKTVector& w_jmu = w(j1, mu1);
            
            // The Matrix A
            for(LongInt j2=0; j2<r; j2++)
             {
                w_jmu.update(alpha*attr_xx.tensorInnerProductAt(j1, j2, mu1), v(j2, mu1));
             }            
         }
     }
    
   return value;
 }


DGKTSDCG& DGKTSDCG::computeDirection(const DKTS& x, const DKTS& a)
 { 
    (*this)
     .computeAinv()
     .evaluateAinv(1.0, attr_gradient, attr_work)
    ;
     
   return (*this);
 }


DGKTSDCG& DGKTSDCG::computeAinv()
 {
    LongInt d = DGKTSDCG::d();
    LongInt r = DGKTSDCG::r();
    
    LongInt info = 0;
    char const_u = 'U';
    
    for(LongInt mu=0; mu<d; mu++)
     {
        for(LongInt j1=0; j1<r; j1++)
         {
            for(LongInt j2=j1; j2<r; j2++)
             {
                const LongReal value =  attr_xx.tensorInnerProductAt(j1, j2, mu);
                
                Ainv(j1, j2, mu) = value;                
                Ainv(j2, j1, mu) = value;
             }
         }
         
        info = TensorCalculus::Lapack<double>::potrf(const_u, r, &Ainv(0, 0, mu), r);
        info = TensorCalculus::Lapack<double>::potri(const_u, r, &Ainv(0, 0, mu), r);
     }
    
   return (*this);
 }

 
bool DGKTSDCG::evaluateAinv(const LongReal& alpha, const DKTS& v, DKTS& w) const
 {
    bool value = true;
    
    const LongInt d = v.d();
    const LongInt r = v.k();
    
    w().setNull();
    
    for(LongInt mu1=0; mu1<d; mu1++)
     {
        for(LongInt j1=0; j1<r; j1++)
         {
            DKTVector& w_jmu = w(j1, mu1);
            
            // The Matrix A^-1
            for(LongInt j2=0; j2<r; j2++)
             {
                w_jmu.update(alpha*Ainv(j1, j2, mu1), v(j2, mu1));
             }            
         }
     }
    
   return value;
 } 


DGKTSDCG& DGKTSDCG::computeAu (const LongReal& t)
 {
    LongInt k = DGKTSDDataBaseCG::k();
    LongInt r = DGKTSDDataBaseCG::r();
    LongInt d = DGKTSDDataBaseCG::d();
    LongInt nu;
    
    for(LongInt mu=1; mu<d; mu++)
     {
        for(LongInt j=0; j<r; j++)
         {
             for(LongInt i=0; i<k; i++)
             {
                if(mu == 1)
                 {
                    (*this).Au(i, j, 0)=1;
                 }
                nu = mu-1;
                (*this).Au(i, j, mu)  = (*this).Au(i, j, nu) * (attr_xa.innerProductAt(j, i, nu) + t*attr_da.innerProductAt(j, i, nu));
                //cout<<"Au  "<<(*this).Au(i,j,mu)<<"   ";
             }
         }
     }
    //cout<<endl;
   return (*this);
 }
   
         
DGKTSDCG& DGKTSDCG::computeAo (const LongReal& t)
{
    LongInt k = DGKTSDDataBaseCG::k();
    LongInt r = DGKTSDDataBaseCG::r();
    LongInt d = DGKTSDDataBaseCG::d();
    LongInt nu;
    
    //for Schleife R�ckw�rts
    for(LongInt mu=d-2; 0<=mu; mu--)
     {
        for(LongInt j=0; j<r; j++)
         {   
            for(LongInt i=0; i<k; i++)
             {
                if(mu == d-2)
                 {
                    (*this).Ao(i, j, d-1) = 1;
                 }
                nu = mu+1;
                (*this).Ao(i, j, mu)  = (*this).Ao(i, j, nu) * (attr_xa.innerProductAt(j, i, nu) + t*attr_da.innerProductAt(j, i, nu));
             }
         }
     }
   return (*this);
 } 
 
         
DGKTSDCG& DGKTSDCG::computeXu (const LongReal& t)
 {  
    LongInt r = DGKTSDDataBaseCG::r();
    LongInt d = DGKTSDDataBaseCG::d();
    LongInt nu;
    LongReal temp;
    
//    cout<<"t = "<<t<<endl;
    for(LongInt mu=1; mu<d; mu++)
     {
        for(LongInt j=0; j<r; j++)
         {
            for(LongInt i=0; i<r; i++)
             {
                if(mu==1)
                 {  
                    (*this).Xu(i, j, 0)=1;
                 }
                nu = mu-1;
                
                temp  = attr_xx.innerProductAt(j, i, nu);
                temp += t*(attr_xd.innerProductAt(j, i, nu) + attr_xd.innerProductAt(i, j, nu));
                temp += t*t*attr_dd.innerProductAt(j, i, nu);
                
                (*this).Xu(i, j, mu)  = (*this).Xu(i, j, nu) * temp;
                //cout<<"("<<i<<", "<<j<<", "<<mu<<") : "<<(*this).Xu(i,j,mu)<<endl;
             }
            
         }
     }
   return (*this);
 }
 
                   
DGKTSDCG& DGKTSDCG::computeXo (const LongReal& t)
 {
    LongInt r = DGKTSDDataBaseCG::r();
    LongInt d = DGKTSDDataBaseCG::d();
    LongInt nu;
    LongReal temp;
    
    //for Schleife R�ckw�rts 
    for(LongInt mu=d-2; 0<=mu; mu--)
     {
        for(LongInt j=0; j<r; j++)
         {
            for(LongInt i=0; i<r; i++)
             {
                if(mu==d-2)
                 {
                    (*this).Xo(i, j, d-1) = 1;
                 }
                 
                nu = mu+1;
                
                temp  = attr_xx.innerProductAt(j, i, nu);
                temp += t*(attr_xd.innerProductAt(j, i, nu) + attr_xd.innerProductAt(i, j, nu));
                temp += t*t*attr_dd.innerProductAt(j, i, nu);
                
                (*this).Xo(i, j, mu)  = (*this).Xo(i, j, nu) * temp;
             }
         }
     }
   return (*this);
 }      
