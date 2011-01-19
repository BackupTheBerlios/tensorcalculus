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

// DGKTSDNewton.cpp: Implementierung der Klasse DGKTSDNewton.
//
//////////////////////////////////////////////////////////////////////

#include "DGKTSDNewton.hpp"
#include "BlasInterface.hpp"



DGKTSDNewton::DGKTSDNewton()
 {
 
 }


DGKTSDNewton::~DGKTSDNewton()
 {

 }


bool DGKTSDNewton::evaluateHf(const DKTS& a, const DKTS& xi, const DKTS& v, DKTS& w) const
 {
    bool value = true;

    LongInt d = DGKTSDNewton::d();
    LongInt k = DGKTSDNewton::k();
    LongInt l = DGKTSDNewton::l();
    LongInt m = DGKTSDNewton::m();

    LongReal facf1  = attr_kappa1;
    LongReal facg1  = attr_kappa2*attr_lambda1;
    LongReal facg2  = attr_kappa3*attr_lambda2;
    LongReal tempF  = 2.0*(LongReal)(d-1);

    w().setNull();

    for(LongInt mu1=0; mu1<d; mu1++)
     {
        for(LongInt j1=0; j1<k; j1++)
         {
            DKTVector&      w_jmu      = w(j1, mu1);
            const LongReal& preFactorL = W(j1, mu1);            

            for(LongInt mu2=0; mu2<d; mu2++)
             {
                if(mu2==mu1)
                 {
                    // G1
                    if(mu1==0)
                     {
                        w_jmu.update(attr_mu*facg1*tempF*preFactorL*innerProduct(xi(j1,mu1), v(j1, mu1)), xi(j1,mu1));
                     }
                    else
                     {
                        w_jmu.update(attr_mu*facg1*2.0*preFactorL*innerProduct(xi(j1,mu1), v(j1, mu1)), xi(j1,mu1));
                     }
                    // end 
                    
                    // The Matrix A and lambda2 G2(A)
                    for(LongInt j2=0; j2<k; j2++)
                     {                     
                        LongReal preFactor2  =  W(j2, mu2);
                         
                        if(j1==j2)
                         {
                            w_jmu.update((facf1 + attr_mu*facg2)*psi(j1, j2, mu2, mu2)*preFactor2, v(j2, mu2));
                         }
                        else
                         {
                            w_jmu.update(facf1*psi(j1, j2, mu2, mu2)*preFactor2, v(j2, mu2));
                         }
                     }
                 }
                else
                 {

                    for(LongInt j2=0; j2<k; j2++)
                     {
                        LongReal preFactorP  =  W(j2, mu2);
                        LongReal preFactor2  =  facf1*preFactorP;
                        
                        if(j1==j2)
                         {
                            
                            LongReal mPreFactor2 = -1.0*preFactor2;
                            
                            //The Matrix C and lambda2*G2(C)
                            w_jmu.update((facf1 + 2.0*facg2)*attr_mu*psi(j1, j2, mu1, mu2)*innerProduct(xi(j1, mu2), v(j2, mu2))*preFactorP, xi(j2, mu1));


                            if(EPS_NULL<kappa())
                             {                            
                                // The Matrix B = x^mu1 * d_j1,mu1,mu2 * (x^mu2)^t
                                LongReal& vjm = v(j2, mu2)(0);

                                TensorCalculus::Blas<double>::gemv(DGKTSDFullMethod::const_trans, m, k, preFactor2, &xi(mu2)(0), m, &vjm, const_inc,
                                        const_null, &attr_workK[0], const_inc);
                                   
                                for(LongInt j=0; j<k; j++)
                                 {
                                    attr_workK[j] *= (attr_mu*psi(j1, j, mu1, mu2));
                                 }
                             
                                TensorCalculus::Blas<double>::gemv(DGKTSDFullMethod::const_notConjTrans, m, k, const_eins, &xi(mu1)(0), m, &attr_workK[0], const_inc,
                                        const_eins, &w_jmu(0), const_inc);


                                // The Matrix D = a^mu1 * da_j1,mu1,mu2 * (a^mu2)^t
                                TensorCalculus::Blas<double>::gemv(DGKTSDFullMethod::const_trans, m, l, mPreFactor2, &a(mu2)(0), m, &vjm, const_inc,
                                        const_null, &attr_workL[0], const_inc);

                                for(LongInt i=0; i<l; i++)
                                 {
                                    attr_workL[i] *= (attr_mu*phi(j1, i, mu1, mu2));
                                 }

                                TensorCalculus::Blas<double>::gemv(DGKTSDFullMethod::const_notConjTrans, m, l, const_eins, &a(mu1)(0), m, &attr_workL[0], const_inc,
                                        const_eins, &w_jmu(0), const_inc);
                                       
                             }// End if(EPS_NULL<kappa())

                            // G1                               
                            if(mu1==0)
                             {
                                w_jmu.update(facg1*(-2.0)*attr_mu*preFactorP*innerProduct(xi(j1,mu2), v(j1, mu2)), xi(j1, mu1));
                             }
                            else if(mu2==0)
                             {
                                w_jmu.update(facg1*(-2.0)*attr_mu*preFactorP*innerProduct(xi(j1,mu2), v(j1, mu2)), xi(j1, mu1));
                             }                             
                         }
                        else //(j1 != j2)
                         {
                            //The Matrix C
                            w_jmu.update(preFactor2*attr_mu*psi(j1, j2, mu1, mu2)*innerProduct(xi(j1, mu2), v(j2, mu2)), xi(j2, mu1));                         
                         }
                     }
                 }                
             }

            w_jmu *= preFactorL;
         }// End for(LongInt j1=0; j1<k; j1++)
     }// End for(LongInt mu1=0; mu1<d; mu1++)    

   return value;
 }


bool DGKTSDNewton::evaluateA(const DKTS& a, const DKTS& xi, 
                       const DKTS& v, DKTS& w,
                       const LongReal& alpha, const LongReal& beta) const
 {
    bool value = true;
    LongInt d = DGKTSDNewton::d();
    LongInt k = DGKTSDNewton::k();

    w().setNull();

    for(LongInt mu1=0; mu1<d; mu1++)
     {
        for(LongInt j1=0; j1<k; j1++)
         {
            DKTVector& w_jmu = w(j1, mu1);
            const LongReal& preFactor = W(j1, mu1);

            // The Matrix A
            for(LongInt j2=0; j2<k; j2++)
             {
                w_jmu.update(psi(j1, j2, mu1, mu1)*W(j2, mu1), v(j2, mu1));
             }

            w_jmu *= preFactor*beta;
         }// End for(LongInt j1=0; j1<k; j1++)
     }// End for(LongInt mu1=0; mu1<d; mu1++)    

    w().update(alpha, v());

   return value;
 }


bool DGKTSDNewton::evaluateAinv(const DKTS& v,  DKTS& w) const
 {
    bool value = true;

    LongInt d = DGKTSDNewton::d();
    LongInt k = DGKTSDNewton::k();

    w().setNull();

    for(LongInt mu1=0; mu1<d; mu1++)
     {
        for(LongInt j1=0; j1<k; j1++)
         {
            DKTVector& w_jmu = w(j1, mu1);

            for(LongInt j2=0; j2<k; j2++)
             {
                w_jmu.update(A(j1, j2, mu1), v(j2, mu1));
             }
         }// End for(LongInt j1=0; j1<k; j1++)
     }// End for(LongInt mu1=0; mu1<d; mu1++)

   return value;
 }


DGKTSDNewton& DGKTSDNewton::solveDirection(const DKTS& a,  const DKTS& xi)
 {
    cout << setprecision(4);

    (*this)
     .useKrylovMethod(a, xi)
    ;

   return(*this);
 }


DGKTSDNewton& DGKTSDNewton::useKrylovMethod(const DKTS& a,  const DKTS& xi)
 {
    (*this)
     .invertA(0.0, 1.0)
    ;
            
    (*this)
     .setMu(MIN(attr_mu/attr_gamma,1.0))
    ;
    
    attr_direction() += xi();    

    const LongInt maxSteps   = maxStepsCR();
    const LongInt maxStepsL1 = maxStepsCRL1();
    const LongInt maxStepsL2 = maxStepsCRL2();
            
    if(cg(a, xi, maxStepsL1)==false)
     {
        (*this)
         .setMu(attr_mu*attr_gamma)
        ;
     
         while(cg(a, xi, maxStepsL2)==false)
          {
             (*this)
              .setMu(attr_mu*attr_gamma)
             ;
          }
     }

    evaluateW(attr_direction);

   return(*this);
 }


bool DGKTSDNewton::cg(const DKTS& ai, const DKTS& xi, const LongInt& maxSteps)
 {
    bool value = false;

    DKTVector& r     = attr_r();
    DKTVector& p     = attr_p();
    DKTVector& gradV = attr_gradientV();    
    DKTVector& x     = attr_direction();
    
    DKTVector& a     = attr_ap();    
    DKTVector& b     = attr_b();    

    //Init
    // r=grad - H*d
    evaluateHf(ai, xi, attr_direction, attr_r);
    r *= -1.0;
    r += gradV;    

    LongReal errorE = l2(r);

    // p=W^-1*r
    evaluateAinv(attr_r, attr_p);
    
    LongReal error  = innerProduct(r, p);
    LongReal errorO = error;    
    LongReal eps    = epsilonCR()*errorE;            
    
    LongReal lambda = 1.0, a0 = 1.0;
    
    LongInt numberOfIterations = 1;    

    while(eps<errorE && numberOfIterations <= maxSteps)
     {              
        // a = H*p
        evaluateHf(ai, xi, attr_p, attr_ap);
              
        lambda = error/innerProduct(a, p);

        x.update( lambda, p);
        r.update(-lambda, a);

        // b=W^-1*r
        evaluateAinv(attr_r, attr_b);

        error  = innerProduct(r, b);

        a0 = error/errorO;

        p *= a0;
        p.update(1.0, b);
        
        errorO = error;
                
        errorE = l2(r);
    
        numberOfIterations++;
     }

    if(errorE<=eps)
     {
        //! \todo compute beta        
        value = true;        
     }

    if(attr_printCout)
     {
        //cout << numberOfIterations-1 << " ";
     }

   return value; 
 }
