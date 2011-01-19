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

// DGKTSD.cpp: Implementierung der Klasse DGKTSD.
//
//////////////////////////////////////////////////////////////////////

#include "DGKTSD.hpp"



DGKTSD::DGKTSD()
 {
 
 }


DGKTSD::~DGKTSD()
 {

 }


DGKTSD&  DGKTSD::setParameter(const LongReal& norm, const LongReal& error, const LongInt& r)
 {
    //LongReal lambda1 = MIN(0.09*norm, 1.0e6); 
    LongReal lambda1 = 0.5;

    //LongReal lambda2 = error*1.0e-4;
    LongReal lambda2 = 0.0;

    if(r==1)
     {
        lambda2 = 0.0;
     }    

    (*this) 
     .setLambda1(lambda1)
     .setLambda2(lambda2)
    ;
    
   return (*this);
 }


DGKTSD&  DGKTSD::setParameter(const LongReal& norm, const LongInt& r)
 {
    //LongReal lambda1 = MIN(0.09*norm, 1.0e6); 
    LongReal lambda1 = 0.5;
    //LongReal lambda2 = 0.0;
    
    if(r==1)
     {
        (*this) 
         .setLambda2(0.0)
        ;
     }    

    (*this) 
     .setLambda1(lambda1)
     //.setLambda2(lambda2)     
    ;
    
   return (*this);
 }


bool DGKTSD::evaluateWHf(const DKTS& a, const DKTS& xi, const DKTS& v, DKTS& w) const
 {
    bool value = true;

    LongInt d = DGKTSD::d();
    LongInt k = DGKTSD::k();
    LongInt l = DGKTSD::l();
    LongInt m = DGKTSD::m();

    w.setNull();
    attr_w().setNull();

    for(LongInt mu1=0; mu1<d; mu1++)
     {
        for(LongInt j1=0; j1<k; j1++)
         {
            const LongReal& preFactor = W(j1, mu1);

            DKTVector& work_jmu = attr_w(j1, mu1);

            if(mu1!=0)
             {
                //The Matrix G
                const DKTVector& xjm = xi(j1, mu1);
                const DKTVector& vjm = v(j1, mu1);

                work_jmu.update(2.0*innerProduct(xjm, vjm), xjm);

                work_jmu.update((x(j1, j1, mu1) - 1.0), vjm);
                work_jmu *= 4.0*attr_lambda1*preFactor;
             }
            for(LongInt mu2=0; mu2<d; mu2++)
             {
                if(mu2!=mu1)
                 {
                    for(LongInt j2=0; j2<k; j2++)
                     {
                        //The Matrix C and lambda2*G2(C)
                        if(j1==j2)
                         {
                            work_jmu.update((1.0+2.0*attr_lambda2)*psi(j1, j2, mu1, mu2)*innerProduct(xi(j1, mu2), v(j2, mu2))*W(j2, mu2), xi(j2, mu1));

/*                            
                            // The Matrix B = x^mu1 * d_j1,mu1,mu2 * (x^mu2)^t
                            LongReal& vjm = v(j1, mu2)(0);

                            dgemv (&AMatrix::conjTrans, &m, &k, &const_eins, &xi(mu2)(0), &m, &vjm, &const_inc,
                                   &const_null, &attr_workK[0], &const_inc);

                            for(LongInt j=0; j<k; j++)
                             {
                                attr_workK[j] *= psi(j1, j, mu1, mu2);
                             }

                            dgemv (&AMatrix::notConjTrans, &m, &k, &const_eins, &xi(mu1)(0), &m, &attr_workK[0], &const_inc,
                                   &const_eins, &work_jmu(0), &const_inc);

                            // The Matrix D = a^mu1 * da_j1,mu1,mu2 * (a^mu2)^t

                            dgemv (&AMatrix::conjTrans, &m, &l, &const_minus_eins, &a(mu2)(0), &m, &vjm, &const_inc,
                                   &const_null, &attr_workL[0], &const_inc);

                            for(LongInt i=0; i<l; i++)
                             {
                                attr_workL[i] *= phi(j1, i, mu1, mu2);
                             }

                            dgemv (&AMatrix::notConjTrans, &m, &l, &const_eins, &a(mu1)(0), &m, &attr_workL[0], &const_inc,
                                   &const_eins, &work_jmu(0), &const_inc);
*/                                   
                         }
                        else
                         {
                            work_jmu.update(psi(j1, j2, mu1, mu2)*innerProduct(xi(j1, mu2), v(j2, mu2))*W(j2, mu2), xi(j2, mu1));                         
                         }
                     }
                 }                                  
             }  

            work_jmu.update(attr_lambda2*psi(j1, j1, mu1, mu1)*W(j1, mu1), v(j1, mu1));            
            work_jmu *= preFactor;
         }// End for(LongInt j1=0; j1<k; j1++)

        for(LongInt j11=0; j11<k; j11++)
         {
            DKTVector& w_jmu = w(j11, mu1);

            for(LongInt j2=0; j2<k; j2++)
             {
                w_jmu.update(A(j11, j2, mu1), attr_w(j2, mu1));
             }
         }// End for(LongInt j11=0; j11<k; j1++)
     }// End for(LongInt mu1=0; mu1<d; mu1++)
  
    w() += v();

   return value;
 }


bool DGKTSD::evaluateHf(const DKTS& a, const DKTS& xi, const DKTS& v, DKTS& w) const
 {
    bool value = true;

    LongInt d = DGKTSD::d();
    LongInt k = DGKTSD::k();

    w().setNull();

    for(LongInt mu1=0; mu1<d; mu1++)
     {
        for(LongInt j1=0; j1<k; j1++)
         {
            DKTVector& w_jmu = w(j1, mu1);
            const LongReal& preFactor = W(j1, mu1);

            if(mu1!=0)
             {
                //The Matrix G
                const DKTVector& xjm = xi(j1, mu1);
                const DKTVector& vjm = v(j1, mu1);

                w_jmu.update(2.0*innerProduct(xjm, vjm), xjm);

                w_jmu.update((x(j1, j1, mu1) - 1.0), vjm);
                w_jmu *= 4.0*attr_lambda1*preFactor;
             }
            for(LongInt mu2=0; mu2<d; mu2++)
             {
                if(mu2==mu1)
                 {
                    // The Matrix A and lambda2 G2(A)
                    for(LongInt j2=0; j2<k; j2++)
                     {
                        if(j1==j2)
                         {
                            w_jmu.update((1.0+2.0*attr_lambda2)*psi(j1, j2, mu2, mu2)*W(j2, mu2), v(j2, mu2));
                         }
                        else
                         {
                            w_jmu.update(psi(j1, j2, mu2, mu2)*W(j2, mu2), v(j2, mu2));
                         }
                     }
                 }
                else
                 {
                    //The Matrix C and lambda2*G2(C)
                    for(LongInt j2=0; j2<k; j2++)
                     {
                        if(j1==j2)
                         {
                            w_jmu.update((1.0+attr_lambda2)*psi(j1, j2, mu1, mu2)*innerProduct(xi(j1, mu2), v(j2, mu2))*W(j2, mu2), xi(j2, mu1));

/*

                    // The Matrix B = x^mu1 * d_j1,mu1,mu2 * (x^mu2)^t
                    LongReal& vjm = v(j1, mu2)(0);

                    dgemv (&AMatrix::conjTrans, &m, &k, &const_eins, &xi(mu2)(0), &m, &vjm, &const_inc,
                           &const_null, &attr_workK[0], &const_inc);

                    for(LongInt j=0; j<k; j++)
                     {
                        attr_workK[j] *= psi(j1, j, mu1, mu2);
                     }

                    dgemv (&AMatrix::notConjTrans, &m, &k, &const_eins, &xi(mu1)(0), &m, &attr_workK[0], &const_inc,
                           &const_eins, &w_jmu(0), &const_inc);

                    // The Matrix D = a^mu1 * da_j1,mu1,mu2 * (a^mu2)^t

                    dgemv (&AMatrix::conjTrans, &m, &l, &const_minus_eins, &a(mu2)(0), &m, &vjm, &const_inc,
                           &const_null, &attr_workL[0], &const_inc);

                    for(LongInt i=0; i<l; i++)
                     {
                        attr_workL[i] *= phi(j1, i, mu1, mu2);
                     }

                    dgemv (&AMatrix::notConjTrans, &m, &l, &const_eins, &a(mu1)(0), &m, &attr_workL[0], &const_inc,
                           &const_eins, &w_jmu(0), &const_inc);

*/
                            
                         }
                        else
                         {
                            w_jmu.update(psi(j1, j2, mu1, mu2)*innerProduct(xi(j1, mu2), v(j2, mu2))*W(j2, mu2), xi(j2, mu1));                         
                         }
                     }
                 }                
             }

            w_jmu *= preFactor;
         }// End for(LongInt j1=0; j1<k; j1++)
     }// End for(LongInt mu1=0; mu1<d; mu1++)    

   return value;
 }


bool DGKTSD::evaluateA(const DKTS& a, const DKTS& xi, 
                       const DKTS& v, DKTS& w,
                       const LongReal& alpha, const LongReal& beta) const
 {
    bool value = true;
    LongInt d = DGKTSD::d();
    LongInt k = DGKTSD::k();

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


bool DGKTSD::evaluateAinv(const DKTS& v,  DKTS& w) const
 {
    bool value = true;

    LongInt d = DGKTSD::d();
    LongInt k = DGKTSD::k();

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


DGKTSD& DGKTSD::solveDirection(const DKTS& a,  const DKTS& xi)
 {
    cout << setprecision(4);

    (*this)
     .useKrylovMethod(a, xi)
     //.useAlsMethod(a, xi)
    ;        

   return(*this);
 }


DGKTSD& DGKTSD::useKrylovMethod(const DKTS& a,  const DKTS& xi)
 {
    (*this)
     .invertA(0.0, 1.0)
    ;
    
    // ALS - Test
   //evaluateAinv(attr_gradient, attr_direction);
    
    //cr(a, xi);
    crA(a, xi);        
    evaluateW(attr_direction);
    
    /*
    const LongReal c = 1.0e-6;
    
    const DKTVector&  g = attr_gradient();
    DKTVector&        d = attr_direction();
    
    const LongReal l2_g = l2(g);
    
    const LongReal cosDG = innerProduct(d,g)/l2(d)/l2_g;
       
    if( cosDG < c*MIN(1,l2_g))
     {
        cout << "+";
        d = g;
     }        
    */

   return(*this);
 }


bool DGKTSD::cr(const DKTS& ai, const DKTS& xi)
 {
    bool value = true;

    DKTVector& w    = attr_w();
    DKTVector& r    = attr_r();
    DKTVector& x    = attr_direction();
    DKTVector& a    = attr_ap();
    DKTVector& aOld = attr_apo();
    DKTVector& b    = attr_b();
    DKTVector& bOld = attr_bo();
    DKTVector& p    = attr_p();
    DKTVector& pOld = attr_po();
    DKTVector& grad = attr_gradientV();

    //Init
    // r=grad - H*d
    evaluateHf(ai, xi, attr_direction, attr_r);
    r *= -1.0;
    r += grad;    

/*
    // po = 0
    attr_po.setAllEntrysNull();
*/

    // p=W^-1*r
    evaluateAinv(attr_r, attr_p);

    // a = H*p
    evaluateHf(ai, xi, attr_p, attr_ap);

    // b = W^-1*ap
    evaluateAinv(attr_ap, attr_b);
          
    // w = A*b
    evaluateHf(ai, xi, attr_b, attr_w);

    LongReal eps    = epsilonCR();
    LongReal error  = innerProduct(a, b);
    LongReal errorO = error;
    LongReal errorE = l2(r);
    
    LongReal lambda = 1.0, a0 = 1.0, a1 = 1.0;

    lambda = innerProduct(r, b)/error;
    a0     = innerProduct(w, b)/error;

    x.update( lambda, p);
    r.update(-lambda, a);

    bOld = p;  
    p *= -a0;
    p.update(1.0, b);
    pOld = bOld;

    bOld = a;
    a *= -a0;
    a.update(1.0, w);
    aOld = bOld;

    bOld = b;
    evaluateAinv(attr_ap, attr_b);
    evaluateHf(ai, xi, attr_b, attr_w);

    errorO = error;
    error  = innerProduct(a, b);

    const LongInt     maxSteps = maxStepsCR();
    LongInt numberOfIterations = 2;

    while(eps<errorE && numberOfIterations <= maxSteps)
     {
        lambda = innerProduct(r, b)/error;
        a0     = innerProduct(w, b)/error;
        a1     = innerProduct(w, bOld)/errorO;

        x.update( lambda, p);
        r.update(-lambda, a);

        bOld = p;  
        p *= -a0;
        p.update(1.0, b);
        p.update(-a1, pOld);
        pOld = bOld;


        bOld = a;
        a *= -a0;
        a.update(1.0, w);
        a.update(-a1, aOld);
        aOld = bOld;

        bOld = b;
        evaluateAinv(attr_ap, attr_b);
        evaluateHf(ai, xi, attr_b, attr_w);

        errorO = error;
        error  = innerProduct(a, b);
        errorE = l2(r);
        
        numberOfIterations++;
     }    

    cout << "l2(p) = " << l2(p) << " " << l2(r);
    
   return value;
 }


bool DGKTSD::crA(const DKTS& ai, const DKTS& xi)
 {
    bool value = true;

    DKTVector& w    = attr_bo();
    DKTVector& r    = attr_r();
    DKTVector& x    = attr_direction();
    DKTVector& a    = attr_ap();
    DKTVector& c    = attr_po();
    DKTVector& b    = attr_b();
    DKTVector& p    = attr_p();
    DKTVector& grad = attr_gradientV();

    //Init
    // r=grad - H*d
    evaluateHf(ai, xi, attr_direction, attr_r);
    r *= -1.0;
    r += grad;    

    // p=W^-1*r
    evaluateAinv(attr_r, attr_p);

    // a = H*p
    evaluateHf(ai, xi, attr_p, attr_ap);

    // b = W^-1*ap
    evaluateAinv(attr_ap, attr_b);
          
    // w = p
    w = p;
    
    LongReal error  = innerProduct(a, b);
    LongReal errorO = error;
    LongReal errorE = l2(r);
    LongReal eps    = epsilonCR();//*errorE;
    
    LongReal lambda = 1.0, a0 = 1.0;

    const LongInt     maxSteps = maxStepsCR();
    LongInt numberOfIterations = 1;    

    while(eps<errorE && numberOfIterations <= maxSteps)
     {        
        lambda = innerProduct(r, b)/error;

        x.update( lambda, p);
        r.update(-lambda, a);
        w.update(-lambda, b);        

        evaluateWHf(ai, xi, attr_b, attr_po);

        a0 = innerProduct(r, c)/error;

        p *= -a0;
        p.update(1.0, w);
        
        // a = H*p
        evaluateHf(ai, xi, attr_p, attr_ap);
        // b = W^-1*ap
        evaluateAinv(attr_ap, attr_b);

        errorO = error;
        error  = innerProduct(a, b);
        errorE = l2(r);
    
        numberOfIterations++;
     }

    if(attr_printCout)
     {
        cout << "crItr = " << setw(4) << numberOfIterations-1 << " ";
     }
/*
    if(maxSteps < numberOfIterations)
     {
        x = grad;
     }
*/
   return value;
 }
