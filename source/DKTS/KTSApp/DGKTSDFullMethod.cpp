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

// DGKTSDFullMethod.cpp: Implementierung der Klasse DGKTSDFullMethod.
//
//////////////////////////////////////////////////////////////////////

#include "DGKTSDFullMethod.hpp"
#include "LapackInterface2.hpp"
#include "BlasInterface.hpp"


char     DGKTSDFullMethod::const_u            =  'U';
char     DGKTSDFullMethod::const_jobz         =  'N';
char     DGKTSDFullMethod::const_notConjTrans =  'n';
char     DGKTSDFullMethod::const_trans        =  't';
LongInt  DGKTSDFullMethod::const_inc          =  1;
LongReal DGKTSDFullMethod::const_eins         =  1.0;
LongReal DGKTSDFullMethod::const_minus_eins   = -1.0;
LongReal DGKTSDFullMethod::const_null         =  0.0;


DGKTSDFullMethod::DGKTSDFullMethod()
 {
 }


DGKTSDFullMethod::~DGKTSDFullMethod()
 {

 }


DKTSDIterationInfo DGKTSDFullMethod::startIteration(const DKTS& a, DKTS& xi, const LongReal& normA, DKTSDDataBlock& dataBlock)
 {
    DKTSDIterationInfo infoEntry;

    const LongInt maxSteps      = DGKTSDFullMethod::maxSteps(); 
    const LongInt maxStepsAmijo = DGKTSDFullMethod::maxStepsAmijo();

    LongInt numberOfIterations = 1, numberOfIterationsAmijo = 1;

    DKTVector& grad = attr_gradient();
    DKTVector& d    = attr_direction();
    DKTVector& x    = xi(); 

    const LongReal eps   = epsilon();
    const LongReal sigma = DGKTSDFullMethod::sigma();
    const LongReal prec  = DGKTSDFullMethod::precision();
    
    const LongReal invNormA = 1.0/normA;

    const LongInt order    = a.d();
    const LongInt kRankR   = a.k();
    const LongInt kRankr   = xi.k();
    const LongInt spanSize = a.n();
   
    bool errorFlag = true;    
   
    IString caption = (IString("d ") + IString("$ = ") + IString(order)    + IString(",$") + IString("\\hspace{15pt}")
                     + IString("R ") + IString("$ = ") + IString(kRankR)   + IString(",$") + IString("\\hspace{15pt}")
	                    + IString("r ") + IString("$ = ") + IString(kRankr)   + IString(",$") + IString("\\hspace{15pt}")
	                    + IString("m ") + IString("$ = ") + IString(spanSize) + IString("$"));
    
    dataBlock.setCaptionString(caption);
    
    if(a.k()==xi.k())
     {
        xi = a;
        errorFlag = false;
     }

    xi.balanced();
    (*this)     
     .initialValuesOfax(a, xi)
    ;        
            
    LongReal dist    = sqrt(fabs(normA*normA + 2.0*f()));
    LongReal errOld  = dist*invNormA;
    LongReal errNew  = errOld;            
        
    (*this)
     .setLambda1(attr_error*invNormA)
    ;

    if(xi.k()==1)
     {
        (*this)
         .setLambda2(0.0)
        ;          
     }
    else
     {
        (*this)
         .setLambda2(errOld*1.0e-4)
        ;     
     }

    (*this)
     .solveGradient(a, xi)
     .evaluateW(attr_gradient, attr_gradientV);
    ;

    LongReal delta   = l2(x);
    LongReal alpha   = 2.0;        
    LongReal dg      = 1.0;
    LongReal fOld    = F();
    LongReal fNew    = 0.5*fOld;
    LongReal decr    = 1.0;
    
    attr_error   = l2(grad);

    if(attr_error<1.0e-3)
     {     
        (*this)
         .setMaxStepsCRL1(30)
         .setMaxStepsCRL2(5)
         .setGamma(0.5)
        ;
     }
    else if(attr_error<1.0e-4)
     {
        (*this)
         .setMaxStepsCRL1(20)
         .setMaxStepsCRL2(8)
         .setGamma(0.01)
        ;     
     }

    if(attr_printCout)
     {
        cout << setprecision(2) << scientific;
        cout << "Working Memory      = " << setw(4) << (LongReal)(memory()*sizeof(LongReal))/(LongReal)1048576 << " MByte" << endl;
								cout << endl;
        cout << setprecision(6);
        cout << "|A-X_0|/|A|         = " << dist/normA << '\t' << "|f'(x_0)| = " << attr_error << endl; //<< '\t' << "|A| = " << normA << endl;    
        cout << endl;    
        cout << setprecision(6)          << setiosflags(ios::left);     
     }    

    DKTSDData data0(0, errNew, attr_error, delta, 1.0, fabs(0.5*normA*normA + fOld)/normA);

    dataBlock.addEntry(data0);

    while((eps<attr_error || prec <= decr) && numberOfIterations <= maxSteps && errorFlag)
     {               
        LongReal temp = 1;
	
	       if(attr_printCout)
         {
            cout << "itr. = " << setw(5) << numberOfIterations;
         }
        
        solveDirection(a, xi);

        dg = innerProduct(d, grad);

        if(dg < 0.0)
         {
            temp = -1;
            dg *= -1.0;
            d  *= -1.0;
         }
        
        prepareValuesOfadx(a, xi);
        fNew = F(alpha);         
         
        numberOfIterationsAmijo = 0;
        alpha = 2.0;

         do
          {
             alpha *= 0.5;          
             fNew   = F(alpha);          
             numberOfIterationsAmijo++;
             
          } while(fOld-fNew < sigma*alpha*dg && numberOfIterationsAmijo < maxStepsAmijo);         

         if(fOld < fNew)
          { 
             errorFlag = false;
          }
         else
          {
             x.update(-alpha, d);
          }
         
         (*this)
          .initialAllValues(a, xi)
         ;         

         attr_error = l2(grad);
         delta      = l2(d);

         if((fabs(alpha-1.0)<EPS_NULL || decr<prec) && attr_useNewtonMethode)
          {
             setKappa(1.0);
             setLambda1(0.0);
             
             if(attr_printCout)
              {
                 cout << " n ";
              }
          }
         else
          {
             setKappa(0.0);
             setLambda1(attr_error*invNormA);
             
             if(attr_printCout)
              {
                 cout << " l ";
              }
          }                  
         
         if(errorFlag==true)
          {
             fOld   = fNew;
             errNew = sqrt(fabs(normA*normA + 2.0*f()))*invNormA;
             decr   = (errOld - errNew)/errOld;         
             errOld = errNew;
          }

         DKTSDData data(temp*numberOfIterations, errNew, attr_error, delta, alpha, fabs(0.5*normA*normA + fNew)/normA);
        
         if(attr_printCout)
          {
             cout <<         "|f'(x_k)| = "  << scientific << setw(11) << attr_error; 
             cout << '\t' << "|delta| = "    << setw(10)   << delta;
             cout << resetiosflags( ::std::ios::scientific );
													cout << '\t' << "a_k = "            << setw(6)    << alpha;
													cout << '\t' << "|A-X_itr|/|A| = "  << scientific << setw(10) << errNew;
													cout << endl;
												
             cout << resetiosflags(::std::ios::scientific);												 
          }
         
         numberOfIterations++;
	 
	        dataBlock.addEntry(data);
     }    
    
    numberOfIterations--;
 
    if(!errorFlag)
     {
        cout << "Warning, f(x_(k-1)) < f(x_k) !!!" << endl;
     }	

    const LongReal dOld = dist*invNormA;
    const LongReal dNew = errNew;
    const LongReal diff = dOld - dNew;

    if(attr_printCout)
     {
        cout << setprecision(attr_preC);
        cout << endl;
        cout << "----------------------------------" << endl;
        cout << "oldAppError = " << dOld << endl;
        cout << "newAppError = " << dNew << endl;
        cout << "----------------------------------" << endl;
        cout << "diff        = " << diff << endl;
        cout << setprecision(4);
        cout << "diff[%]     = " << diff/dOld*100 << "%" << endl;
        cout << "----------------------------------" << endl;
        cout << endl;
     }
     
    const LongReal quality = 1.0e-2;

    if(attr_error<=eps)
     {
        if(delta <=quality)
         {
            (*this)
             .setQuality(1)
            ;
         }
        else
         {
            (*this)
             .setQuality(2)
            ;

         }
     }
    else
     {
        if(delta <=quality && attr_error <= 1.0e-1)
         {            
            (*this)
             .setQuality(3)
            ;
         }
        else 
         {
            (*this)
             .setQuality(5)
            ;

         }
     }         

    //Werte werden ins Datenfeld geschrieben
    infoEntry.setStartError(dOld);
    infoEntry.setNumberOfNewtonSteps(numberOfIterations);
    infoEntry.setRelativeDifferenz(diff/dOld*100);
    infoEntry.setError(dNew);
				infoEntry.setGradient(attr_error);
    
   return infoEntry;
 }


DGKTSDFullMethod& DGKTSDFullMethod::initialValuesOfa(const DKTS& a, const DKTS& xi)
 {
    const LongInt d  = DGKTSDFullMethod::d();
    const LongInt k  = DGKTSDFullMethod::k();
    const LongInt l  = DGKTSDFullMethod::l();

    LongInt index = 0;

    for(LongInt j=0; j<k; j++)
     {
        for(LongInt i=0; i<l; i++)
         {
            index = isIndexOffsetOfa(j, i);

            for(LongInt mu=0; mu<d; mu++)
             {             
                attr_a[index+mu] = innerProduct(a(i, mu), xi(j, mu));
             } 

            LongInt index_ji = (j*l+i);
            LongInt index2 = 0;

            for(LongInt mu1=0; mu1<d; mu1++)
             {
                for(LongInt mu2=0; mu2<=mu1; mu2++)
                 {
                    index2 = isIndexOffsetOfmuA(mu1, mu2);
                    LongReal& wert = attr_phi[index_ji+index2];

                    wert = 1.0;

                    for(LongInt nu=0; nu<d; nu++)
                     {
                        if(nu!=mu1 && nu!=mu2)
                         {
                            wert *= attr_a[index+nu];
                         }                         
                     }
                 }//End : for(LongInt mu2=0; mu2<=mu1; mu2++)
             }//End : for(LongInt mu1=0; mu1<d; mu1++)
         }
     }

   return (*this);
 }


DGKTSDFullMethod& DGKTSDFullMethod::initialValuesOfax(const LongReal& alpha)
 {
    const LongInt d  = DGKTSDFullMethod::d();
    const LongInt k  = DGKTSDFullMethod::k();
    const LongInt l  = DGKTSDFullMethod::l();

    LongInt index = 0;

    LongReal wert1 = 0.0;

    for(LongInt j1=0; j1<k; j1++)
     {
        for(LongInt i=0; i<l; i++)
         {
            index = isIndexOffsetOfa(j1, i);

            for(LongInt mu=0; mu<d; mu++)
             {             
                attr_a[index+mu] += -alpha*attr_ad[index+mu];
             } 

            LongInt index_ji = (j1*l+i);
            LongInt index2 = 0;

            for(LongInt mu1=0; mu1<d; mu1++)
             {
                for(LongInt mu2=0; mu2<=mu1; mu2++)
                 {
                    index2 = isIndexOffsetOfmuA(mu1, mu2);
                    LongReal& wert = attr_phi[index_ji+index2];

                    wert = 1.0;

                    for(LongInt nu=0; nu<d; nu++)
                     {
                        if(nu!=mu1 && nu!=mu2)
                         {
                            wert *= attr_a[index+nu];
                         }                         
                     }
                 }//End : for(LongInt mu2=0; mu2<=mu1; mu2++)
             }//End : for(LongInt mu1=0; mu1<d; mu1++)
         }// End : for(LongInt i=0; i<l; i++)

        for(LongInt j2=0; j2<=j1; j2++)
         {
            index = isIndexOffsetOfx(j1, j2);

            for(LongInt mu=0; mu<d; mu++)
             {
                attr_x[index+mu] += (- alpha*(xd(j1, j2, mu) + xd(j2, j1, mu)) + alpha*alpha*dd(j1, j2, mu));
             }

            LongInt index_jj1 = j1*k + j2;
            LongInt index_jj2 = j2*k + j1;

            LongInt index2 = 0;

            for(LongInt mu1=0; mu1<d; mu1++)
             {                
                for(LongInt mu2=0; mu2<=mu1; mu2++)
                 {
                    index2 = isIndexOffsetOfmuX(mu1, mu2);

                    wert1 = 1.0;

                    for(LongInt nu=0; nu<d; nu++)
                     {
                        if(nu!=mu1 && nu!=mu2)
                         {
                            wert1 *= attr_x[index+nu];
                         }                         
                     }

                    attr_psi[index_jj1+index2] = wert1;                       
                    attr_psi[index_jj2+index2] = wert1;

                 }//End : for(LongInt mu2=0; mu2<=mu1; mu2++)

                const LongReal wert = 1.0/sqrt(wert1);

                if(j1==j2)
                 {                    
                    W(j1, mu1) = wert;
                 }

                A(j1, j2, mu1) = wert1; // (wert1=psi(j1, j2, mu1, mu1);)

             }//End : for(LongInt mu1=0; mu1<d; mu1++)
         }

        for(LongInt mu1=0; mu1<d; mu1++)
         {
            const LongReal& wj1m1 = W(j1, mu1);
           
            for(LongInt j2=0; j2<=j1; j2++)
             {
                A(j1, j2, mu1) *= (wj1m1*W(j2, mu1));
             }
         }
     }    

   return (*this);
 }


DGKTSDFullMethod& DGKTSDFullMethod::initialValuesOfax(const DKTS& a, const DKTS& xi)
 {
    const LongInt d  = DGKTSDFullMethod::d();
    const LongInt k  = DGKTSDFullMethod::k();
    const LongInt l  = DGKTSDFullMethod::l();

    LongInt index = 0;

    LongReal wert1 = 0.0;

    for(LongInt j1=0; j1<k; j1++)
     {
        for(LongInt i=0; i<l; i++)
         {
            index = isIndexOffsetOfa(j1, i);

            for(LongInt mu=0; mu<d; mu++)
             {             
                attr_a[index+mu] = innerProduct(a(i, mu), xi(j1, mu));
             } 

            LongInt index_ji = (j1*l+i);
            LongInt index2 = 0;

            for(LongInt mu1=0; mu1<d; mu1++)
             {
                for(LongInt mu2=0; mu2<=mu1; mu2++)
                 {
                    index2 = isIndexOffsetOfmuA(mu1, mu2);
                    LongReal& wert = attr_phi[index_ji+index2];

                    wert = 1.0;

                    for(LongInt nu=0; nu<d; nu++)
                     {
                        if(nu!=mu1 && nu!=mu2)
                         {
                            wert *= attr_a[index+nu];
                         }                         
                     }
                 }//End : for(LongInt mu2=0; mu2<=mu1; mu2++)
             }//End : for(LongInt mu1=0; mu1<d; mu1++)
         }// End : for(LongInt i=0; i<l; i++)

        for(LongInt j2=0; j2<=j1; j2++)
         {
            index = isIndexOffsetOfx(j1, j2);

            for(LongInt mu=0; mu<d; mu++)
             {
                attr_x[index+mu] = innerProduct(xi(j1, mu), xi(j2, mu));
             }

            LongInt index_jj1 = j1*k + j2;
            LongInt index_jj2 = j2*k + j1;

            LongInt index2 = 0;

            for(LongInt mu1=0; mu1<d; mu1++)
             {                
                for(LongInt mu2=0; mu2<=mu1; mu2++)
                 {
                    index2 = isIndexOffsetOfmuX(mu1, mu2);

                    wert1 = 1.0;

                    for(LongInt nu=0; nu<d; nu++)
                     {
                        if(nu!=mu1 && nu!=mu2)
                         {
                            wert1 *= attr_x[index+nu];
                         }                         
                     }

                    attr_psi[index_jj1+index2] = wert1;                       
                    attr_psi[index_jj2+index2] = wert1;

                 }//End : for(LongInt mu2=0; mu2<=mu1; mu2++)

                const LongReal wert = 1.0/sqrt(wert1);

                if(j1==j2)
                 {                    
                    W(j1, mu1) = wert;
                 }
                
                A(j1, j2, mu1) = wert1; // (wert1=psi(j1, j2, mu1, mu1);)                

             }//End : for(LongInt mu1=0; mu1<d; mu1++)
         }// End : for(LongInt j2=0; j2<=j1; j2++)

        for(LongInt mu1=0; mu1<d; mu1++)
         {
            const LongReal& wj1m1 = W(j1, mu1);
                      
            for(LongInt j2=0; j2<=j1; j2++)
             {
                A(j1, j2, mu1) *= (wj1m1*W(j2, mu1));
             }
         }
     }// End : for(LongInt j1=0; j1<k; j1++)

   return (*this);
 }


DGKTSDFullMethod& DGKTSDFullMethod::prepareValuesOfadx(const DKTS& a, const DKTS& xi)
 {
    const LongInt d  = DGKTSDFullMethod::d();
    const LongInt k  = DGKTSDFullMethod::k();
    const LongInt l  = DGKTSDFullMethod::l();

    LongInt index1=0, index2=0, index3=0;

    LongReal wert1 = 0.0;

    for(LongInt j1=0; j1<k; j1++)
     {
        LongInt j2 = 0;

        for(j2=0; j2<=j1; j2++)
         {
            index1 = isIndexOffsetOfx(j1, j2);
            index2 = isIndexOffsetOfxd(j1, j2);

            for(LongInt mu=0; mu<d; mu++)
             {
                attr_dd[index1+mu] = innerProduct(attr_direction(j1, mu), attr_direction(j2, mu));
                attr_xd[index2+mu] = innerProduct(xi(j1, mu), attr_direction(j2, mu));
             }
         }

        for(j2=j1+1; j2<k; j2++)
         {
            index2 = isIndexOffsetOfxd(j1, j2);

            for(LongInt mu=0; mu<d; mu++)
             {
                attr_xd[index2+mu] = innerProduct(xi(j1, mu), attr_direction(j2, mu));
             }
         }

        for(LongInt i=0; i<l; i++)
         {
            index3 = isIndexOffsetOfa(j1, i);

            for(LongInt mu=0; mu<d; mu++)
             {             
                attr_ad[index3+mu] = innerProduct(a(i, mu), attr_direction(j1, mu));
             } 
         }
     }
    
 
   return (*this);
 }


DGKTSDFullMethod& DGKTSDFullMethod::invertA(const LongReal& alpha, const LongReal& beta)
 {
    LongInt d    = DGKTSDFullMethod::d();
    LongInt k    = DGKTSDFullMethod::k();
    LongInt inc  = 1;
    LongInt n    = d*k*k;

    LongInt info = 0;

    TensorCalculus::Blas<LongReal>::scal(n, beta, &A(0), inc);

    if(fabs(alpha) < EPS_NULL)
     {
        for(LongInt mu=0; mu<d; mu++)
         {
        	info = TensorCalculus::Lapack<LongReal>::potrf(const_u, k, &A(mu), k);
        	info = TensorCalculus::Lapack<LongReal>::potri(const_u, k, &A(mu), k);
         }
     }
    else
     {
        for(LongInt mu=0; mu<d; mu++)
         {
            for(LongInt j1=0; j1<k; j1++)
             {
                A(j1, j1, mu) += alpha;
             }

            info = TensorCalculus::Lapack<LongReal>::potrf(const_u, k, &A(mu), k);
            info = TensorCalculus::Lapack<LongReal>::potri(const_u, k, &A(mu), k);
         }
     }

   return (*this);
 }


LongReal DGKTSDFullMethod::optimalOmega()
 {
    LongReal value = 0.0;

    LongInt d     = DGKTSDFullMethod::d();
    LongInt k     = DGKTSDFullMethod::k();
    LongInt info  = 0;
    LongInt dim   = k*k;
    LongInt inc   = 1;
    LongInt lwork = MAX(1, 3*k-1);

    LongReal lmin =  upperBound();
    LongReal lmax = -lmin;

    for(LongInt mu=0; mu<d; mu++)
     {
        TensorCalculus::Blas<LongReal>::copy(dim, &A(mu), inc, &attr_eigA[0], inc);
        info = TensorCalculus::Blas<LongReal>::syev(const_jobz, const_u, k, &attr_eigA[0], k, &attr_eigValue[0], &attr_eigWork[0], lwork);
 
        lmin = MIN(lmin, attr_eigValue[0]);
        lmax = MAX(lmax, attr_eigValue[k-1]);
     }

    value = sqrt(fabs(lmax*lmin));

   return value;
 }


LongReal DGKTSDFullMethod::f()const
 {
    const LongInt d = DGKTSDFullMethod::d();
    const LongInt k = DGKTSDFullMethod::k();
    const LongInt l = DGKTSDFullMethod::l();

    LongReal value=0.0;

    for(LongInt j1=0; j1<k; j1++)
     {
        for(LongInt i=0; i<l; i++)
         {
            value -= phi(j1, i, 0, 0)*a(j1, i, 0);            
         }

        value += 0.5*psi(j1, j1, 0, 0)*x(j1, j1, 0);

        for(LongInt j2=0; j2<j1; j2++)
         {
            value += psi(j1, j2, 0, 0)*x(j1, j2, 0);
         }
     }

   return value;
 }


LongReal DGKTSDFullMethod::F()const
 {
    const LongInt d = DGKTSDFullMethod::d();
    const LongInt k = DGKTSDFullMethod::k();
    const LongInt l = DGKTSDFullMethod::l();

    LongReal f1=0.0, g2=0.0, value=0.0;

    for(LongInt j1=0; j1<k; j1++)
     {
        // f1 
        for(LongInt i=0; i<l; i++)
         {
            f1 -= phi(j1, i, 0, 0)*a(j1, i, 0);            
         }

        f1 += 0.5*psi(j1, j1, 0, 0)*x(j1, j1, 0);

        for(LongInt j2=0; j2<j1; j2++)
         {
            f1 += psi(j1, j2, 0, 0)*x(j1, j2, 0);
         }                  

        // g2(x)
        g2 += x(j1, j1, 0)*psi(j1, j1, 0, 0);
        
     }// End for(LongInt j1=0; j1<k; j1++)
  
    value = attr_kappa1*f1 + 0.5*attr_lambda2*attr_kappa3*g2;

   return value;
 }


LongReal DGKTSDFullMethod::g2()const
 {
    const LongInt k = DGKTSDFullMethod::k();

    LongReal value=0.0, g2=0.0;

    for(LongInt j1=0; j1<k; j1++)
     {
        g2 += x(j1, j1, 0)*psi(j1, j1, 0, 0);        
     }
 
    value = 0.5*attr_kappa3*g2;

   return value;
 }


LongReal DGKTSDFullMethod::F(const LongReal& alpha)const
 {
    const LongInt d = DGKTSDFullMethod::d();
    const LongInt k = DGKTSDFullMethod::k();
    const LongInt l = DGKTSDFullMethod::l();

    LongReal value=0.0, val=1.0, f1=0.0, temp=1.0, g2=0.0;

    for(LongInt j1=0; j1<k; j1++)
     {
        for(LongInt i=0; i<l; i++)
         {
            val = 1.0;
            for(LongInt mu=0; mu<d; mu++)
             {
                val *= (a(j1, i, mu) - alpha*ad(j1, i, mu));
             }
            f1 -= val;
         }
        
        LongInt mu = 0;
        
        val   = 1.0;                
        for(mu=0; mu<d; mu++)
         {
            val   *= (x(j1, j1, mu) - alpha*2.0*xd(j1, j1, mu) + alpha*alpha*dd(j1, j1, mu));
         }
        f1 += 0.5*val;
        
        g2    += val;

        for(LongInt j2=0; j2<j1; j2++)
         {
            val = 1.0;
            for(LongInt mu=0; mu<d; mu++)
             {
                val *= (x(j1, j2, mu) - alpha*(xd(j1, j2, mu) + xd(j2, j1, mu)) + alpha*alpha*dd(j1, j2, mu));
             }
            f1 += val;
         }
                          
     }// End for(LongInt j1=0; j1<k; j1++)

    value = attr_kappa1*f1 + 0.5*attr_lambda2*attr_kappa3*g2;

   return value;
 }


DGKTSDFullMethod& DGKTSDFullMethod::solveGradient(const DKTS& a, const DKTS& xi)
 {
    LongInt d = DGKTSDFullMethod::d();
    LongInt k = DGKTSDFullMethod::k();
    LongInt l = DGKTSDFullMethod::l();
    LongInt m = DGKTSDFullMethod::m();

    LongReal fac1a = -1.0*attr_kappa1;
    LongReal fac1b =  attr_kappa1;
    LongReal fac2  =  attr_lambda2*attr_kappa3;    

    for(LongInt mu1=0; mu1<d; mu1++)
     {
        for(LongInt j1=0; j1<k; j1++)
         {
            DKTVector& gradV = attr_gradient(j1, mu1);
            LongReal& grad = gradV(0);
            
            TensorCalculus::Blas<LongReal>::gemv(DGKTSDFullMethod::const_notConjTrans, m, l, fac1a, &a(mu1)(0), m, &phi(j1, 0, mu1, mu1), const_inc,
                    const_null, &grad, const_inc);


            TensorCalculus::Blas<LongReal>::gemv(DGKTSDFullMethod::const_notConjTrans, m, k, fac1b, &xi(mu1)(0), m, &psi(j1, 0, mu1, mu1), const_inc,
                    const_eins, &grad, const_inc);
                               
            // grad(g2(x))
            gradV.update(fac2*psi(j1, j1, mu1, mu1), xi(j1, mu1));
         }
     }

   return (*this);
 }


DGKTSDFullMethod& DGKTSDFullMethod::initialAllValues(const DKTS& a, DKTS& xi)
 {  
    xi.balanced();
    (*this)     
     .initialValuesOfax(a, xi)
     //.checkSummands()
     .solveGradient(a, xi)
     .evaluateW(attr_gradient, attr_gradientV);
    ;

   return (*this);
 }


DGKTSDFullMethod& DGKTSDFullMethod::initialAllValues(const DKTS& a, DKTS& xi, const LongReal& alpha)
 {
    xi.balanced();
    (*this)
     .initialValuesOfax(alpha)
     //.checkSummands()
     .solveGradient(a, xi)
    ;
   return (*this);
 }

bool DGKTSDFullMethod::evaluateW(const DKTS& v, DKTS& w) const
 {
    bool value = true;

    LongInt d = DGKTSDFullMethod::d();
    LongInt k = DGKTSDFullMethod::k();

    w() = v();

    for(LongInt mu1=0; mu1<d; mu1++)
     {
        for(LongInt j1=0; j1<k; j1++)
         {
            w(j1, mu1) *= W(j1, mu1);

         }// End for(LongInt j1=0; j1<k; j1++)
     }// End for(LongInt mu1=0; mu1<d; mu1++)

   return value;
 }


bool DGKTSDFullMethod::evaluateW(DKTS& w) const
 {
    bool value = true;

    LongInt d = DGKTSDFullMethod::d();
    LongInt k = DGKTSDFullMethod::k();

    for(LongInt mu1=0; mu1<d; mu1++)
     {
        for(LongInt j1=0; j1<k; j1++)
         {
            w(j1, mu1) *= W(j1, mu1);

         }// End for(LongInt j1=0; j1<k; j1++)
     }// End for(LongInt mu1=0; mu1<d; mu1++)

   return value;
 }
