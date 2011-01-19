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

// DKTS.cpp: Implementierung der Klasse DKTS.
//
//////////////////////////////////////////////////////////////////////

#include "DKTS.hpp"
#include "LapackInterface2.hpp"
#include "BlasInterface.hpp"

char     DKTS::const_u            =  'U';
char     DKTS::const_jobz         =  'N';
char     DKTS::const_notConjTrans =  'n';
char     DKTS::const_trans        =  't';
LongInt  DKTS::const_inc          =  1;
LongReal DKTS::const_eins         =  1.0;
LongReal DKTS::const_minus_eins   = -1.0;
LongReal DKTS::const_null         =  0.0;




DKTS::DKTS(const LongInt& d, const LongInt& k, const LongInt& n)
:KTS(d, k, n)
 {    
    (*this)
     .allocateDataSpace(d, k, n)
    ;
 }

    
DKTS::DKTS(const DKTS& A)
:KTS(A.d(), A.k(), A.n())
 {
    (*this)
     .allocateDataSpace(A.d(), A.k(), A.n())
    ;

    (*this) = A;
 }


DKTS::~DKTS()
 {
    (*this)
     .deleteDataSpace()
    ;
 }


DKTS& DKTS::resize(const LongInt& d, const LongInt& k, const LongInt& n)
 {    
    const LongInt dOld = DKTS::d();
    const LongInt kOld = DKTS::k();
    const LongInt nOld = DKTS::n();
                
    if(dOld!=d || kOld!=k || nOld!=n)
     {
        (*this)
         .deleteDataSpace()
         .allocateDataSpace(d, k, n)
        ;
     }
    else
     {
        (*this)
         .setNull()
        ;
     }

   return (*this);
 }
 
 
DKTS& DKTS::allocateDataSpace(const LongInt& d, const LongInt& k, const LongInt& n)
 {
    (*this)
     .setDimension(d*k*n)
     .setKN(k*n)
     .setD(d)
     .setK(k)
     .setN(n)     
    ; 

    const LongInt dim_values = dimension();
    const LongInt dim_kn     = kn();

/*    
     try
      {             
*/      
         attr_values                    = (LongRealPointer)  new LongReal[dim_values];    
         attr_vecOfTensor               = (DKTVectorPointer) new DKTVector[1];
         attr_vecOfDimension            = (DKTVectorPointer) new DKTVector[d];
         attr_vecOfRepresentationVector = (DKTVectorPointer) new DKTVector[d*k];
/*         
      }
     catch (bad_alloc & exc)
      {
        cout << "Warning in DKTS::allocateDataSpace(const LongInt& d, const LongInt& k, const LongInt& n) !!!" << endl;            
        exit(1);
      }
*/
    //set
    (*this)().set(attr_dimension, attr_values[0]);
    
    for(LongInt mu=0; mu<d; mu++)
     {
        (*this)(mu).set(attr_kn, (*this)(0, mu, 0));
        
        for(LongInt j=0; j<k; j++)
         {
            (*this)(j, mu).set(attr_n, (*this)(j, mu, 0));
         }
     }

    // set all entrys
    (*this).setNull();
    
   return (*this);
 }
 
 
DKTS& DKTS::deleteDataSpace()
 {
    (*this)
     .setDimension(0)
     .setKN(0)
     .setD(0)
     .setK(0)
     .setN(0)
    ; 
        
    delete [] attr_values;
    delete [] attr_vecOfTensor;    
    delete [] attr_vecOfDimension;
    delete [] attr_vecOfRepresentationVector;
   
   return (*this);
 }


DKTS& DKTS::operator = (const DKTS& A)
 {
    setCopyOf(A);
      
   return (*this);
 }


DKTS& DKTS::operator *= (const LongReal& alpha)
 {
    const LongInt d = DKTS::d();
    const LongInt k = DKTS::k();

    if(fabs(alpha)<EPS_NULL)
     {
        (*this)().setAllEntrys(0.0);
     }
    else
     {
        LongReal sign = 1.0;
        
        if(alpha<0.0)
         {
            sign = -1.0; 
         }        
        
        const LongReal factor = pow(fabs(alpha), 1.0/((LongReal)d));
    
        for(LongInt mu=0; mu<d; mu++)
         {
            for(LongInt j=0; j<k; j++)
             {
                if(mu==0)
                 {
                    (*this)(j, mu) *= sign*factor;
                 }
                else
                 {
                    (*this)(j, mu) *= factor;
                 }
             }
         }
     }
     
   return (*this);
 }


ostream& operator << (ostream& s, const DKTS& A)
 {
    const LongInt d = A.d();
    const LongInt k = A.k();
    const LongInt n = A.n();    
        
    s << "(" << d << ", " << k << ", " << n << ")"<< endl;
    
    s << setprecision(3) << scientific << setiosflags(ios::left);
    
    for(LongInt i=0; i<n; i++)
     {    
//         for(LongInt j=0; j<k; j++)    
        for(LongInt mu=0; mu<d; mu++)        
         {
            //for(LongInt mu=0; mu<d; mu++)        
            for(LongInt j=0; j<k; j++)
             {
                s << setw(8);
                s << A(j, mu, i) << " ";
             }  
            s << " | ";          
         }
        s << endl; 
     }     
    
   return s;
 }
 
 DKTS& DKTS::setCopyOf(const DKTS& A)
 {
    const LongInt d = A.d();
    const LongInt k = A.k();
    const LongInt n = A.n();

    (*this)
     .resize(d, k, n)     
    ; 
    
    ((*this)()) = (A());

   return (*this);
 }


DKTS& DKTS::isPartialCopyOf(const DKTS& A)
 {
    const LongInt d  = A.d();
    const LongInt k  = DKTS::k();
    const LongInt n  = A.n();
    const LongInt k1 = MIN(A.k(), k);

    (*this)
     .resize(d, k, n)
     .copyFromTo(A, 0, 0, k1-1)
    ;    

   return (*this);
 }


DKTS& DKTS::copyFromTo(const DKTS& A, const LongInt& os, const LongInt& j0, const LongInt& j1)
 {
    const LongInt d = DKTS::d();
    
    for(LongInt mu=0; mu<d; mu++)
     {
        for(LongInt j=j0; j<=j1; j++)
         {
            (*this)(os+j, mu) = A(j, mu);
         }
     }

   return (*this);
 }


DKTS& DKTS::setSummandOf(const DKTS& A, const LongInt& index)
 {
    const LongInt d = A.d();
    const LongInt n = A.n();
    
    (*this)
     .resize(d, 1, n)
    ;
        
    for(LongInt mu=0; mu<d; mu++)
     {
        (*this)(0, mu) = A(index, mu);
     }        
 
   return (*this);
 }


LongReal DKTS::setMainPartOf(const DKTS& A, const LongReal& eps)
 { 
    const LongInt d = A.d();
    const LongInt r = A.k();
    const LongInt n = A.n();

    LongReal normA = frobeniusNorm(A);
 
    LongReal normR = 0.0;
    
    LongInt newR = r;
    
    for(newR=r; (1<newR) && (normR<eps*normA); newR--)
     {
        const LongReal temp = A.frobeniusNormOfSummand(newR-1);
        
        normR += temp;
     }
    
    newR++;
    
    newR = MIN(r, newR);
    
    (*this)
     .resize(d, newR, n)
    ;                
        
    for(LongInt i=0; i<newR; i++)
     {
        for(LongInt mu=0; mu<d; mu++)
         {
            (*this)(i,mu) = A(i,mu);
         }
     }
 
   return normA;
 }

DKTS& DKTS::setHadamardProductOf(const LongReal& a, const DKTS& A, const DKTS& B)
 { 
    const LongInt dA = A.d(); 
    const LongInt rA = A.k(); 
    const LongInt nA = A.n(); 
    
    const LongInt dB = B.d(); 
    const LongInt rB = B.k(); 
    const LongInt nB = B.n(); 
    
    const LongInt d  = MIN(dA, dB);
    const LongInt r  = rA*rB;
    const LongInt n  = MIN(nA, nB);
    
    (*this)
     .resize(d, r, n)
    ;
        
    for(LongInt i1=0; i1<rA; i1++)
     {
        for(LongInt i2=0; i2<rB; i2++)
         {
            const LongInt index = i1*rB + i2;
            
            for(LongInt mu=0; mu<d; mu++)
             {
                (*this)(index,mu).pwProduct(1.0, A(i1,mu), B(i2,mu));
             }            
         }
     }
    
    (*this) *= a;
    
   return (*this);
 }


DKTS& DKTS::successiveApproximation(const DKTS& A, const LongReal& eps, const LongInt& maxSteps)
 {
    const LongInt d = A.d();
    const LongInt n = A.n();
    
    const LongReal eps2 = eps*frobeniusNorm(A);
    
    (*this)
     .setAlsApproximationOf(A, maxSteps)
    ;            
    
    DKTS xi, temp, rho;
    
    rho.setSumOf(1.0, A, -1.0, (*this));
    
    LongReal error = frobeniusNorm(rho);
    
    while(error > eps2)
     {
        xi.setAlsApproximationOf(rho, maxSteps);

        temp = (*this);
        (*this).setSumOf(1.0, temp, 1.0, xi);        
        (*this).improveApproximation(A, 10, 2);
        
        temp = rho;
        rho.setSumOf(1.0, temp, -1.0, xi);
        
        error = frobeniusNorm(rho);        
     }
 
   return (*this);
 }


DKTS& DKTS::sortIndexSet()
 {
    const LongInt d = DKTS::d();
    const LongInt R = DKTS::k();
               
    LongIntPointer  indexSet = (LongIntPointer)  new LongInt[R];
    LongRealPointer values   = (LongRealPointer) new LongReal[R];

    for(LongInt i=0; i<R; i++)
     {
        values[i]   = frobeniusNormOfSummand(i);        
        indexSet[i] = i;
     }
    
    quickSort(0, R-1, values, indexSet);
    
    DKTS A((*this));

    for(LongInt i=0; i<R; i++)
     {
        for(LongInt mu=0; mu<d; mu++)
         {
            (*this)(i,mu) = A(indexSet[i], mu);
         }
     }

    delete [] indexSet;
    delete [] values;

   return (*this); 
 }

 
LongInt DKTS::partition(LongInt low, LongInt high, LongRealPointer f, LongIntPointer indexSet)
 {
    LongInt i,j;
    LongReal pivot = f[indexSet[low]];
    i = low;

    for(j=low+1; j<=high; j++) 
     {
       if (pivot<=f[indexSet[j]])
        {
           i++;
           swapIndex(i, j, f, indexSet);
        }
     }

    swapIndex(i, low, f, indexSet);

   return i;
 }

 
 
DKTS& DKTS::quickSort(LongInt low, LongInt high, LongRealPointer f, LongIntPointer indexSet)
 {
    LongInt m = 0;

    if (low < high) 
     {
        m = partition(low, high, f, indexSet);
        quickSort(low, m-1, f, indexSet);
        quickSort(m+1, high, f, indexSet);
     }

   return (*this); 
 }
 
 
DKTS& DKTS::swapIndex(const LongInt& i, const LongInt& j, LongRealPointer f, LongIntPointer indexSet)
 {
    LongInt& a_i = indexSet[i];
    LongInt& a_j = indexSet[j];

    const LongInt temp = a_i;

    a_i = a_j;
    a_j = temp;

   return (*this);  
 }


bool DKTS::writeAllNorms(ostream& s) const
 {
    bool value = true;
    
    const LongInt d = DKTS::d();
    const LongInt r = DKTS::k();
    
    for(LongInt mu=0; mu<d; mu++)
     {
        s << "[" << mu << "] : ";
        
        for(LongInt j=0; j<r; j++)
         {
            s << l2((*this)(j, mu)) << ", ";
         }
        s << endl;
     }
 
   return value;
 }
 

LongReal frobeniusNorm (const DKTS& A)
 {
    LongReal value   = 0.0;
    LongReal product = 1.0;
    LongReal sum     = 0.0;
    
    const LongInt d = A.d();
    const LongInt r = A.k();
    
    for(LongInt i1=0; i1<r; i1++)
     {
        product = 1.0;
        for(LongInt mu=0; mu<d; mu++)
         {
            const DKTVector& a_imu = A(i1, mu);
            
            product *= innerProduct(a_imu, a_imu);
         }
        
        value += product;
             
        sum = 0.0;
        for(LongInt i2=i1+1; i2<r; i2++)
         {
            product = 1.0;
            for(LongInt mu=0; mu<d; mu++)
             {
                product *= innerProduct(A(i1, mu), A(i2, mu));
             }
            sum += product;
         }
        
        value += 2.0*sum;
     }             

   return sqrt(fabs(value));
 }


LongReal DKTS::frobeniusNormOfSummand(const LongInt& i) const
 {
    const LongInt d = DKTS::d();

    LongReal value = 1.0;

    for(LongInt mu=0; mu<d; mu++)
     {
        const DKTVector& v = (*this)(i, mu);
        value *= l2(v);
     }   

   return value;
 }


LongReal innerProduct (const DKTS& A, const DKTS& B)
 {
    const LongInt k = A.k();
    const LongInt l = B.k();
    const LongInt d = A.d();

    LongReal value = 0.0, factor_ij = 1.0;
     
    for(LongInt i=0; i<k; i++)
     {
        for(LongInt j=0; j<l; j++)
         {
            factor_ij = 1.0;
            
            for(LongInt nu=0; nu<d; nu++)
             {
                factor_ij *= innerProduct(A(i,nu), B(j,nu));
             }
            value += factor_ij;
         }
     }      

   return value;
 }


DKTS& DKTS::readDataFrom(const IString& file)
 {
    ifstream from(file);

    LongInt d=0, k=0, n=0, n1=0;

    char t[64];    

    from.setf(ios::scientific, ios::floatfield);
    
    from >> t >> t >> t >> t;
    from >> d;

    for(LongInt mu=0; mu<d; mu++)
     {
       from >> t >> t;
       from >> n1;
       n = MAX(n, n1);
     }

    from >> t >> t >> t;
    from >> k;     

cout << "file = " << file << endl;

cout << "d = " << d << endl;
cout << "k = " << k << endl;
cout << "n = " << n << endl;


    resize(d, k, n);

    for(LongInt j=0; j<k; j++)
     {
        for(LongInt mu=0; mu<d; mu++)
         {
            DKTVector& v = (*this)(j,mu);

            for(LongInt i=0; i<n; i++)
             {
                LongReal value;
                from >> value;
                
                v(i) = value;
             }
         }
     }

   return (*this);
 }


bool DKTS::writeDataTo(const IString& file) const
 {
    const LongInt d = DKTS::d();
    const LongInt k = DKTS::k();
    const LongInt n = DKTS::n();

    ofstream to(file);

    to << "Tensor in d = " << d << endl;
    
    for(LongInt mu=0; mu<d; mu++)
     {
        to << "n[" << mu << "] = " << n << endl;
     }

    to << "Rank k = " << k << endl;

    to.setf(ios::scientific, ios::floatfield);
    to.precision(20);
    
    for(LongInt j=0; j<k; j++)
     {
        for(LongInt mu=0; mu<d; mu++)
         {
            DKTVector& v = (*this)(j,mu);

            for(LongInt i=0; i<n; i++)
             {
                to << v(i) << endl;
             }
         }
     }

   return true;
 }


DKTS& DKTS::readDataFrom(const IString& file, const LongInt& k1)
 {
    ifstream from(file);

    LongInt d=0, k=0, n=0, n1=0;
  
    char t[64];

    from >> t >> t >> t >> t;
    from >> d;

    for(LongInt mu=0; mu<d; mu++)
     {    
        from >> t >> t;
        from >> n1;
        n = MAX(n, n1);
     }

    from >> t >> t >> t;
    from >> k;     

    resize(d, MAX(k,k1), n);

    k = MIN(k,k1);

    for(LongInt j=0; j<k; j++)
     {
        for(LongInt mu=0; mu<d; mu++)
         {
            DKTVector& v = (*this)(j,mu);

            for(LongInt i=0; i<n; i++)
             {
                from >> v(i);
             }
         }
     }

   return (*this);
 }
 

DKTS& DKTS::reScaled()
 {
    const LongInt d = DKTS::d();
    const LongInt k = DKTS::k();

    LongReal value = 1.0;

    for(LongInt j=0; j<k; j++)
     {
        value = frobeniusNormOfSummand(j);

        for(LongInt mu=0; mu<d; mu++)
         {
            DKTVector& x = (*this)(j,mu);
            
            x.normalized();

            if(mu==0)
             {
                x *= value;
             }
         }
     }

   return (*this);
 }
 
 
DKTS& DKTS::setRand(const LongReal eps)
 {
    const LongInt d = DKTS::d();
    const LongInt k = DKTS::k();    
    
    (*this)().setRand(eps);     
    
   return (*this);
 }


DKTS& DKTS::setNull()
 {
    (*this)().setNull();
 
   return (*this);
 }

 
//edit 2010-01-23
DKTS& DKTS::setAllEntries(const LongReal& a)
 {
    const LongInt d = (*this).d();
    const LongInt n = (*this).n();
    
    const LongReal d_inv = 1.0/((LongReal)d);
    
    (*this)
     .resize(d, 1, n)
    ;
    
    (*this)().setAllEntrys(pow(a, d_inv));
    
   return(*this);
 }
//end edit
 
 
DKTS& DKTS::scale(const LongReal eps, const LongReal& alpha)
 {        
    const LongInt d = DKTS::d();
    const LongInt k = DKTS::k();    

    LongReal factor = eps;

    for(LongInt j=0; j<k; j++)    
     {        
        const LongReal factorM = exp(log(alpha*factor)/LongReal(d));

        for(LongInt mu=0; mu<d; mu++)
         {
            (*this)(j,mu) *= factorM;
         }
        factor *= eps;
     }
        
   return (*this);
 }

DKTS& DKTS::balanced()
 {
    const LongInt d = DKTS::d();
    const LongInt k = DKTS::k();    

    for(LongInt j=0; j<k; j++)
     {
        LongReal normA = 1.0;
        LongInt mu = 0;
        
        for(mu=0; mu<d; mu++)
         {
            DKTVector& a = (*this)(j,mu);
            normA *= a.normalized();            
         }
        
        LongReal fac = pow(normA, 1.0/((LongReal)d));
        
        for(mu=0; mu<d; mu++)
         {
            (*this)(j,mu) *= fac;
         }        
     }
    
   return (*this);
 }
 
 
DKTS& DKTS::regeneratedBy(const DKTS& x, const DKTS& Z)
 {
    const LongInt d = Z.d();
    const LongInt k = x.k();
    const LongInt l = Z.k();
    const LongInt n = Z.n();

    resize(d, k, n);

    for(LongInt mu=0; mu<d; mu++)
     {
        for(LongInt j=0; j<k; j++)
         {
            DKTVector& v = (*this)(j,mu);
            DKTVector& a = x(j,mu);
        
            for(LongInt i=0; i<l; i++)
             {
                v.add(a(i), Z(i, mu), 1.0);
             }           
         }
     }  
 
   return (*this);
 } 
 
 
DKTS& DKTS::setOrthogonal2(const DKTS& A, const LongReal& error)
 { 
    DKTS B(A);
          
    integer m     = B.n();
    integer n     = B.k(); 
    integer dim   = MIN(n,m);
    integer lda   = m;
    integer lwork = 3*n+1;
    integer info  = 0;
        
    integer* jpvt = (integer*) new integer[n];
    
    B.balanced();

    for(LongInt i1=0; i1<n; i1++)
     {
        jpvt[i1] = 0;
        //B(i1,0).normalized();
     }

    LongRealPointer work = (LongRealPointer) new LongReal[lwork]; 
    LongRealPointer tau  = (LongRealPointer) new LongReal[dim]; 

    const LongInt d = A.d();

    LongInt R = 0;
    for(LongInt mu=0; mu<d; mu++)
     {        
        LongRealPointer A_mu = &(B(0, mu, 0));        

        info = TensorCalculus::Lapack<double>::geqp3(m, n, A_mu, lda, jpvt, tau, work, lwork);
        
        LongInt muR = 0;                
  
        for(LongInt i=0; i<n; i++)
         {
            jpvt[i] = 0;

            if(i<dim)
             {
                const LongReal ve = A_mu[i + i*m];
                
                if(fabs(ve)>error)
                 {
                    muR ++;
                 }                
             }
            
         }  
                
        R = MAX(R, muR);
        
        integer mQ   = m;
        integer nQ   = dim;
        integer dimQ = nQ;

        info = TensorCalculus::Lapack<double>::orgqr(mQ, nQ, dimQ, A_mu, lda, tau, work, lwork);
        
     }
         
    (*this)
     .resize(d, R, m)
    ;
    
    dim = m*R;
    integer inc = 1;
    for(LongInt mu=0; mu<d; mu++)
     {        
        LongRealPointer B_mu = &(B(0, mu, 0));
        LongRealPointer A_mu = &((*this)(0, mu, 0));

        TensorCalculus::Blas<double>::copy(dim, B_mu, inc, A_mu, inc);

     }
       
    delete [] work;
    delete [] tau;
    
    delete [] jpvt;    
    
   return (*this);
 }


DKTS& DKTS::setCoefficientsSystemOf(const DKTS& A, const DKTS& Z)
 {
    const LongInt d = A.d();
    const LongInt k = A.k();
    const LongInt l = Z.k();

    resize(d, k, l);

    for(LongInt mu=0; mu<d; mu++)
     {
        for(LongInt j=0; j<k; j++)
         {
            DKTVector& a  = A(j,mu);
            DKTVector& al = (*this)(j,mu);

            for(LongInt i=0; i<l; i++)
             {
                al(i) = innerProduct(Z(i,mu), a);
             }
         }
     }

   return (*this);
 }


LongReal DKTS::distanceTo(const DKTS& A) const
 {
    LongReal value = 0.0;

    value = innerProduct((*this), (*this)) - 2.0*innerProduct((*this), A) + innerProduct(A, A);

   return sqrt(fabs(value));
 }


LongReal DKTS::distanceTo(const LongInt& j1, const LongInt& j2) const
 {
    const LongInt d = DKTS::d();
    
    LongReal value = -1.0;
    
    LongReal normj1 = frobeniusNormOfSummand(j1);
    LongReal normj2 = frobeniusNormOfSummand(j2);
    LongReal skp    = 1.0;
    
    for(LongInt mu=0; mu<d; mu++)
     {        
        skp    *= innerProduct((*this)(j1, mu), (*this)(j2, mu));        
     }
    
    value = sqrt(fabs(normj1*normj1 - 2.0*skp + normj2*normj2));
    
    cout << value << " : " << sqrt(normj1*normj1 + 2.0*skp + normj2*normj2) << endl;
    
   return value;
 }
 

DKTS& DKTS::setProductOf(const LongReal& alpha, const DKTS& A1, const DKTS& A2, const LongReal& beta)
 {
    const LongInt d  = A1.d();
    const LongInt n  = A1.n();

    const LongInt k1 = A1.k();
    const LongInt k2 = A2.k();

    const LongInt  k = k1*k2;
    
    LongInt m = (LongInt) sqrt((double)n);
        
    resize(d, k, n);

    LongInt j = 0;

    for(LongInt mu=0; mu<d; mu++)
     {
        for(LongInt j1=0; j1<k1; j1++)
         {
            for(LongInt j2=0; j2<k2; j2++)
             {
                j = j1 + j2*k1;
    
                TensorCalculus::Blas<double>::gemm(const_notConjTrans, const_notConjTrans,
                      m, m, m,
                      alpha, &A1(j1, mu)(0),
                      m, &A2(j2, mu)(0),
                      m, beta, &(*this)(j, mu)(0), m);
             }
         }
     }    

    
   return (*this); 
 } 


bool DKTS::evaluateAt(const LongReal& alpha, const DKTS& v, const LongReal& beta, DKTS& w) const
 {
    bool value = true;

    const LongInt d  = v.d();    
    const LongInt k1 = (*this).k();    
    const LongInt k2 = v.k();

    const LongInt  k = k1*k2;

    LongInt m  = v.n();            
    
    w.resize(d, k, m);

    LongInt j = 0;

    for(LongInt mu=0; mu<d; mu++)
     {
        for(LongInt j1=0; j1<k1; j1++)
         {
            for(LongInt j2=0; j2<k2; j2++)
             {
                j = j1 + j2*k1;
         
                TensorCalculus::Blas<double>::gemv( const_notConjTrans, m, m,
                       alpha, &(*this)(j1, mu)(0), m,
                       &v(j2, mu)(0), const_inc, beta, &w(j, mu)(0), const_inc);
             }
         }
     }    
 
   return value;
 }

DKTS& DKTS::setSumOf (const LongReal& a, const DKTS& A1, const LongReal& b, const DKTS& A2, const LongReal& c, const DKTS& A3)
 {
    const LongInt d1 = A1.d();
    const LongInt d2 = A2.d();
    const LongInt d3 = A3.d();

    const LongInt d  = MAX(MAX(d1, d2),d3);
    const LongInt n  = MAX(MAX(A1.n(), A2.n()), A3.n());

    const LongInt k1 = A1.k();
    const LongInt k2 = A2.k();
    const LongInt k3 = A3.k();

    const LongInt  k = k1 + k2 + k3;

    resize(d, k, n);
    
    LongInt j=0;

    // mu=0
    for(j=0; j<k1; j++)
     {
        ((*this)(j, 0)).setCopyOf(a, A1(j, 0));
     }

    for(j=0; j<k2; j++)
     {
        ((*this)(j+k1, 0)).setCopyOf(b, A2(j, 0));
     }

    for(j=0; j<k3; j++)
     {
        ((*this)(j+k1+k2, 0)).setCopyOf(c, A3(j, 0));
     }


    for(LongInt mu=1; mu<d; mu++)
     {
        for(j=0; j<k1; j++)
         {
            (*this)(j,mu)       = A1(j,mu);
         }

        for(j=0; j<k2; j++)
         {
            (*this)(j+k1,mu)    = A2(j,mu);
         }
         
        for(j=0; j<k3; j++)
         {
            (*this)(j+k1+k2,mu) = A3(j,mu);
         }         
     }
 
   return (*this); 
 }

DKTS& DKTS::setSumOf(const LongReal& alpha, const DKTS& A1, const LongReal& beta, const DKTS& A2)
 {
    const LongInt d1 = A1.d();
    const LongInt d2 = A2.d();

    const LongInt d  = MAX(d1, d2);
    const LongInt n  = MAX(A1.n(), A2.n());

    const LongInt k1 = A1.k();
    const LongInt k2 = A2.k();

    const LongInt  k = k1 + k2;

    resize(d, k, n);
    
    LongInt j=0;

    // mu=0
    for(j=0; j<k1; j++)
     {
        ((*this)(j, 0)).setCopyOf(alpha, A1(j, 0));
     }

    for(j=0; j<k2; j++)
     {
        ((*this)(j+k1, 0)).setCopyOf(beta, A2(j, 0));
     }

    for(LongInt mu=1; mu<d; mu++)
     {
        for(j=0; j<k1; j++)
         {
            (*this)(j,mu) = A1(j,mu);
         }

        for(j=0; j<k2; j++)
         {
            (*this)(j+k1,mu) = A2(j,mu);
         }
     }

   return (*this);
 }


bool DKTS::writeDiagonalDataTo(const IString& file)const
 {
    bool value = true;
    
    const LongInt d = DKTS::d();
    const LongInt r = DKTS::k();
    const LongInt n = DKTS::n();

    ofstream to(file);
    
    to << n << endl;
    
    LongReal factor=1.0, sum=0.0;
    
    for(LongInt i=0; i<n; i++)
     {
        sum = 0.0;
        for(LongInt j=0; j<r; j++)
         {
            factor = 1.0;
            for(LongInt mu=0; mu<d; mu++)
             {
                factor *= (*this)(j, mu)(i);
             }
            sum += factor;
         }
        to << i << '\t' << sum << endl;
     }
        
   return value;
 }


bool DKTS::writeIndexSet(ostream& s) const
 {
    const LongInt k = DKTS::k();

    LongReal cond = 0.0, ni = 0.0;
 
    for(LongInt i=0; i<k; i++)
     {
        ni = frobeniusNormOfSummand(i);
        s << i << " : " << ni << endl;
        
        cond += ni;
     }

    cond = sqrt(cond);
    
    const LongReal normA = frobeniusNorm((*this));
    s << endl << "|A|       = " <<  normA << endl;
    s         << "cond of A = " << cond/frobeniusNorm((*this)) << endl <<endl;
   return true;
 }


bool DKTS::writeIndexSet(DescriptionData& s) const
 {
    const LongInt k = DKTS::k();

    LongReal cond = 0.0, ni = 0.0;

    for(LongInt i=0; i<k; i++)
     {
        ni = frobeniusNormOfSummand(i);
        s.addString(i);
        s.addString(" : ");
        s.addString(ni);
        s.addString("\n");

        cond += ni;
     }

    cond = sqrt(cond);

    const LongReal normA = frobeniusNorm((*this));

    s.addString("\n");
    s.addString("|A|       = ");
    s.addString(normA);
    s.addString("cond of A = ");
    s.addString(cond/frobeniusNorm((*this)));
    s.addString("\n\n");
   return true;
 }


bool DKTS::writeParameter(ostream& s) const
 {
    const LongInt d = DKTS::d();
    const LongInt k = DKTS::k();
    const LongInt n = DKTS::n();

    s << "d = " << d << endl;
    s << "k = " << k << endl;
    s << "n = " << n << endl;

   return true;
 }


LongReal DKTS::operator () (const IntegerArray& l) const
 {
    const LongInt d = DKTS::d();
    const LongInt r = DKTS::k();        
    
    LongReal value = 0.0;
    LongReal factor = 1.0;    

    for(LongInt j=0; j<r; j++)
     {
        factor = 1.0;
        for(LongInt mu=0; mu<d; mu++)
         {  
            const LongReal entry = (*this)(j, mu, l(mu));
        
            factor *= entry;
         }        
        value += factor;
     }    

   return value;
 }


LongReal DKTS::operator () (const IntegerArray& l, const LongInt& jMax, const LongInt& jMin) const
 {
    const LongInt d = DKTS::d();
    const LongInt r = DKTS::k();        
    
    LongReal value = 0.0;
    LongReal factor = 1.0;
            
    for(LongInt j=jMin; j<jMax; j++)
     {
        factor = 1.0;
        for(LongInt mu=0; mu<d; mu++)
         {  
            const LongReal entry = (*this)(j, mu, l(mu));
        
            factor *= entry;
         }
        value += factor;
     }    

   return value;
 }


DKTS& DKTS::setCrossApproximationOf(const DKTS& A, const LongInt& maxSteps)
 {
    const LongInt d = A.d();
    const LongInt n = A.n();
    
    IntegerArray i(d, n);
    
    setCrossApproximationOf(A, i, maxSteps);
    
   return (*this);
 }


DKTS& DKTS::setCrossApproximationOf(const DKTS& A, IntegerArray& i, const LongInt& maxSteps)
 {
    const LongInt d = A.d();
    const LongInt R = A.k();
    const LongInt n = A.n();
    
    (*this).resize(d, 1, n);
  
    LongReal pivot  = fabs(A(i));

    // find a good pivot
    for(LongInt l=0; l<R; l++)
     {  
        IntegerArray j(d);
        
        LongInt maxI = -1;
        
        for(LongInt mu=0; mu<d; mu++)
         {                        
            j(mu) = (A(l, mu)).indexOfMax();            
         }
        
        const LongReal newValue = fabs(A(j));
        
        if(pivot<newValue)
         {
            pivot = newValue;
            i=j;
         }
     }
       
    // improving Pivot    
    for(LongInt l=0; l<maxSteps; l++)
     {
        for(LongInt mu=0; mu<d; mu++)
         {        
            LongInt& i_m = i(mu);
            LongInt iMuP = i_m;        
        
            for(LongInt m=0; m<n; m++)
             {
                i_m  = m;
                
                const LongReal value = fabs(A(i));
                
                if(fabs(pivot) < value)
                 {
                    iMuP = m;
                    pivot = value;
                 }
             }
             
            i_m = iMuP;
         }
     }

    if(fabs(pivot)<EPS_NULL)
     {
        cout << "pW";
     }
    // computing Approximation
    
    for(LongInt mu=0; mu<d; mu++)
     {
        DKTVector& x = (*this)(0, mu);
        
        LongInt& i_m = i(mu);
        const LongInt iMu = i_m;        
        
        for(LongInt m=0; m<n; m++)
         {
            i_m  = m;            
            
            if(0<mu)
             {
                x(m) = A(i)/pivot;
             }
            else
             {
                x(m) = A(i);
             }
         }
        i_m = iMu;     
     }          
     
     const LongReal g1 = innerProduct((*this), A);
     const LongReal g2 = innerProduct((*this), (*this));
          
     (*this) *= g1/g2;     
 
   return (*this);
 }


DKTS& DKTS::setAlsApproximationOf(const DKTS& A, const LongInt& maxSteps, const bool& useCross)
 {
    const LongInt d = A.d();
    const LongInt R = A.k();
    const LongInt n = A.n();    

    RMatrix S(R,d);
    RVector s(R);

    if(useCross)
     {
        (*this)
         .setCrossApproximationOf(A, 5)
        ;
     }    
        
    for(LongInt mu=0; mu<d; mu++)
     {
        DKTVector& x = (*this)(0,mu);
        x.normalized();
        
        for(LongInt i=0; i<R; i++)
         {
            S(i, mu) = innerProduct(x, A(i,mu));
         }
     }
    
    LongReal temp = 1.0;
        
    for(LongInt k=0; k<maxSteps; k++)
     {
        for(LongInt mu=0; mu<d; mu++)
         {
            DKTVector& x = (*this)(0,mu);
                                    
            LongInt nu = 0;
            
            temp = 1.0;
            for(nu=0; nu<d; nu++)
             {
                if(nu!=mu)
                 {
                    temp *= S(0,nu);
                 }
             }
            
            x  = A(0, mu);
            x *= temp;
                                    
            LongInt i = 0;
            
            for(i=1; i<R; i++)
             {
                temp = 1.0;
                for(LongInt nu=0; nu<d; nu++)
                 {
                    if(nu!=mu)
                     {
                        temp *= S(i,nu);
                     }
                 }
                 
                x.update(temp, A(i,mu));
             }
            
            temp = 1.0;
            for(nu=0; nu<d; nu++)
             {
                if(nu!=mu)
                 {
                    temp *= innerProduct((*this)(0,nu), (*this)(0,nu));
                 }
             }            
            
            x /= temp;
            
            for(i=0; i<R; i++)
             {
                S(i, mu) = innerProduct(x, A(i,mu));
             }
         }     
     }     
    
    const LongReal g1 = innerProduct((*this), A);
    const LongReal g2 = innerProduct((*this), (*this));
          
    (*this) *= g1/g2;     
      
   return (*this);
 }
 
 
DKTS& DKTS::improveApproximation(const DKTS& A, const LongInt& maxStepsAdd, const LongInt& maxStepsAls)
 {
    const LongInt d = A.d();
    const LongInt R = A.k();
    const LongInt r = (*this).k();
    const LongInt n = A.n();    

    RMatrix S(R,d), T(r,d);
    RVector s(R);

    for(LongInt l=0; l<maxStepsAdd; l++)
     {    
        for(LongInt j1=0; j1<r; j1++)
         {        
            for(LongInt mu=0; mu<d; mu++)
             {
                DKTVector& x = (*this)(j1,mu);                
        
                for(LongInt i=0; i<R; i++)
                 {
                    S(i, mu) = innerProduct(x, A(i,mu));
                 }
                 
                for(LongInt j=0; j<r; j++)
                 {
                    T(j, mu) = innerProduct(x, (*this)(j,mu));
                 }                
             }
    
             LongReal temp = 1.0;
        
             for(LongInt k=0; k<maxStepsAls; k++)
              {
                 for(LongInt mu=0; mu<d; mu++)
                  {
                     DKTVector& x = (*this)(j1,mu);
                                    
                     LongInt nu = 0;
            
                     temp = 1.0;
                     for(nu=0; nu<d; nu++)
                      {
                         if(nu!=mu)
                          {
                             temp *= S(0,nu);
                          }
                      }
            
                         x  = A(0, mu);
                         x *= temp;
                                    
                         LongInt i = 0, j = 0;
            
                         for(i=1; i<R; i++)
                          {
                             temp = 1.0;
                             for(LongInt nu=0; nu<d; nu++)
                              {
                                 if(nu!=mu)
                                  {
                                     temp *= S(i,nu);
                                  }
                              }
                 
                             x.update(temp, A(i,mu));
                          }
            
                         for(j=0; j<r; j++)
                          {
                             if(j1!=j)
                              {
                                 temp = 1.0;
                                 for(LongInt nu=0; nu<d; nu++)
                                  {
                                     if(nu!=mu)
                                      {
                                         temp *= T(j,nu);
                                      }
                                  }
                                                                    
                                 x.update(-1.0*temp, (*this)(j,mu)); 
                              }                                             
                          }            
            
            
                         temp = 1.0;
                         for(nu=0; nu<d; nu++)
                          {
                             if(nu!=mu)
                              {
                                 temp *= innerProduct((*this)(j1,nu), (*this)(j1,nu));
                              }
                          }                        
                         x /= temp;
            
                         for(i=0; i<R; i++)
                          {
                             S(i, mu) = innerProduct(x, A(i,mu));
                          }
                          
                         for(LongInt j=0; j<r; j++)
                          {
                             T(j, mu) = innerProduct(x, (*this)(j,mu));
                          }
                          
                  }// End for(LongInt mu=0; mu<d; mu++)     
              }// End for(LongInt k=0; k<maxSteps2; k++)     
         }// End for(LongInt j1=0; j1<r; j1++)
     } // End for(LongInt l=0; l<maxSteps1; l++)
       
   return (*this);
 }
 
 
DKTS& DKTS::setModellExample(const LongInt& d, const LongInt& n)
 {
    (*this)
     .resize(d, d, n)
    ;        
    
    (*this)().setAllEntrys(1.0);
    
    DKTVector& a = (*this)(0,0);
    
    for(LongInt i=0; i<n; i++)
     {
        a(i) = (LongReal)(i+1);
     }
    
    for(LongInt mu=1; mu<d; mu++)
     {
        (*this)(mu,mu) = a;
     }    
 
   return (*this);
 }
 
 
bool DKTS::analyseVectorSystem() const
 {
    bool value=true;

    char job  = 'N';
    char uplo = 'L';
    
    LongInt d = DKTS::d(); 
    LongInt r = DKTS::k(); 
    LongInt n = DKTS::n(); 
  
    LongReal eins = 1.0;

    LongInt lwork = 3*r-1;
    int     info  = 0;    

    cout << setprecision(3) << scientific;
    
    for(LongInt mu=0; mu<d; mu++)
     {
        RMatrix G(r,r);
        RVector eig(r);
        RVector work(lwork);
         
        LongReal& x_mu = (*this)(mu)(0);
        
        TensorCalculus::Blas<double>::gemm(const_trans, const_notConjTrans, r, r, n, eins, &x_mu, n, &x_mu, n, eins, &G(0,0), r);
        
        info = TensorCalculus::Blas<double>::syev (job, uplo, r, &G(0,0), r, &eig(0), &work(0), lwork);
        
        cout << "[" << mu << "] : ";
                        
        for(LongInt j=0; j<r; j++)
         {
            cout << eig(j) << ", ";
         }
        cout << endl;
     }    
   
   return value; 
 }
 
 
bool DKTS::analyseGSystem() const 
 {
    bool value=true;

    char   job  = 'N';
    char   uplo = 'L';
    
    LongInt d = DKTS::d(); 
    LongInt r = DKTS::k(); 
  
    LongReal eins = 1.0;

    LongInt lwork = 3*r-1;
    int     info  = 0;    

    cout << setprecision(3) << scientific;
    
    for(LongInt mu1=0; mu1<d; mu1++)
     {
        RMatrix G(r,r);
        RVector eig(r);
        RVector work(lwork);

        for(LongInt j1=0; j1<r; j1++)
         {
            for(LongInt j2=0; j2<r; j2++)
             {
                LongReal& g=G(j1,j2);
                
                g = 1.0;
                for(LongInt mu=0; mu<d; mu++)
                 {
                    if(mu!=mu1)
                     {
                       g *= innerProduct((*this)(j1,mu), (*this)(j2,mu));
                     }
                 }
             }
         }         
        
        info = TensorCalculus::Blas<double>::syev(job, uplo, r, &G(0,0), r, &eig(0), &work(0), lwork);
        
        cout << "[" << mu1 << "] : ";
                        
        for(LongInt j=0; j<r; j++)
         {
            cout << eig(j) << ", ";
         }
        cout << endl;
     }    
   
   return value; 
 }


bool DKTS::evaluateOrthProjAt(const DKTS& v, DKTS& w) const
 {
    bool value = true;
    
    LongInt d = v.d();
    LongInt r = v.k();
    
    LongInt n = (*this).n();
    LongInt k = (*this).k();
    
    w.resize(d, r, n);
    
    LongReal eins  =  1.0;
    LongReal nuLL  =  0.0;
    
    for(LongInt mu=0; mu<d; mu++)
     {
        RVector temp(k);
        
        LongRealPointer U = &(*this)(mu)(0);
        LongRealPointer t = &temp(0);
                
        for(LongInt j=0; j<r; j++)
         {                        
            LongRealPointer vV = &v(j,mu)(0);
            LongRealPointer wV = &w(j,mu)(0);            

            TensorCalculus::Blas<double>::gemv(const_trans, n, k, eins, U, n, vV, const_inc, nuLL, t, const_inc);
            TensorCalculus::Blas<double>::gemv(const_notConjTrans, n, k, eins, U, n, t, const_inc, nuLL, wV, const_inc);
         }
     }    
    
   return value;
 }


bool DKTS::evaluateOrthCompProjAt (const DKTS& v, DKTS& w) const
 {
    bool value = true;
    
    LongInt d = v.d();
    LongInt r = v.k();

    LongInt n = (*this).n();
    LongInt k = (*this).k();    
    
    w.resize(d, r, n);

    w()=v();
     
    LongReal eins  =  1.0;
    LongReal mEins = -1.0;
    LongReal nuLL  =  0.0;
    
    for(LongInt mu=0; mu<d; mu++)
     {
        RVector temp(k);
        
        LongRealPointer U = &(*this)(0,mu)(0);
        LongRealPointer t = &temp(0);
                
        for(LongInt j=0; j<r; j++)
         {                        
            LongRealPointer vV = &v(j,mu)(0);
            LongRealPointer wV = &w(j,mu)(0);                        
            
            TensorCalculus::Blas<double>::gemv(const_trans, n, k, eins, U, n, vV, const_inc, nuLL, t, const_inc);

            TensorCalculus::Blas<double>::gemv(const_notConjTrans, n, k, mEins, U, n, t, const_inc, eins, wV, const_inc);
         }
     }    
           
   return value;
 }


DKTS& DKTS::normalizedAddens()
 {
    const LongInt d = DKTS::d();
    const LongInt r = DKTS::k();
 
    for(LongInt mu=0; mu<d; mu++)
     {
        for(LongInt j=0; j<r; j++)
         {
            (*this)(j,mu).normalized();
         }
     }
 
   return (*this);
 }
 
 
LongInt DKTS::addR1Addent(const DKTS& x) 
 {
    LongInt value = -1;
     
    const LongInt r =  DKTS::k();
    const LongInt d =  DKTS::d();
      
    const LongReal normX = frobeniusNorm(x);
        
    LongReal dotProd = 1.0;
    LongReal normA   = 1.0;    
    LongReal waMin   = 1.0;    
    LongReal lambda  = 1.0;
    
    for(LongInt j=0; j<r; j++)
     {             
        normA = (*this).frobeniusNormOfSummand(j);
        
        dotProd = 1.0;        
        for(LongInt mu=0; mu<d; mu++)
         {
            dotProd *= innerProduct((*this)(j,mu), x(0,mu));
         }

        LongReal wa = fabs(normA*normA*normX*normX - dotProd*dotProd);

        if(wa<1.0e-10 && wa<=waMin)
         {
            waMin  = wa;
            value  = j;
            lambda = dotProd/normA;
         }
     }

    if(value!=-1)
     {
        (*this)(value, 0) *= 1.0+lambda;
     }
    else
     {
        DKTS xt((*this));        
	       (*this).setSumOf(1.0, xt, 1.0, x);
     }
  
   return value;
 }
 
 
LongReal maximumValueOf(const DKTS& a, IntegerArray& indexOfMaximum)
 {  
    LongReal value = 1.0;
    const LongInt d = a.d();
    const LongInt m = a.n();
            
    indexOfMaximum.resize(d);
    
    const LongReal a1 = lowerBoundOfMinimumValue(a);
    const LongReal a2 = upperBoundOfMaximumValue(a);
    
    cout << "a1 = " << a1 << endl;
    cout << "a2 = " << a2 << endl;
    
    DKTS b(a);
    
    if(a1 < 0.0)
     {
        
        DKTS e(d, 1, m);
                
        e().setAllEntrys(1.0);
        
        b.setSumOf(1.0, a, (fabs(a1)+100), e);
     }
    
    b.balanced();
    computeMaximumValueOf(b, indexOfMaximum);
    
    value = a(indexOfMaximum);
       
   return value;
 }

 
LongReal computeMaximumValueOf (const DKTS& a, IntegerArray& indexOfMaximum) 
 {  
    LongReal value = 1.0;
        
    LongReal factor = 1.0;
    
    const LongInt d = a.d();
    const LongInt r = a.k();
    const LongInt m = a.n();
    
    DKTS x(d, 1, m), w(d, 1, m);
    
    x().setAllEntrys(1.0/sqrt((LongReal)m));        

    RVector v(r);
    
for(LongInt runs=0; runs<16; runs++)
{
    for(LongInt mu1=0; mu1<d; mu1++)
     {        
        for(LongInt j=0; j<r; j++)
         {
            LongReal& v_j = v(j);
            
            v_j = 1.0;
            for(LongInt mu=0; mu<d; mu++)
             {
                if(mu!=mu1)
                 {
                        LongReal innerValue = 0.0;
                        
                        for(LongInt l=0; l<m; l++)
                         {
                            innerValue += a(j,mu,l)*x(0,mu,l)*x(0,mu,l);
                         }
                    
                    v_j *= innerValue; //innerProduct(a(j,mu), x(0,mu));
                 }
             }             
         }
        
        TensorCalculus::Blas<double>::gemv (DKTS::const_notConjTrans, m, r, factor, &a(mu1)(0), m, &v(0), DKTS::const_inc,
               DKTS::const_null, &w(0, mu1)(0), DKTS::const_inc);
        
        
        for(LongInt l=0; l<m; l++)
         {
            x(0, mu1, l) *= w(0, mu1, l);
         }
        
        x(0,mu1).normalized();
        LongInt index_mu1;
        maximumValueOf(x(0, mu1), index_mu1);        
        /*        
        x(0, mu1).setAllEntrys(0.0);
        x(0, mu1, index_mu1) = 1.0;
        */
        indexOfMaximum(mu1) = index_mu1;
     }                 
}    
   // cout << x;    
   return value;
 }

 
LongReal minimumValueOf(const DKTS& a, IntegerArray& indexOfMinimum)
 {  
    LongReal value = 1.0;
    
    DKTS b(a);
    
    b *= -1.0;
     
    maximumValueOf(b, indexOfMinimum);
    
    value = a(indexOfMinimum);
    
   return value;
 } 

 
LongReal upperBoundOfMaximumValue(const DKTS& a)
 {         
    const LongInt d = a.d();
    const LongInt r = a.k();
    
    LongReal value = 0.0;    
    for(LongInt j=0; j<r; j++) 
     {
        LongReal valueJ = 1.0;        
        for(LongInt mu=0; mu<d; mu++)
         {
            LongInt indexMu;
            
            valueJ *= maximumValueOf(a(j, mu), indexMu);            
         }
        
        value += valueJ;
     }
 
   return value;
 } 


LongReal lowerBoundOfMinimumValue(const DKTS& a)
 {         
    const LongInt d = a.d();
    const LongInt r = a.k();
    
    LongReal value = 0.0;    
    for(LongInt j=0; j<r; j++) 
     {
        LongReal valueJ = 1.0;        
        for(LongInt mu=0; mu<d; mu++)
         {
            LongInt indexMu;
            
            valueJ *= minimumValueOf(a(j, mu), indexMu);            
         }
        
        value += valueJ;
     }
 
   return value;
 } 

