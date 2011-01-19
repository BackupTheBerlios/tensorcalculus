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

// RkRMatrix.cpp: Implementierung der Klasse RkRMatrix.
//
//////////////////////////////////////////////////////////////////////

#include "RkRMatrix.hpp"
#include "HMatrixInterface.hpp"
#include "BlasInterface.hpp"

RkRMatrix::RkRMatrix(const LongInt m, const LongInt n, const LongInt k)
:ARMatrix(m, n), attr_A(m, MIN(n, k)), attr_B(n, MIN(n, k))
 {
    (*this)
     .setRank(MIN(n, k))
     .setType(isRkRMatrix)
    ; 
 }


RkRMatrix::RkRMatrix (const LongInt n, const LongInt m, const LongInt k, const bool& flag)
:ARMatrix(m, n), attr_A(m, k), attr_B(n, k)
 {
    (*this)
     .setRank(k)
     .setType(isRkRMatrix)
    ; 
 }


RkRMatrix::RkRMatrix(const RkRMatrix& aRkMatrix)
:ARMatrix(aRkMatrix.numberOfRows(), aRkMatrix.numberOfColumns()), 
 attr_A(aRkMatrix.A()), attr_B(aRkMatrix.B())
 {
    (*this)
     .setRank(aRkMatrix.rank())
     .setType(isRkRMatrix)
    ;
 }


RkRMatrix::RkRMatrix (const RMatrix& Z, const LongInt k)
:ARMatrix(Z.numberOfRows(), Z.numberOfColumns()), 
 attr_A(Z.numberOfRows(), MIN(Z.numberOfColumns(), k)), attr_B(Z.numberOfColumns(), MIN(Z.numberOfColumns(), k))
 {
    (*this)
     .setRank(MIN(Z.numberOfColumns(), k))
     .setType(isRkRMatrix)
    ;

    RMatrix Temp(Z);

    LongRealSvd svd(Temp, rank());
 
    svd.multiplyUSigma();
    attr_A = svd.attr_U;
    attr_B = svd.attr_V;

 }

RkRMatrix::~RkRMatrix()
 {
 }


RkRMatrix& RkRMatrix::operator = (const RkRMatrix& R)
 {
    (*this)
     .setRank(MIN(R.numberOfColumns(), R.rank()))
     .setNumberOfRows(R.numberOfRows())
     .setNumberOfColumns(R.numberOfColumns())
     .setType(isRkRMatrix)
    ; 

    attr_A = R.attr_A;
    attr_B = R.attr_B; 

   return (*this);
 }


RkRMatrix& RkRMatrix::setA (const RMatrix& aRMatrix)
 {
    LongInt k  = aRMatrix.numberOfColumns();
    LongInt rk = rank();
 
    if(rk!=k)
     {
        throw SimpleException(IString("Warning In RkRMatrix.setA, Bad Rank !!!"));
     }
    else
     {
        attr_A = aRMatrix;
     }

   return (*this);
 }


RkRMatrix& RkRMatrix::setB (const RMatrix& aRMatrix)
 {
    LongInt k  = aRMatrix.numberOfColumns();
    LongInt rk = rank();
 
    if(rk!=k)
     {
        throw SimpleException(IString("Warning In RkRMatrix.setB, Bad Rank !!!"));
     }
    else
     {
        attr_B = aRMatrix;
     }

   return (*this);
 }


RkRMatrix& RkRMatrix::setAB(const RMatrix& aRMatrix1, const RMatrix& aRMatrix2)
 {
    LongInt k1 = aRMatrix1.numberOfColumns();
    LongInt k2 = aRMatrix2.numberOfColumns();

    if(k1!=k2)
     {
        throw SimpleException(IString("Warning In RkRMatrix.setAB, A.numberOfColumns() != BnumberOfRows() !!!"));
     }
    else if(k1>rank())
     {
        throw SimpleException(IString("Warning In RkRMatrix.setAB, A.numberOfColumns()> rank() !!!"));
     }
    else
     {
        attr_A = aRMatrix1;
        attr_B = aRMatrix2;        

        setRank(k1);
     }

   return (*this);
 }


RkRMatrix& RkRMatrix::turncation2(const LongInt newRank)
 { 
    //! \todo  Hier k�nnte man die Multiplikation von qrA.Q*svd.U und analog f�r B mit ZUNMQR probieren.
    const LongInt k = attr_A.numberOfColumns();  
    //! QR Zerlegung von A  

    LongRealQR qrA(attr_A, k);

    //! QR Zerlegung von B
    LongRealQR qrB(attr_B, k);

    RMatrix R(qrA.attr_R*((RMatrix&) qrB.attr_R.setTranspose()));

    LongRealSvd svd(R, newRank);

    //! neues A setzen
    svd.multiplyUSigma();
    qrA.attr_Q *= svd.attr_U;    
    attr_A      = qrA.attr_Q;

    //! neues B setzen
    qrB.attr_Q *= svd.attr_V;
    attr_B      = qrB.attr_Q;
  
    setRank(newRank);

   return (*this);
 }


RkRMatrix& RkRMatrix::operator += (const RkRMatrix& aRkMatrix)
 {
    const LongInt rank1   = rank();
    const LongInt rank2   = aRkMatrix.rank();
    const LongInt newRank = MAX(rank1, rank2);   

    //! A und B aufweiten  
    attr_A.attachColumns(aRkMatrix.attr_A);
    attr_B.attachColumns(aRkMatrix.attr_B);

    //! Die formale Summe (Rang = k1+k2) auf newRank=max{k1,k2} k�rzen.
    turncation2(newRank);

   return (*this);
 }


RkRMatrix& RkRMatrix::operator += (const RMatrix& aMatrix)
 {
    RMatrix Temp(aMatrix);

    Temp.update(std::real(COMPLEX_UNIT), attr_A, (RMatrix&) attr_B.setTranspose(), std::real(COMPLEX_UNIT));

    LongRealSvd svd(Temp, rank());
 
    svd.multiplyUSigma();
    attr_A = svd.attr_U;

    attr_B = svd.attr_V;

   return (*this);
 }


RkRMatrix& RkRMatrix::addMatrix (const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex)
 {
    const AMatrixType theMatrixA = A.type();

    LongInt m  = numberOfRows();
    LongInt n  = numberOfColumns();
    LongInt m1 = A.numberOfRows();
    LongInt n1 = A.numberOfColumns();

    if(theMatrixA==isRkRMatrix)
     {
        const RkRMatrix& R = (const RkRMatrix&)A;
        LongInt k = R.rank();
        RMatrix At(m, k);
        RMatrix Bt(n, k); 

        int nx = 1;
            
        if(m1<=m && n1<=n)
         {
            for(LongInt i=0; i<k; i++)
             {
	       // dcopy (&MIN(m1, m), (LongReal*) &(R.attr_A(0, i)), &nx, (LongReal*) &(At(rowIndex, i)), &nx);
               // dcopy (&MIN(n1, n), (LongReal*) &(R.attr_B(0, i)), &nx, (LongReal*) &(Bt(colIndex, i)), &nx);
               TensorCalculus::Blas<double>::copy (MIN(m1, m), &(R.attr_A(0, i)), nx, &(At(rowIndex, i)), nx);
               TensorCalculus::Blas<double>::copy (MIN(n1, n), &(R.attr_B(0, i)), nx, &(Bt(colIndex, i)), nx);
             }
         }

        else if(m<m1 && n<n1)
         {
            for(LongInt i=0; i<k; i++)
             {
	       // dcopy (&MIN(m1, m), (LongReal*) &(R.attr_A(rowIndex, i)), &nx, (LongReal*) &(At(0, i)), &nx);
               // dcopy (&MIN(n1, n), (LongReal*) &(R.attr_B(colIndex, i)), &nx, (LongReal*) &(Bt(0, i)), &nx);
               TensorCalculus::Blas<double>::copy (MIN(m1, m), &(R.attr_A(rowIndex, i)), nx, &(At(0, i)), nx);
               TensorCalculus::Blas<double>::copy (MIN(n1, n), &(R.attr_B(colIndex, i)), nx, &(Bt(0, i)), nx);
              }
         }
        else
         {
            throw SimpleException(IString("Warning In RkRMatrix::addMatrix (const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex), Case Not Implemented or Bad!!!"));
         }
             
        //(*this) += T;
        //! A und B aufweiten  
        attr_A.attachColumns(At);
        attr_B.attachColumns(Bt);

        //! Die formale Summe k�rzen.
        turncation2(MAX(rank(),k));
     }
    else if(theMatrixA==isRMatrix)
     {
        const RMatrix& C = (const RMatrix&)A;
        int nx = 1;

        RMatrix T(m, n);

        if(m1<m && n1<n)
         {
            for(LongInt i=0; i<n1; i++)
             {
	       // dcopy (&m, (LongReal*) &(C(0, i)), &nx, (LongReal*) &(T(rowIndex, colIndex+i)), &nx);
               TensorCalculus::Blas<double>::copy (m, &(C(0, i)), nx, &(T(rowIndex, colIndex+i)), nx);
             }
         }
        else if(m<m1 && n<n1)
         {        
            for(LongInt i=0; i<n; i++)
             {
	       // dcopy (&m, (LongReal*) &(C(rowIndex, colIndex+i)), &nx, (LongReal*) &(T(0, i)), &nx);
               TensorCalculus::Blas<double>::copy (m, &(C(rowIndex, colIndex+i)), nx, &(T(0, i)), nx);
              }
         }
        else if(m==m1 && n==n1)
         {        
            for(LongInt i=0; i<n; i++)
             {
	       // dcopy (&m, (LongReal*) &(C(0, i)), &nx, (LongReal*) &(T(0, i)), &nx);
               TensorCalculus::Blas<double>::copy (m, &(C(0, i)), nx, &(T(0, i)), nx);
             }
         }
        else
         {
            throw SimpleException(IString("Warning In RkRMatrix::addMatrix (const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex), Case Not Implemented or Bad!!!"));
         }

        //(*this) += T;
        T.update(attr_A, (attr_B.setTranspose()));

        LongRealSvd svd(T, rank());

        svd.multiplyUSigma();
        attr_A = svd.attr_U;
        attr_B = svd.attr_V;
     }
    else if(theMatrixA==isHMatrixInterface)
     {
        RkRMatrix R(m1, n1, rank());

        R.isConversionOf(A);

        if( rowIndex==0 && colIndex==0 && m==m1 && n==n1)
         {
            (*this) += R;
         }
        else
         {    
            (*this).addMatrix(R, rowIndex, colIndex);
         }
     }
    else
     {
        throw SimpleException(IString("Warning In RkRMatrix::addMatrix (const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex), Type Not Implemented !!!"));
     }    

   return (*this);
 }


RkRMatrix& RkRMatrix::addScalMatrix (const LongReal& z, const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex)
 {
    const AMatrixType theMatrixA = A.type();

    LongInt m  = numberOfRows();
    LongInt n  = numberOfColumns();
    LongInt m1 = A.numberOfRows();
    LongInt n1 = A.numberOfColumns();

    if(theMatrixA==isRkRMatrix)
     {
        const RkRMatrix& R = (const RkRMatrix&)A;
        LongInt k = R.rank();
        RMatrix At(m, k);
        RMatrix Bt(n, k); 

        int nx = 1;
            
        if(m1<=m && n1<=n)
         {
            for(LongInt i=0; i<k; i++)
             {
	       // dcopy (&MIN(m1, m), (LongReal*) &(R.attr_A(0, i)), &nx, (LongReal*) &(At(rowIndex, i)), &nx);
               // dcopy (&MIN(n1, n), (LongReal*) &(R.attr_B(0, i)), &nx, (LongReal*) &(Bt(colIndex, i)), &nx);
               TensorCalculus::Blas<double>::copy (MIN(m1, m), &(R.attr_A(0, i)), nx, &(At(rowIndex, i)), nx);
               TensorCalculus::Blas<double>::copy (MIN(n1, n), &(R.attr_B(0, i)), nx, &(Bt(colIndex, i)), nx);
             }
         }

        else if(m<m1 && n<n1)
         {
            for(LongInt i=0; i<k; i++)
             {
	       // dcopy (&MIN(m1, m), (LongReal*) &(R.attr_A(rowIndex, i)), &nx, (LongReal*) &(At(0, i)), &nx);
               // dcopy (&MIN(n1, n), (LongReal*) &(R.attr_B(colIndex, i)), &nx, (LongReal*) &(Bt(0, i)), &nx);
	       TensorCalculus::Blas<double>::copy (MIN(m1, m), &(R.attr_A(rowIndex, i)), nx, &(At(0, i)), nx);
	       TensorCalculus::Blas<double>::copy (MIN(n1, n), &(R.attr_B(colIndex, i)), nx, &(Bt(0, i)), nx);
               }
         }
        else
         {
            throw SimpleException(IString("Warning In RkRMatrix::addMatrix (const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex), Case Not Implemented or Bad!!!"));
         }
             
        At *= z;  
        //! A und B aufweiten
        attr_A.attachColumns(At);
        attr_B.attachColumns(Bt);

        //! Die formale Summe k�rzen.
        turncation2(MAX(rank(),k));
     }
    else if(theMatrixA==isRMatrix)
     {
        const RMatrix& C = (const RMatrix&)A;
        int nx = 1;

        RMatrix T(m, n);

        if(m1<m && n1<n)
         {
            for(LongInt i=0; i<n1; i++)
             {
	       // dcopy (&m, (LongReal*) &(C(0, i)), &nx, (LongReal*) &(T(rowIndex, colIndex+i)), &nx);
	       TensorCalculus::Blas<double>::copy (m, &(C(0, i)), nx, &(T(rowIndex, colIndex+i)), nx);
             }
         }
        else if(m<m1 && n<n1)
         {        
            for(LongInt i=0; i<n; i++)
             {
	       // dcopy (&m, (LongReal*) &(C(rowIndex, colIndex+i)), &nx, (LongReal*) &(T(0, i)), &nx);
               TensorCalculus::Blas<double>::copy (m, &(C(rowIndex, colIndex+i)), nx, &(T(0, i)), nx);
              }
         }
        else if(m==m1 && n==n1)
         {        
            for(LongInt i=0; i<n; i++)
             {
	       // dcopy (&m, (LongReal*) &(C(0, i)), &nx, (LongReal*) &(T(0, i)), &nx);
               TensorCalculus::Blas<double>::copy (m, &(C(0, i)), nx, &(T(0, i)), nx);
             }
         }
        else
         {
            throw SimpleException(IString("Warning In RkRMatrix::addMatrix (const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex), Case Not Implemented or Bad!!!"));
         }

        //(*this) += z*T;
        T *= z;
        T.update(attr_A, (attr_B.setTranspose()));

        LongRealSvd svd(T, rank());

        svd.multiplyUSigma();
        attr_A = svd.attr_U;
        attr_B = svd.attr_V;
     }
    else if(theMatrixA==isHMatrixInterface)
     {
        RkRMatrix R(m1, n1, rank());

        R.isConversionOf(A);
        R *= z;

        if( rowIndex==0 && colIndex==0 && m==m1 && n==n1)
         {
            (*this) += R;
         }
        else
         {    
            (*this).addMatrix(R, rowIndex, colIndex);
         }
     }
    else
     {
        throw SimpleException(IString("Warning In RkRMatrix::addScalMatrix (const LongReal& z, const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex), Type Not Implemented !!!"));
     }    

   return (*this);
 }


RkRMatrix& RkRMatrix::addIdentityScaledBy(const LongReal& z)
 {
    const LongInt m   = numberOfRows();
    const LongInt n   = numberOfColumns();
    const LongInt dim = MIN(m,n);

    RMatrix zId(m, n);
 
    for(LongInt i=0; i<dim; i++)
     {
        zId(i,i) = z;
     }

    addMatrix(zId);

   return(*this);
 }


RkRMatrix& RkRMatrix::isConversionOf(const AMatrix& A)
 {
    const AMatrixType theMatrixA = A.type();

    const LongInt m    = A.numberOfRows();
    const LongInt n    = A.numberOfColumns();
    const LongInt rank = RkRMatrix::rank();

    setDimension(m, n, rank);

    if(theMatrixA==isRMatrix)
     {
        RMatrix T((RMatrix&) A);
        LongRealSvd svd(T, rank);

        svd.multiplyUSigma();
        attr_A = svd.attr_U;
        attr_B = svd.attr_V;
     }
    else if(theMatrixA==isRkRMatrix)
     {        
        (*this) = (RkRMatrix&) A;
        turncation2(rank);   
     }
    else if(theMatrixA==isHMatrixInterface)
     {
        const HMatrixInterface& Ha = (const HMatrixInterface&)A;

        const LongInt blockM = Ha.numberOfRowBlocks();
        const LongInt blockN = Ha.numberOfColumnBlocks();

        HMatrixInterface H(m, n, blockM, blockN);

        LongInt maxRank = 0;
        for(LongInt i=0; i<blockM; i++)
         {
            for(LongInt j=0; j<blockN; j++)
             { 
                AMatrix* Ha_ij = Ha(i,j);
                
                H(i,j) = new RkRMatrix(Ha_ij->numberOfRows(), Ha_ij->numberOfColumns(), rank);
                maxRank += MIN(rank, Ha_ij->numberOfColumns());
                H(i,j)->isConversionOf(*Ha_ij);
             }            
         }                
        isPreparedConversionOf(H, rank, maxRank);
     }
    else
     {
        throw SimpleException(IString("Warning In RkRMatrix::isConversionOf(const AMatrix& A), Type Not Implemented !!!"));
     }

   return (*this);
 }


//! (A1,B1)<-(A1,B1)*(A2,B2)
RkRMatrix& RkRMatrix::operator *= (RkRMatrix& aRkMatrix)
 {
    const LongInt k1 = rank();
    const LongInt k2 = aRkMatrix.rank();

    ((RMatrix&) attr_B.setTranspose()) *= aRkMatrix.attr_A;      //! B1 <- B1*A2 (abgespeichert in B1)

    if(k2<=k1)
     {
        attr_A *= (attr_B);         //! A1 <- A1*(B1*A2)
        attr_B  = aRkMatrix.attr_B;
        setRank(k1);
     }
    else
     {        
        ((RMatrix&) attr_B.setTranspose()).leftMultiplied(aRkMatrix.attr_B);
        setRank(k2);
     }
   return (*this);
 }


RkRMatrix& RkRMatrix::operator *= (RMatrix& aMatrix)
 {
    attr_B.leftMultiplied((RMatrix&) (aMatrix.setTranspose()));
   return (*this);
 }


RkRMatrix& RkRMatrix::operator *= (const LongReal& alpha)
 {
    attr_A *= alpha; 

   return (*this);
 }


RkRMatrix& RkRMatrix::operator /= (const LongReal& alpha)
 {
    attr_A *= 1.0/alpha; 

   return (*this);
 }


RVector RkRMatrix::operator * (const RVector& x)
 {
    RVector value(numberOfColumns()), temp(rank());
 
    temp  = ((RMatrix&) attr_B.setTranspose())*x;
    value = attr_A*temp;
    
   return value;
 }


bool RkRMatrix::operator () (const LongInt& i, RVector& w, const LongInt& j, const RVector& v)
 {
    bool value = true;

    RVector temp(rank());

    (attr_B.setTranspose())(0, temp, j, v);
    attr_A(i, w, 0, temp);  
   return value;
 }


bool RkRMatrix::partialEvaluateTransposeAt (const LongInt& i, RVector& w, const LongInt& j, const RVector& v)
 {    
    RVector temp(rank());

    attr_A.partialEvaluateTransposeAt(0, temp, j, v);

    attr_B(i, w, 0, temp);      

   return true;
 }


bool RkRMatrix::partialevaluateTransposeAt (const LongInt& i, RVector& w, const LongInt& j, const RVector& v)
 {    
    RVector temp(rank());

    attr_A.partialEvaluateTransposeAt(0, temp, j, v);

    attr_B(i, w, 0, temp);      

   return true;
 }


RkRMatrix operator + (const RkRMatrix& aRk, const RkRMatrix& aRkMatrix)
 {
    RkRMatrix value(aRk);

    value += aRkMatrix;

   return value;
 }


RkRMatrix operator * (const RkRMatrix& aRk, RkRMatrix& aRkMatrix)
 {
    RkRMatrix value(aRk);

    value *= aRkMatrix;

   return value;
 }


RkRMatrix operator * (const RkRMatrix& aRk, RMatrix&   aRMatrix)
 {
    RkRMatrix value(aRk);

    value *= aRMatrix;

   return value;
 }


RkRMatrix operator * (const RMatrix& aRMatrix, RkRMatrix& aRk)
 {
    RkRMatrix value(aRk);

    value.attr_A.leftMultiplied((RMatrix&)aRMatrix);

   return value;
 }

RkRMatrix operator * (const RkRMatrix& aRk, const LongReal& alpha) 
 {
    RkRMatrix value(aRk);

    value *= alpha;

   return value;
 }


RkRMatrix operator / (const RkRMatrix& aRk, const LongReal& alpha) 
 {
    RkRMatrix value(aRk);

    value /= alpha;

   return value;
 }


RMatrix RkRMatrix::convert2RMatrix() 
 {
    RMatrix Z(attr_A*(attr_B.setTranspose()));
   return Z; 
 }


AMatrix&  RkRMatrix::setNullMatrix ()
 {
    attr_A.setNullMatrix();
    attr_B.setNullMatrix();

   return (*this);
 }


RkRMatrix& RkRMatrix::update(const AMatrix& A, const AMatrix& B)
 {
    const AMatrixType theMatrixA = A.type();
    const AMatrixType theMatrixB = B.type();

    if(theMatrixA==isRMatrix || theMatrixB==isRMatrix)
     {
        const LongInt m = A.numberOfRows();
        const LongInt n = B.numberOfColumns();

        RMatrix C(m,n);
        C.update(A,B);

        C += (convert2RMatrix());

        LongRealSvd svd(C, rank());
        svd.multiplyUSigma();

        attr_A = svd.attr_U;
        attr_B = svd.attr_V;
     }
    else if (theMatrixA==isRkRMatrix || theMatrixB==isRkRMatrix)
     {
        if( theMatrixA==isRkRMatrix && theMatrixB!=isRkRMatrix )
         {
            const LongInt m = A.numberOfRows();
            const LongInt n = B.numberOfColumns();

            RkRMatrix& RA = (RkRMatrix&) A;

            const LongInt k = RA.rank();

            RkRMatrix R(m, n, k);
            
            R.setA(RA.attr_A);

            for(LongInt j=0; j<k; j++)
             {
                RVector helpVp_v(n, &RA.attr_B(0,j)); 
                RVector helpVp_w(n, &R.attr_B(0,j));

                B.evaluateTransposeAt(helpVp_v, helpVp_w);
             }

            (*this) += R;            
         }
        else if(theMatrixA!=isRkRMatrix && theMatrixB==isRkRMatrix)
         {
            const LongInt m = A.numberOfRows();
            const LongInt n = B.numberOfColumns();

            RkRMatrix& RB = (RkRMatrix&) B; 
            RkRMatrix R(m, n, RB.rank());
            
            R.attr_A.update(A, RB.attr_A);
            R.setB(RB.attr_B);

            (*this) += R;
         }
        else
         {
            RkRMatrix R((RkRMatrix&)A);
            R *= (RkRMatrix&)B;

            (*this) += R;
         }
     }
    else if (theMatrixA==isHMatrixInterface && theMatrixB==isHMatrixInterface)
     {
        const LongInt m = A.numberOfRows();
        const LongInt n = B.numberOfColumns();

        const HMatrixInterface& H1 = (HMatrixInterface&) A;
        const HMatrixInterface& H2 = (HMatrixInterface&) B;

        const LongInt blockM = H1.numberOfRowBlocks();
        const LongInt blockK = H1.numberOfColumnBlocks();
        const LongInt blockN = H2.numberOfColumnBlocks();

        HMatrixInterface H(m, n, blockM, blockN);

        const LongInt rank = RkRMatrix::rank();

        LongInt m1 = 0;
        LongInt n1 = 0; 

        LongInt maxRank = 0;

        for(LongInt i=0; i<blockM; i++)
         {
            for(LongInt j=0; j<blockN; j++)
             { 
                m1 = H1(i,0)->numberOfRows();
                n1 = H2(0,j)->numberOfColumns();

                maxRank += MIN(n1, rank);

                H(i,j) = new RkRMatrix(m1, n1, rank);
                for(LongInt k=0; k<blockK; k++)
                 {
                    H(i,j)->update(H1(i,k), H2(k,j));
                 }
             }            
         }
        RkRMatrix R(m, n, maxRank);
        R.isPreparedConversionOf(H, rank, maxRank);
      (*this) += R;              
     }
    else 
     {
        throw SimpleException(IString("Warning In RkRMatrix::update(const AMatrix& A, const AMatrix& B), Matrix Type Not Considered !!!"));
     }

   return (*this);
 }


RkRMatrix& RkRMatrix::isMultiplicationOf(const AMatrix& A, const AMatrix& B)
 {
    const AMatrixType theMatrixA = A.type();
    const AMatrixType theMatrixB = B.type();

    if(theMatrixA==isRMatrix || theMatrixB==isRMatrix)
     {
        const LongInt m = A.numberOfRows();
        const LongInt n = B.numberOfColumns();

        RMatrix C(m,n);

        C.update(A,B);

        LongRealSvd svd(C, rank());
        svd.multiplyUSigma();

        attr_A = svd.attr_U;
        attr_B = svd.attr_V;
     }
    else if (theMatrixA==isRkRMatrix || theMatrixB==isRkRMatrix)
     {
        if( theMatrixA==isRkRMatrix && theMatrixB!=isRkRMatrix )
         {
            const LongInt m = A.numberOfRows();
            const LongInt n = B.numberOfColumns();

            RkRMatrix& RA = (RkRMatrix&) A;

            const LongInt k = RA.rank();

            setDimension(m, n, k);
            
            setA(RA.attr_A);

            for(LongInt j=0; j<k; j++)
             {
                RVector helpVp_v(n, &RA.attr_B(0,j)); 
                RVector helpVp_w(n, &attr_B(0,j));

                B.evaluateTransposeAt(helpVp_v, helpVp_w);
             }                        
         }
        else if(theMatrixA!=isRkRMatrix && theMatrixB==isRkRMatrix)
         {
            const LongInt m = A.numberOfRows();
            const LongInt n = B.numberOfColumns();

            RkRMatrix& RB = (RkRMatrix&) B; 
            setDimension(m, n, RB.rank());
            
            attr_A.isMultiplicationOf(&A, &RB.attr_A);
            setB(RB.attr_B);            
         }
        else
         {
            (*this)  = (RkRMatrix&)A;            
            (*this) *= (RkRMatrix&)B;
         }
     }
    else if (theMatrixA==isHMatrixInterface && theMatrixB==isHMatrixInterface)
     {
        const LongInt m = A.numberOfRows();
        const LongInt n = B.numberOfColumns();

        const HMatrixInterface& H1 = (HMatrixInterface&) A;
        const HMatrixInterface& H2 = (HMatrixInterface&) B;

        const LongInt blockM = H1.numberOfRowBlocks();
        const LongInt blockK = H1.numberOfColumnBlocks();
        const LongInt blockN = H2.numberOfColumnBlocks();

        HMatrixInterface H(m, n, blockM, blockN);

        const LongInt rank = RkRMatrix::rank();

        LongInt m1 = 0;
        LongInt n1 = 0; 

        LongInt maxRank = 0;

        for(LongInt i=0; i<blockM; i++)
         {
            for(LongInt j=0; j<blockN; j++)
             { 
                m1 = H1(i,0)->numberOfRows();
                n1 = H2(0,j)->numberOfColumns();

                maxRank += MIN(n1, rank);

                H(i,j) = new RkRMatrix(m1, n1, rank);
                for(LongInt k=0; k<blockK; k++)
                 {
                    H(i,j)->update(H1(i,k), H2(k,j));
                 }
             }            
         }

        setDimension(m, n, maxRank);

        isPreparedConversionOf(H, rank, maxRank);
     }
    else 
     {
        throw SimpleException(IString("Warning In RkRMatrix::isMultiplicationOf(const AMatrix& A, const AMatrix& B), Matrix Type Not Considered !!!"));
     }
        
   return (*this);
 }


bool RkRMatrix::isPreparedConversionOf (const HMatrixInterface& H, const LongInt& newRank, const LongInt& maxRank)
 {
    const LongInt m = H.numberOfRows();
    const LongInt n = H.numberOfColumns();    

    setDimension(m, n, maxRank);

    const LongInt blockM = H.numberOfRowBlocks();
    const LongInt blockN = H.numberOfColumnBlocks();

    LongInt rowBlockIndex = 0;
    LongInt colBlockIndex = 0;
    LongInt rank          = 0;
    LongInt eins          = 1;
    LongInt indexRank     = 0;

    for(int i=0; i<blockM; i++)
     {        
        for(int j=0; j<blockN; j++)
         {
            // Achtung, eigentlich sollte man hier den Type pr�fen,
            // mach ich jetzt nicht, da  isPreparedConversionOf private ist !!!
            RkRMatrix& R = (RkRMatrix&) *H(i,j);
            rank = R.rank();            
            for(int k=0; k<rank; k++)
             {
	       // dcopy (&R.attr_numberOfRows, (LongReal*) &R.attr_A(0, k), &eins, 
               //                               (LongReal*) &attr_A(rowBlockIndex, k+indexRank), &eins
               //         );
               //  dcopy (&R.attr_numberOfColumns, (LongReal*) &R.attr_B(0, k), &eins, 
               //                                  (LongReal*) &attr_B(colBlockIndex, k+indexRank), &eins
               //         ); 
               TensorCalculus::Blas<double>::copy (R.attr_numberOfRows, &R.attr_A(0, k), eins, 
                                             &attr_A(rowBlockIndex, k+indexRank), eins
                       );
               TensorCalculus::Blas<double>::copy (R.attr_numberOfColumns, &R.attr_B(0, k), eins, 
                                                &attr_B(colBlockIndex, k+indexRank), eins
                       ); 
             }
            colBlockIndex += R.numberOfColumns();
            indexRank     += rank;
         }
        colBlockIndex = 0;
        // Achtung, eigentlich sollte man hier den Type pr�fen,
        // mach ich jetzt nicht, da  isPreparedConversionOf private ist !!!
        RkRMatrix& R1 = (RkRMatrix&) *H(i,0);
        rowBlockIndex += R1.numberOfRows();
     }
    turncation2(newRank);   

   return true;
 }


RkRMatrix& RkRMatrix::setDimension(const LongInt& m, const LongInt& n, const LongInt& k)
 {
    attr_A.setDimension(m, k);
    attr_B.setDimension(n, k);

    (*this)
     .setRank(k)
     .setNumberOfRows(m)
     .setNumberOfColumns(n)
    ;

   return (*this);
 }


ostream& operator << (ostream& s, const RkRMatrix& R)
 {
    s << R.attr_A << endl << R.attr_B << endl;
   return s;
 }


istream& operator >> (istream& s, RkRMatrix& R)
 {
    s >> R.attr_A >> R.attr_B;

    R.setNumberOfRows(R.attr_A.numberOfRows());
    R.setRank(R.attr_A.numberOfColumns());
    R.setNumberOfColumns(R.attr_B.numberOfRows());
   return s;
 }


LongReal RkRMatrix::frobeniusNorm () const
 {
    LongInt m = numberOfRows();
    LongInt n = numberOfColumns();
    LongInt k = rank();

    LongReal norm = 0.0, temp = 0.0;

    for(LongInt i=0; i<m; i++)
     {
        RVector a(k, &attr_A(i,0));

        for(LongInt j=0; j<n; j++)
         {            
            RVector b(k, &attr_B(j,0));
            
            temp =  innerProduct(a,b);
            norm += temp*temp;
         }
     }

   return sqrt(norm);
 }
