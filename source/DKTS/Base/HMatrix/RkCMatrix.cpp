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

// RkCMatrix.cpp: Implementierung der Klasse RkCMatrix.
//
//////////////////////////////////////////////////////////////////////

#include "RkCMatrix.hpp"
#include "HMatrixInterface.hpp"
#include "BlasInterface.hpp"
// #include "LapackInterface2.hpp"

RkCMatrix::RkCMatrix(const LongInt m, const LongInt n, const LongInt k)
:ACMatrix(m, n), attr_A(m, MIN(n, k)), attr_B(n, MIN(n, k))
 {
    (*this)
     .setRank(MIN(n, k))
     .setType(isRkCMatrix)
    ; 
 }


RkCMatrix::RkCMatrix(const RkCMatrix& aRkMatrix)
:ACMatrix(aRkMatrix.numberOfRows(), aRkMatrix.numberOfColumns()), 
 attr_A(aRkMatrix.A()), attr_B(aRkMatrix.B())
 {
    (*this)
     .setRank(aRkMatrix.rank())
     .setType(isRkCMatrix)
    ;
 }


RkCMatrix::RkCMatrix (const CMatrix& Z, const LongInt k)
:ACMatrix(Z.numberOfRows(), Z.numberOfColumns()), 
 attr_A(Z.numberOfRows(), MIN(Z.numberOfColumns(), k)), attr_B(Z.numberOfColumns(), MIN(Z.numberOfColumns(), k))
 {
    (*this)
     .setRank(MIN(Z.numberOfColumns(), k))
     .setType(isRkCMatrix)
    ;

    CMatrix Temp(Z);

    ComplexSvd svd(Temp, rank());
 
    svd.multiplyUSigma();
    attr_A = svd.attr_U;
    attr_B = svd.attr_V;

 }

RkCMatrix::~RkCMatrix()
 {
 }


RkCMatrix& RkCMatrix::operator = (const RkCMatrix& R)
 {
    (*this)
     .setRank(MIN(R.numberOfColumns(), R.rank()))
     .setNumberOfRows(R.numberOfRows())
     .setNumberOfColumns(R.numberOfColumns())
     .setType(isRkCMatrix)
    ; 

    attr_A = R.attr_A;
    attr_B = R.attr_B; 

   return (*this);
 }


RkCMatrix& RkCMatrix::setA (const CMatrix& aCMatrix)
 {
    LongInt k  = aCMatrix.numberOfColumns();
    LongInt rk = rank();
 
    if(rk!=k)
     {
        throw SimpleException(IString("Warning In RkCMatrix.setA, Bad Rank !!!"));
     }
    else
     {
        attr_A = aCMatrix;
     }

   return (*this);
 }


RkCMatrix& RkCMatrix::setB (const CMatrix& aCMatrix)
 {
    LongInt k  = aCMatrix.numberOfColumns();
    LongInt rk = rank();
 
    if(rk!=k)
     {
        throw SimpleException(IString("Warning In RkCMatrix.setB, Bad Rank !!!"));
     }
    else
     {
        attr_B = aCMatrix;
     }

   return (*this);
 }


RkCMatrix& RkCMatrix::setAB(const CMatrix& aCMatrix1, const CMatrix& aCMatrix2)
 {
    LongInt k1 = aCMatrix1.numberOfColumns();
    LongInt k2 = aCMatrix2.numberOfColumns();

    if(k1!=k2)
     {
        throw SimpleException(IString("Warning In RkCMatrix.setAB, A.numberOfColumns() != BnumberOfRows() !!!"));
     }
    else if(k1>rank())
     {
        throw SimpleException(IString("Warning In RkCMatrix.setAB, A.numberOfColumns()> rank() !!!"));
     }
    else
     {
        attr_A = aCMatrix1;
        attr_B = aCMatrix2;        

        setRank(k1);
     }

   return (*this);
 }


RkCMatrix& RkCMatrix::turncation2(const LongInt newRank)
 { 
    //! \todo  Hier k�nnte man die Multiplikation von qrA.Q*svd.U und analog f�r B mit ZUNMQR probieren.
    const LongInt k = attr_A.numberOfColumns();  
    //! QR Zerlegung von A  
    ComplexQR qrA(attr_A, k);

    //! QR Zerlegung von B
    ComplexQR qrB(attr_B, k);

    CMatrix R(qrA.attr_R*(qrB.attr_R.setConjugateTranspose()));

    ComplexSvd svd(R, newRank);

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


RkCMatrix& RkCMatrix::operator += (const RkCMatrix& aRkMatrix)
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


RkCMatrix& RkCMatrix::operator += (const CMatrix& aMatrix)
 {
    CMatrix Temp(aMatrix);

    Temp.update(COMPLEX_UNIT, attr_A, attr_B.setConjugateTranspose(), COMPLEX_UNIT);

    ComplexSvd svd(Temp, rank());
 
    svd.multiplyUSigma();
    attr_A = svd.attr_U;

    attr_B = svd.attr_V;

   return (*this);
 }


RkCMatrix& RkCMatrix::addMatrix (const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex)
 {
    const AMatrixType theMatrixA = A.type();

    LongInt m  = numberOfRows();
    LongInt n  = numberOfColumns();
    LongInt m1 = A.numberOfRows();
    LongInt n1 = A.numberOfColumns();

    if(theMatrixA==isRkCMatrix)
     {
        const RkCMatrix& R = (const RkCMatrix&)A;
        LongInt k = R.rank();
        CMatrix At(m, k);
        CMatrix Bt(n, k); 

        int nx = 1;
            
        if(m1<=m && n1<=n)
         {
            for(LongInt i=0; i<k; i++)
             {
	       TensorCalculus::Blas< std::complex<double> >::copy (MIN(m1, m), &(R.attr_A(0, i)), nx, &(At(rowIndex, i)), nx);
	       TensorCalculus::Blas< std::complex<double> >::copy (MIN(n1, n), &(R.attr_B(0, i)), nx, &(Bt(colIndex, i)), nx);
                // zcopy (&MIN(m1, m), (_MKL_Complex16*) &(R.attr_A(0, i)), &nx, (_MKL_Complex16*) &(At(rowIndex, i)), &nx);
                // zcopy (&MIN(n1, n), (_MKL_Complex16*) &(R.attr_B(0, i)), &nx, (_MKL_Complex16*) &(Bt(colIndex, i)), &nx);
             }
         }

        else if(m<m1 && n<n1)
         {
            for(LongInt i=0; i<k; i++)
             {
	       TensorCalculus::Blas< std::complex<double> >::copy (MIN(m1, m), &(R.attr_A(rowIndex, i)), nx, &(At(0, i)), nx);
               TensorCalculus::Blas< std::complex<double> >::copy (MIN(n1, n), &(R.attr_B(colIndex, i)), nx, &(Bt(0, i)), nx);
		// zcopy (&MIN(m1, m), (_MKL_Complex16*) &(R.attr_A(rowIndex, i)), &nx, (_MKL_Complex16*) &(At(0, i)), &nx);
                // zcopy (&MIN(n1, n), (_MKL_Complex16*) &(R.attr_B(colIndex, i)), &nx, (_MKL_Complex16*) &(Bt(0, i)), &nx);
	     }
         }
        else
         {
            throw SimpleException(IString("Warning In RkCMatrix::addMatrix (const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex), Case Not Implemented or Bad!!!"));
         }
             
        //(*this) += T;
        //! A und B aufweiten  
        attr_A.attachColumns(At);
        attr_B.attachColumns(Bt);

        //! Die formale Summe k�rzen.
        turncation2(MAX(rank(),k));
     }
    else if(theMatrixA==isCMatrix)
     {
        const CMatrix& C = (const CMatrix&)A;
        int nx = 1;

        CMatrix T(m, n);

        if(m1<m && n1<n)
         {
            for(LongInt i=0; i<n1; i++)
             {
	       // zcopy (&m, (_MKL_Complex16*) &(C(0, i)), &nx, (_MKL_Complex16*) &(T(rowIndex, colIndex+i)), &nx);
               TensorCalculus::Blas< std::complex<double> >::copy(m, &(C(0, i)), nx, &(T(rowIndex, colIndex+i)), nx);
             }
         }
        else if(m<m1 && n<n1)
         {        
            for(LongInt i=0; i<n; i++)
             {
	       // zcopy (&m, (_MKL_Complex16*) &(C(rowIndex, colIndex+i)), &nx, (_MKL_Complex16*) &(T(0, i)), &nx);
               TensorCalculus::Blas< std::complex<double> >::copy (m, &(C(rowIndex, colIndex+i)), nx, &(T(0, i)), nx);
             }
         }
        else if(m==m1 && n==n1)
         {        
            for(LongInt i=0; i<n; i++)
             {
	       // zcopy (&m, (_MKL_Complex16*) &(C(0, i)), &nx, (_MKL_Complex16*) &(T(0, i)), &nx);
               TensorCalculus::Blas< std::complex<double> >::copy (m, &(C(0, i)), nx, &(T(0, i)), nx);
             }
         }
        else
         {
            throw SimpleException(IString("Warning In RkCMatrix::addMatrix (const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex), Case Not Implemented or Bad!!!"));
         }

        //(*this) += T;
        T.update(attr_A, (attr_B.setConjugateTranspose()));

        ComplexSvd svd(T, rank());

        svd.multiplyUSigma();
        attr_A = svd.attr_U;
        attr_B = svd.attr_V;
     }
    else if(theMatrixA==isHMatrixInterface)
     {
        RkCMatrix R(m1, n1, rank());

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
        throw SimpleException(IString("Warning In RkCMatrix::addMatrix (const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex), Type Not Implemented !!!"));
     }    

   return (*this);
 }


RkCMatrix& RkCMatrix::addScalMatrix (const std::complex<double>& z, const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex)
 {
    const AMatrixType theMatrixA = A.type();

    LongInt m  = numberOfRows();
    LongInt n  = numberOfColumns();
    LongInt m1 = A.numberOfRows();
    LongInt n1 = A.numberOfColumns();

    if(theMatrixA==isRkCMatrix)
     {
        const RkCMatrix& R = (const RkCMatrix&)A;
        LongInt k = R.rank();
        CMatrix At(m, k);
        CMatrix Bt(n, k); 

        int nx = 1;
            
        if(m1<=m && n1<=n)
         {
            for(LongInt i=0; i<k; i++)
             {
	       // zcopy (&MIN(m1, m), (_MKL_Complex16*) &(R.attr_A(0, i)), &nx, (_MKL_Complex16*) &(At(rowIndex, i)), &nx);
               // zcopy (&MIN(n1, n), (_MKL_Complex16*) &(R.attr_B(0, i)), &nx, (_MKL_Complex16*) &(Bt(colIndex, i)), &nx);
               TensorCalculus::Blas< std::complex<double> >::copy(MIN(m1, m), &(R.attr_A(0, i)), nx, &(At(rowIndex, i)), nx);
               TensorCalculus::Blas< std::complex<double> >::copy(MIN(n1, n), &(R.attr_B(0, i)), nx, &(Bt(colIndex, i)), nx);
             }
         }

        else if(m<m1 && n<n1)
         {
            for(LongInt i=0; i<k; i++)
             {
	       // zcopy (&MIN(m1, m), (_MKL_Complex16*) &(R.attr_A(rowIndex, i)), &nx, (_MKL_Complex16*) &(At(0, i)), &nx);
               // zcopy (&MIN(n1, n), (_MKL_Complex16*) &(R.attr_B(colIndex, i)), &nx, (_MKL_Complex16*) &(Bt(0, i)), &nx);
               TensorCalculus::Blas< std::complex<double> >::copy(MIN(m1, m), &(R.attr_A(rowIndex, i)), nx, &(At(0, i)), nx);
               TensorCalculus::Blas< std::complex<double> >::copy(MIN(n1, n), &(R.attr_B(colIndex, i)), nx, &(Bt(0, i)), nx);
             }
         }
        else
         {
            throw SimpleException(IString("Warning In RkCMatrix::addMatrix (const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex), Case Not Implemented or Bad!!!"));
         }
             
        At *= z;  
        //! A und B aufweiten
        attr_A.attachColumns(At);
        attr_B.attachColumns(Bt);

        //! Die formale Summe k�rzen.
        turncation2(MAX(rank(),k));
     }
    else if(theMatrixA==isCMatrix)
     {
        const CMatrix& C = (const CMatrix&)A;
        int nx = 1;

        CMatrix T(m, n);

        if(m1<m && n1<n)
         {
            for(LongInt i=0; i<n1; i++)
             {
	       // zcopy (&m, (_MKL_Complex16*) &(C(0, i)), &nx, (_MKL_Complex16*) &(T(rowIndex, colIndex+i)), &nx);
	       TensorCalculus::Blas< std::complex<double> >::copy (m, &(C(0, i)), nx, &(T(rowIndex, colIndex+i)), nx);
             }
         }
        else if(m<m1 && n<n1)
         {        
            for(LongInt i=0; i<n; i++)
             {
	       // zcopy (&m, (_MKL_Complex16*) &(C(rowIndex, colIndex+i)), &nx, (_MKL_Complex16*) &(T(0, i)), &nx);
               TensorCalculus::Blas< std::complex<double> >::copy (m, &(C(rowIndex, colIndex+i)), nx, &(T(0, i)), nx);
             }
         }
        else if(m==m1 && n==n1)
         {        
            for(LongInt i=0; i<n; i++)
             {
	       // zcopy (&m, (_MKL_Complex16*) &(C(0, i)), &nx, (_MKL_Complex16*) &(T(0, i)), &nx);
               TensorCalculus::Blas< std::complex<double> >::copy (m, &(C(0, i)), nx, &(T(0, i)), nx);
             }
         }
        else
         {
            throw SimpleException(IString("Warning In RkCMatrix::addMatrix (const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex), Case Not Implemented or Bad!!!"));
         }

        //(*this) += z*T;
        T *= z;
        T.update(attr_A, (attr_B.setConjugateTranspose()));

        ComplexSvd svd(T, rank());

        svd.multiplyUSigma();
        attr_A = svd.attr_U;
        attr_B = svd.attr_V;
     }
    else if(theMatrixA==isHMatrixInterface)
     {
        RkCMatrix R(m1, n1, rank());

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
        throw SimpleException(IString("Warning In RkCMatrix::addScalMatrix (const Complex& z, const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex), Type Not Implemented !!!"));
     }    

   return (*this);
 }


RkCMatrix& RkCMatrix::addIdentityScaledBy(const std::complex<double>& z)
 {
    const LongInt m   = numberOfRows();
    const LongInt n   = numberOfColumns();
    const LongInt dim = MIN(m,n);

    CMatrix zId(m, n);
 
    for(LongInt i=0; i<dim; i++)
     {
        zId(i,i) = z;
     }

    addMatrix(zId);

   return(*this);
 }


RkCMatrix& RkCMatrix::isConversionOf(const AMatrix& A)
 {
    const AMatrixType theMatrixA = A.type();

    const LongInt m    = A.numberOfRows();
    const LongInt n    = A.numberOfColumns();
    const LongInt rank = RkCMatrix::rank();

    setDimension(m, n, rank);

    if(theMatrixA==isCMatrix)
     {
        CMatrix T((CMatrix&) A);
        ComplexSvd svd(T, rank);

        svd.multiplyUSigma();
        attr_A = svd.attr_U;
        attr_B = svd.attr_V;
     }
    else if(theMatrixA==isRkCMatrix)
     {        
        (*this) = (RkCMatrix&) A;
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
                
                H(i,j) = new RkCMatrix(Ha_ij->numberOfRows(), Ha_ij->numberOfColumns(), rank);
                maxRank += MIN(rank, Ha_ij->numberOfColumns());
                H(i,j)->isConversionOf(*Ha_ij);
             }            
         }                
        isPreparedConversionOf(H, rank, maxRank);
     }
    else
     {
        throw SimpleException(IString("Warning In RkCMatrix::isConversionOf(const AMatrix& A), Type Not Implemented !!!"));
     }

   return (*this);
 }


//! (A1,B1)<-(A1,B1)*(A2,B2)
RkCMatrix& RkCMatrix::operator *= (RkCMatrix& aRkMatrix)
 {
    const LongInt k1 = rank();
    const LongInt k2 = aRkMatrix.rank();

    (attr_B.setConjugateTranspose()) *= aRkMatrix.attr_A;      //! B1 <- B1*A2 (abgespeichert in B1)

    if(k2<=k1)
     {
        attr_A *= (attr_B);         //! A1 <- A1*(B1*A2)
        attr_B  = aRkMatrix.attr_B;
        setRank(k1);
     }
    else
     {        
        (attr_B.setConjugateTranspose()).leftMultiplied(aRkMatrix.attr_B);
        setRank(k2);
     }
   return (*this);
 }


RkCMatrix& RkCMatrix::operator *= (CMatrix& aMatrix)
 {
    attr_B.leftMultiplied((aMatrix.setConjugateTranspose()));
   return (*this);
 }


RkCMatrix& RkCMatrix::operator *= (const std::complex<double>& alpha)
 {
    attr_A *= alpha; 

   return (*this);
 }


RkCMatrix& RkCMatrix::operator /= (const std::complex<double>& alpha)
 {
    attr_A *= 1.0/alpha; 

   return (*this);
 }


CVector RkCMatrix::operator * (const CVector& x)
 {
    CVector value(numberOfColumns()), temp(rank());
 
    temp  = (attr_B.setConjugateTranspose())*x;
    value = attr_A*temp;
    
   return value;
 }


bool RkCMatrix::operator () (const LongInt& i, CVector& w, const LongInt& j, const CVector& v)
 {
    bool value = true;

    CVector temp(rank());

    ((attr_B).setConjugateTranspose())(0, temp, j, v);

    attr_A(i, w, 0, temp);  

   return value;
 }


bool RkCMatrix::partialEvaluateTransposeAt (const LongInt& i, CVector& w, const LongInt& j, const CVector& v)
 {
    
    CVector temp(rank());

    attr_A.partialEvaluateTransposeAt(0, temp, j, v);

    attr_B.conjugated();
    attr_B(i, w, 0, temp);      
    attr_B.conjugated();

   return true;
 }


bool RkCMatrix::partialEvaluateConjugateTransposeAt (const LongInt& i, CVector& w, const LongInt& j, const CVector& v)
 {
    
    CVector temp(rank());

    attr_A.partialEvaluateConjugateTransposeAt(0, temp, j, v);

    attr_B(i, w, 0, temp);      

   return true;
 }


RkCMatrix operator + (const RkCMatrix& aRk, const RkCMatrix& aRkMatrix)
 {
    RkCMatrix value(aRk);

    value += aRkMatrix;

   return value;
 }


RkCMatrix operator * (const RkCMatrix& aRk, RkCMatrix& aRkMatrix)
 {
    RkCMatrix value(aRk);

    value *= aRkMatrix;

   return value;
 }


RkCMatrix operator * (const RkCMatrix& aRk, CMatrix&   aCMatrix)
 {
    RkCMatrix value(aRk);

    value *= aCMatrix;

   return value;
 }


RkCMatrix operator * (const CMatrix& aCMatrix, RkCMatrix& aRk)
 {
    RkCMatrix value(aRk);

    value.attr_A.leftMultiplied((CMatrix&)aCMatrix);

   return value;
 }

RkCMatrix operator * (const RkCMatrix& aRk, const std::complex<double>& alpha)
 {
    RkCMatrix value(aRk);

    value *= alpha;

   return value;
 }


RkCMatrix operator / (const RkCMatrix& aRk, const std::complex<double>& alpha)
 {
    RkCMatrix value(aRk);

    value /= alpha;

   return value;
 }


CMatrix RkCMatrix::convert2CMatrix() 
 {
    CMatrix Z(attr_A*(attr_B.setConjugateTranspose()));
   return Z; 
 }


AMatrix&  RkCMatrix::setNullMatrix ()
 {
    attr_A.setNullMatrix();
    attr_B.setNullMatrix();

   return (*this);
 }


RkCMatrix& RkCMatrix::update(const AMatrix& A, const AMatrix& B)
 {
    const AMatrixType theMatrixA = A.type();
    const AMatrixType theMatrixB = B.type();

    if(theMatrixA==isCMatrix || theMatrixB==isCMatrix)
     {
        const LongInt m = A.numberOfRows();
        const LongInt n = B.numberOfColumns();

        CMatrix C(m,n);
        C.update(A,B);

        C += (convert2CMatrix());

        ComplexSvd svd(C, rank());
        svd.multiplyUSigma();

        attr_A = svd.attr_U;
        attr_B = svd.attr_V;
     }
    else if (theMatrixA==isRkCMatrix || theMatrixB==isRkCMatrix)
     {
        if( theMatrixA==isRkCMatrix && theMatrixB!=isRkCMatrix )
         {
            const LongInt m = A.numberOfRows();
            const LongInt n = B.numberOfColumns();

            RkCMatrix& RA = (RkCMatrix&) A;

            const LongInt k = RA.rank();

            RkCMatrix R(m, n, k);
            
            R.setA(RA.attr_A);

            for(LongInt j=0; j<k; j++)
             {
                CVector helpVp_v(n, &RA.attr_B(0,j)); 
                CVector helpVp_w(n, &R.attr_B(0,j));
  
                B.evaluateConjugateTransposeAt(helpVp_v, helpVp_w);
             }

            (*this) += R;            
         }
        else if(theMatrixA!=isRkCMatrix && theMatrixB==isRkCMatrix)
         {
            const LongInt m = A.numberOfRows();
            const LongInt n = B.numberOfColumns();

            RkCMatrix& RB = (RkCMatrix&) B; 
            RkCMatrix R(m, n, RB.rank());
            
            R.attr_A.update(A, RB.attr_A);
            R.setB(RB.attr_B);

            (*this) += R;
         }
        else
         {
            RkCMatrix R((RkCMatrix&)A);
            R *= (RkCMatrix&)B;

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

        const LongInt rank = RkCMatrix::rank();

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

                H(i,j) = new RkCMatrix(m1, n1, rank);
                for(LongInt k=0; k<blockK; k++)
                 {
                    H(i,j)->update(H1(i,k), H2(k,j));
                 }
             }            
         }
        RkCMatrix R(m, n, maxRank);
        R.isPreparedConversionOf(H, rank, maxRank);
      (*this) += R;              
     }
    else 
     {
        throw SimpleException(IString("Warning In RkCMatrix::update(const AMatrix& A, const AMatrix& B), Matrix Type Not Considered !!!"));
     }

   return (*this);
 }


RkCMatrix& RkCMatrix::isMultiplicationOf(const AMatrix& A, const AMatrix& B)
 {
    const AMatrixType theMatrixA = A.type();
    const AMatrixType theMatrixB = B.type();

    if(theMatrixA==isCMatrix || theMatrixB==isCMatrix)
     {
        const LongInt m = A.numberOfRows();
        const LongInt n = B.numberOfColumns();

        CMatrix C(m,n);

        C.update(A,B);

        ComplexSvd svd(C, rank());
        svd.multiplyUSigma();

        attr_A = svd.attr_U;
        attr_B = svd.attr_V;
     }
    else if (theMatrixA==isRkCMatrix || theMatrixB==isRkCMatrix)
     {
        if( theMatrixA==isRkCMatrix && theMatrixB!=isRkCMatrix )
         {
            const LongInt m = A.numberOfRows();
            const LongInt n = B.numberOfColumns();

            RkCMatrix& RA = (RkCMatrix&) A;

            const LongInt k = RA.rank();

            setDimension(m, n, k);
            
            setA(RA.attr_A);

            for(LongInt j=0; j<k; j++)
             {
                CVector helpVp_v(n, &RA.attr_B(0,j)); 
                CVector helpVp_w(n, &attr_B(0,j));
  
                B.evaluateConjugateTransposeAt(helpVp_v, helpVp_w);
             }                        
         }
        else if(theMatrixA!=isRkCMatrix && theMatrixB==isRkCMatrix)
         {
            const LongInt m = A.numberOfRows();
            const LongInt n = B.numberOfColumns();

            RkCMatrix& RB = (RkCMatrix&) B; 
            setDimension(m, n, RB.rank());
            
            attr_A.isMultiplicationOf(&A, &RB.attr_A);
            setB(RB.attr_B);            
         }
        else
         {
            (*this)  = (RkCMatrix&)A;            
            (*this) *= (RkCMatrix&)B;
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

        const LongInt rank = RkCMatrix::rank();

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

                H(i,j) = new RkCMatrix(m1, n1, rank);
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
        throw SimpleException(IString("Warning In RkCMatrix::isMultiplicationOf(const AMatrix& A, const AMatrix& B), Matrix Type Not Considered !!!"));
     }
        
   return (*this);
 }


bool RkCMatrix::isPreparedConversionOf (const HMatrixInterface& H, const LongInt& newRank, const LongInt& maxRank)
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
            RkCMatrix& R = (RkCMatrix&) *H(i,j);
            rank = R.rank();            
            for(int k=0; k<rank; k++)
             {
	       // zcopy (&R.attr_numberOfRows, (_MKL_Complex16*) &R.attr_A(0, k), &eins, 
               //                               (_MKL_Complex16*) &attr_A(rowBlockIndex, k+indexRank), &eins
               //         );
               //  zcopy (&R.attr_numberOfColumns, (_MKL_Complex16*) &R.attr_B(0, k), &eins, 
               //                                  (_MKL_Complex16*) &attr_B(colBlockIndex, k+indexRank), &eins
               //         ); 
               TensorCalculus::Blas< std::complex<double> >::copy (R.attr_numberOfRows, &R.attr_A(0, k), eins,
                                             &attr_A(rowBlockIndex, k+indexRank), eins
                       );
               TensorCalculus::Blas< std::complex<double> >::copy (R.attr_numberOfColumns, &R.attr_B(0, k), eins,
                                                &attr_B(colBlockIndex, k+indexRank), eins
                       ); 
              }
            colBlockIndex += R.numberOfColumns();
            indexRank     += rank;
         }
        colBlockIndex = 0;
        // Achtung, eigentlich sollte man hier den Type pr�fen,
        // mach ich jetzt nicht, da  isPreparedConversionOf private ist !!!
        RkCMatrix& R1 = (RkCMatrix&) *H(i,0);
        rowBlockIndex += R1.numberOfRows();
     }
    turncation2(newRank);   

   return true;
 }


RkCMatrix& RkCMatrix::setDimension(const LongInt& m, const LongInt& n, const LongInt& k)
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


ostream& operator << (ostream& s, const RkCMatrix& R)
 {
    s << R.attr_A << endl << R.attr_B << endl;
   return s;
 }


istream& operator >> (istream& s, RkCMatrix& R)
 {
    s >> R.attr_A >> R.attr_B;

    R.setNumberOfRows(R.attr_A.numberOfRows());
    R.setRank(R.attr_A.numberOfColumns());
    R.setNumberOfColumns(R.attr_B.numberOfRows());
   return s;
 }


LongReal RkCMatrix::frobeniusNorm () const
 {
    LongInt m = numberOfRows();
    LongInt n = numberOfColumns();
    LongInt k = rank();

    LongReal norm = 0.0;

    for(LongInt i=0; i<m; i++)
     {
        CVector a(k, &attr_A(i,0));

        for(LongInt j=0; j<n; j++)
         {            
            CVector b(k, &attr_B(j,0));

            norm += std::norm(innerProduct(a,b));
         }
     }

   return sqrt(norm);
 }
