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

// HMatrixInterface.cpp: Implementierung der Klasse HMatrixInterface.
//
//////////////////////////////////////////////////////////////////////

#include "HMatrixInterface.hpp"

HMatrixInterface::HMatrixInterface ( const LongInt& m, const LongInt& n,
                                     const LongInt& block_m, const LongInt& block_n)
:AMatrix(m,n)
 {
    (*this)
     .setNumberOfRowBlocks(block_m)
     .setNumberOfColumnBlocks(block_n)
     .setComplexArithmetic()
     .setNumberOfRows(m)
     .setNumberOfColumns(n)
     .setType(isHMatrixInterface)
    ;

    attr_block = new AMatrixPointer [block_m*block_n];
 }
 

HMatrixInterface::~HMatrixInterface ()
 {
    delete [] attr_block;
 }


HMatrixInterface& HMatrixInterface::operator += (const AMatrix& A)
 {
    addMatrix(A, 0, 0);

   return (*this);
 }


HMatrixInterface& HMatrixInterface::addMatrix (const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex)
 {
    const LongInt blockM = numberOfRowBlocks();
    const LongInt blockN = numberOfColumnBlocks();

    const AMatrixType theMatrixA = A.type();
    
    if(theMatrixA == isHMatrixInterface)
     {
        const HMatrixInterface& H = (const HMatrixInterface&)A; 

        for(LongInt i=0; i<blockM; i++)
         {
            for(LongInt j=0; j<blockN; j++)
             {
                (*this)(i,j)->operator += (H(i,j));                
             }
         }
     }
    else 
     {
        LongInt indexM = rowIndex;
        LongInt indexN = colIndex;

        for(LongInt i=0; i<blockM; i++)
         {
            indexN  = colIndex;
            for(LongInt j=0; j<blockN; j++)
             {
                AMatrixPointer H = (*this)(i,j);

                H->addMatrix(A, indexM, indexN);
                indexN += H->numberOfColumns();
             }        
            indexM += ((*this)(i,0))->numberOfRows(); 
         }
     }

   return (*this);
 }


HMatrixInterface& HMatrixInterface::addScalMatrix (const LongReal& z, const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex)
 {
    const LongInt blockM = numberOfRowBlocks();
    const LongInt blockN = numberOfColumnBlocks();

    const AMatrixType theMatrixA = A.type();
    
    if(theMatrixA == isHMatrixInterface)
     {
        const HMatrixInterface& H = (const HMatrixInterface&)A; 
        const LongInt m = numberOfRowBlocks();
        const LongInt n = numberOfColumnBlocks();

        for(LongInt i=0; i<m; i++)
         {
            for(LongInt j=0; j<n; j++)
             {
                (*this)(i,j)->addScalMatrix(z, *H(i,j), 0, 0);                
             }
         }
     }
    else 
     {
        LongInt indexM = rowIndex;
        LongInt indexN = colIndex;

        for(LongInt i=0; i<blockM; i++)
         {
            indexN  = colIndex;
            for(LongInt j=0; j<blockN; j++)
             {
                AMatrixPointer H = (*this)(i,j);

                H->addScalMatrix(z, A, indexM, indexN);
                indexN += H->numberOfColumns();
             }        
            indexM += ((*this)(i,0))->numberOfRows(); 
         }
     }

   return (*this);
 }


HMatrixInterface& HMatrixInterface::addScalMatrix (const std::complex<double>& z, const AMatrix& A, const LongInt& rowIndex, const LongInt& colIndex)
 {
    const LongInt blockM = numberOfRowBlocks();
    const LongInt blockN = numberOfColumnBlocks();

    const AMatrixType theMatrixA = A.type();
    
    if(theMatrixA == isHMatrixInterface)
     {
        const HMatrixInterface& H = (const HMatrixInterface&)A; 
        const LongInt m = numberOfRowBlocks();
        const LongInt n = numberOfColumnBlocks();

        for(LongInt i=0; i<m; i++)
         {
            for(LongInt j=0; j<n; j++)
             {
                (*this)(i,j)->addScalMatrix(z, *H(i,j), 0, 0);                
             }
         }
     }
    else 
     {
        LongInt indexM = rowIndex;
        LongInt indexN = colIndex;

        for(LongInt i=0; i<blockM; i++)
         {
            indexN  = colIndex;
            for(LongInt j=0; j<blockN; j++)
             {
                AMatrixPointer H = (*this)(i,j);

                H->addScalMatrix(z, A, indexM, indexN);
                indexN += H->numberOfColumns();
             }        
            indexM += ((*this)(i,0))->numberOfRows(); 
         }
     }

   return (*this);
 }


HMatrixInterface& HMatrixInterface::addIdentityScaledBy (const LongReal& z)
 {
    const LongInt dim = MIN(numberOfRowBlocks(), numberOfColumnBlocks());

    for(LongInt i=0; i<dim; i++)
     {
        (*this)(i,i)->addIdentityScaledBy(z);
     }

   return (*this);
 }


HMatrixInterface& HMatrixInterface::addIdentityScaledBy (const std::complex<double>& z)
 {
    const LongInt dim = MIN(numberOfRowBlocks(), numberOfColumnBlocks());

    for(LongInt i=0; i<dim; i++)
     {
        (*this)(i,i)->addIdentityScaledBy(z);
     }

   return (*this);
 }
 

HMatrixInterface& HMatrixInterface::isMultiplicationOf(const AMatrix& A, const AMatrix& B)
 {
    //! \todo besser machen !!!
    setNullMatrix();
    update(A, B);
    
   return (*this);
 }


HMatrixInterface& HMatrixInterface::update(const AMatrix& A, const AMatrix& B)
 {
    const AMatrixType theMatrixA = A.type();
    const AMatrixType theMatrixB = B.type();

    if(theMatrixA==isHMatrixInterface && theMatrixB==isHMatrixInterface)
     {
        const HMatrixInterface& H1 = (const HMatrixInterface&)A;
        const HMatrixInterface& H2 = (const HMatrixInterface&)B;
 
        const LongInt m = H1.numberOfRowBlocks();
        const LongInt k = H1.numberOfColumnBlocks();
        const LongInt n = H2.numberOfColumnBlocks();

        for(LongInt i=0; i<m; i++)
         {
            for(LongInt j=0; j<n; j++)
             {
    
                for(int l=0; l<k; l++)
                 {
                    //H_ij += H1(i,l)*H2(l,j) 
                    (*this)(i,j)->update(H1(i,l), H2(l,j));
                 } 
             }
         }
     }
    else if(isComplexArithmetic()==true)
     { 
        if(theMatrixA==isRkCMatrix || theMatrixA==isRkCMatrix)
         { 
            const LongInt m = A.numberOfRows();
            const LongInt n = B.numberOfColumns();         
            LongInt rank = 0;

            if(theMatrixA==isRkCMatrix)
             {
                rank = ((const RkCMatrix&) A).rank(); 
             }
            else
             {
                rank = ((const RkCMatrix&) B).rank();
             }

            RkCMatrix R(m, n, rank);

            R.isMultiplicationOf(A, B);
            (*this) += R;
         }
        else
         {
            const LongInt m = A.numberOfRows();
            const LongInt n = B.numberOfColumns();
                 
            CMatrix C(m, n);
            //! \todo Multiplikation ohne Update
            C.update(A, B);

            (*this) += C;
         }
     }
    else
     {
        if(theMatrixA==isRkRMatrix || theMatrixB==isRkRMatrix)
         { 
            const LongInt m = A.numberOfRows();
            const LongInt n = B.numberOfColumns();         
            LongInt rank = 0;

            if(theMatrixA==isRkRMatrix)
             {
                rank = ((const RkRMatrix&) A).rank(); 
             }
            else
             {
                rank = ((const RkRMatrix&) B).rank();
             }

            RkRMatrix R(m, n, rank);

            R.isMultiplicationOf(A, B);
            (*this) += R;
         }
        else
         {
            const LongInt m = A.numberOfRows();
            const LongInt n = B.numberOfColumns();
                 
            RMatrix C(m, n);
            //! \todo Multiplikation ohne Update
            C.update(A, B);

            (*this) += C;
         }
     }

   return (*this);
 }

HMatrixInterface& HMatrixInterface::isInverseOf(const AMatrix& s, AMatrix& w)
 {
    const AMatrixType theMatrixS = s.type();
    const AMatrixType theMatrixW = w.type();
 
    if(theMatrixS==isHMatrixInterface && theMatrixW==isHMatrixInterface)
     {
        const LongInt blockRows = numberOfRowBlocks();
        const LongInt blockCols = numberOfColumnBlocks();

        const HMatrixInterface& S = (const HMatrixInterface&)s;
              HMatrixInterface& W = (HMatrixInterface&)w;        

        LongInt l, j;
        for(l=0; l<blockRows; l++)
         {
            AMatrixPointer Si_ll = (*this)(l,l);

            Si_ll->isInverseOf(*S(l,l), *W(l,l));

            for(j=0; j<l; j++)
             {
                W(l,j)->isMultiplicationOf(Si_ll, (*this)(l,j));
                (*this)(l,j)->isConversionOf(*W(l,j)); //(l,j) = W(l,j), (same type)
             }                                                          
            for(j=l+1; j<blockCols; j++)
             {
                W(l,j)->isMultiplicationOf(Si_ll, S(l,j));
                S(l,j)->isConversionOf(*W(l,j)); // S(l,j) = W(l,j) (same type)
             }                                                       
            for(LongInt i=l+1; i<blockRows; i++)
             {
                AMatrixPointer S_il = S(i,l);                 
                for(j=0; j<=l; j++)
                 {
                     W(i,j)->isMultiplicationOf(S_il, (*this)(l,j));
                    *W(i,j) *= -1.0;
                    (*this)(i,j)->addMatrix(*W(i,j));
                 }
                for(j=l+1; j<blockCols; j++)
                 {
                     W(i,j)->isMultiplicationOf(S_il, S(l,j));
                    *W(i,j) *= -1.0;
                    S(i,j)->addMatrix(*W(i,j));
                 }
             }//End for(LongInt i=l+1; i<blockRows; i++)
         }//End for(LongInt l=0; l<blockRows; l++)
        for(l=blockRows-1; 0<=l; l--)
         {
            for(LongInt i=l-1; 0<=i; i--)
             {
                for(LongInt j=0; j<blockCols; j++)
                 {
                    W(i,j)->isMultiplicationOf(S(i,l), (*this)(l,j));
                    *W(i,j) *= -1.0;
                    (*this)(i,j)->addMatrix(*W(i,j));
                 }
             }//End for(LongInt i=l-1; 0<=i; i--)
         }//End for(l=blockRows-1; 0<=l; l--)
     }
    else
     {
        throw SimpleException(IString("Warning In HMatrixInterface::isInverseOf, Case Not Implemented!!!"));
     }
   return (*this);
 }


HMatrixInterface& HMatrixInterface::isConversionOf(const AMatrix& A)
 {   
    const AMatrixType theMatrixA = A.type();
 
    if(theMatrixA==isHMatrixInterface)
     {
        const HMatrixInterface& H = (const HMatrixInterface&)A;

        const LongInt blockM = H.numberOfRowBlocks();
        const LongInt blockN = H.numberOfColumnBlocks();

/*
        //! \todo besser machen  
        (*this)
         .setNumberOfRowBlocks(blockM)
         .setNumberOfColumnBlocks(blockN)
         .setNumberOfRows(A.numberOfRows())
         .setNumberOfColumns(A.numberOfColumns())
        ;

        delete [] attr_block;
        attr_block = new AMatrixPointer [blockM*blockN];
*/
        for(LongInt i=0; i<blockM; i++)
         {
            for(LongInt j=0; j<blockN; j++)
             {
                ((*this)(i,j))->isConversionOf(*H(i,j));
             }
         }
     }
    else
     {
        throw SimpleException(IString("Warning In HMatrixInterface::::isConversionOf(const AMatrix& A), Case Not Implemented!!!"));
     }
    
   return(*this);
 }



bool HMatrixInterface::operator()(const LongInt& i, CVector& w, const LongInt& j, const CVector& v) const 
 {
    bool value = true;

    LongInt indexM = i;
    LongInt indexN = j;
   
    const LongInt mB = numberOfRowBlocks();
    const LongInt nB = numberOfColumnBlocks();

    for(LongInt k=0; k<mB; k++)
     {
        indexN  = j;
        for(LongInt l=0; l<nB; l++)
         {
            AMatrixPointer H = (*this)(k,l);

            H->operator ()(indexM, w, indexN, v);
            indexN += H->numberOfColumns();
         }
        indexM += ((*this)(k,0))->numberOfRows();
     }     

   return value;
 }


bool HMatrixInterface::operator()(const LongInt& i, RVector& w, const LongInt& j, const RVector& v) const 
 {
    bool value = true;

    LongInt indexM = i;
    LongInt indexN = j;
   
    const LongInt mB = numberOfRowBlocks();
    const LongInt nB = numberOfColumnBlocks();

    for(LongInt k=0; k<mB; k++)
     {
        indexN  = j;
        for(LongInt l=0; l<nB; l++)
         {
            AMatrixPointer H = (*this)(k,l);

            H->operator ()(indexM, w, indexN, v);
            indexN += H->numberOfColumns();
         }
        indexM += ((*this)(k,0))->numberOfRows();
     }     

   return value;
 }


bool HMatrixInterface::evaluateTransposeAt (const CVector& x, CVector& valueAt) const
 {
    LongInt indexM = 0;
    LongInt indexN = 0;
    
    const LongInt mB = numberOfRowBlocks();
    const LongInt nB = numberOfColumnBlocks();

    for(LongInt i=0; i<nB; i++)
     {
        indexN  = 0;
        for(LongInt j=0; j<mB; j++)
         {
            AMatrixPointer H = (*this)(j,i);

            H->partialEvaluateTransposeAt(indexM, valueAt, indexN, x); 
            indexN += H->numberOfRows();
         }
        indexM += ((*this)(i,0))->numberOfColumns();
     }
 
   return true;
 }


bool HMatrixInterface::partialEvaluateTransposeAt (const LongInt& i, CVector& w, const LongInt& j, const CVector& v) const
 {
    bool value = true;

    LongInt indexM = i;
    LongInt indexN = j;
   
    const LongInt mB = numberOfRowBlocks();
    const LongInt nB = numberOfColumnBlocks();

    for(LongInt k=0; k<nB; k++)
     {
        indexN  = j;
        for(LongInt l=0; l<mB; l++)
         {
            AMatrixPointer H = (*this)(l,k);

            H->partialEvaluateTransposeAt(indexM, w, indexN, v);
            indexN += H->numberOfRows();
         }
        indexM += ((*this)(k,0))->numberOfColumns();
     }     

   return value;
 }


bool HMatrixInterface::evaluateConjugateTransposeAt (const CVector& x, CVector& valueAt) const
 {
    LongInt indexM = 0;
    LongInt indexN = 0;
    
    const LongInt mB = numberOfRowBlocks();
    const LongInt nB = numberOfColumnBlocks();

    for(LongInt i=0; i<nB; i++)
     {
        indexN  = 0;
        for(LongInt j=0; j<mB; j++)
         {
            AMatrixPointer H = (*this)(j,i);

            H->partialEvaluateConjugateTransposeAt(indexM, valueAt, indexN, x); 
            indexN += H->numberOfRows();
         }
        indexM += ((*this)(i,0))->numberOfColumns();
     }
 
   return true;
 }


bool HMatrixInterface::evaluateTransposeAt (const RVector& x, RVector& valueAt) const
 {
    LongInt indexM = 0;
    LongInt indexN = 0;
    
    const LongInt mB = numberOfRowBlocks();
    const LongInt nB = numberOfColumnBlocks();

    for(LongInt i=0; i<nB; i++)
     {
        indexN  = 0;
        for(LongInt j=0; j<mB; j++)
         {
            AMatrixPointer H = (*this)(j,i);

            H->partialEvaluateTransposeAt(indexM, valueAt, indexN, x); 
            indexN += H->numberOfRows();
         }
        indexM += ((*this)(i,0))->numberOfColumns();
     }
 
   return true;
 }


bool HMatrixInterface::partialEvaluateTransposeAt (const LongInt& i, RVector& w, const LongInt& j, const RVector& v) const
 {
    bool value = true;

    LongInt indexM = i;
    LongInt indexN = j;
   
    const LongInt mB = numberOfRowBlocks();
    const LongInt nB = numberOfColumnBlocks();

    for(LongInt k=0; k<nB; k++)
     {
        indexN  = j;
        for(LongInt l=0; l<mB; l++)
         {
            AMatrixPointer H = (*this)(l,k);

            H->partialEvaluateTransposeAt(indexM, w, indexN, v);
            indexN += H->numberOfRows();
         }
        indexM += ((*this)(k,0))->numberOfColumns();
     }     

   return value;
 }


bool HMatrixInterface::partialEvaluateConjugateTransposeAt (const LongInt& i, CVector& w, const LongInt& j, const CVector& v) const
 {
    bool value = true;

    LongInt indexM = i;
    LongInt indexN = j;
   
    const LongInt mB = numberOfRowBlocks();
    const LongInt nB = numberOfColumnBlocks();

    for(LongInt k=0; k<nB; k++)
     {
        indexN  = j;
        for(LongInt l=0; l<mB; l++)
         {
            AMatrixPointer H = (*this)(l,k);

            H->partialEvaluateConjugateTransposeAt(indexM, w, indexN, v);
            indexN += H->numberOfRows();
         }
        indexM += ((*this)(k,0))->numberOfColumns();
     }     

   return value;
 }


HMatrixInterface& HMatrixInterface::operator *= (const std::complex<double>& alpha)
 {
    const LongInt m = numberOfRowBlocks();
    const LongInt n = numberOfColumnBlocks();

    for(LongInt i=0; i<m; i++)
     {
        for(LongInt j=0; j<n; j++)
         {
            (*this)(i,j)->operator *=(alpha);
         }
     }

   return (*this);
 }


HMatrixInterface& HMatrixInterface::operator *= (const LongReal& alpha)
 {
    const LongInt m = numberOfRowBlocks();
    const LongInt n = numberOfColumnBlocks();

    for(LongInt i=0; i<m; i++)
     {
        for(LongInt j=0; j<n; j++)
         {
            (*this)(i,j)->operator *=(alpha);
         }
     }

   return (*this);
 }


HMatrixInterface& HMatrixInterface::operator /= (const std::complex<double>& alpha)
 {
    (*this) *= 1.0/alpha;

   return (*this);
 }


AMatrix&  HMatrixInterface::setNullMatrix ()
 {
    const LongInt m = numberOfRowBlocks();
    const LongInt n = numberOfColumnBlocks();

    for(LongInt i=0; i<m; i++)
     {
        for(LongInt j=0; j<n; j++)
         {
            (*this)(i,j)->setNullMatrix();
         }
     }    

   return (*this);
 }


ostream& operator << (ostream & s, const HMatrixInterface& H)
 {
    const LongInt m = H.numberOfRowBlocks();
    const LongInt n = H.numberOfColumnBlocks();

    s << H.numberOfRows() << '\t' << H.numberOfColumns() << endl;

    for(LongInt i=0; i<m; i++)
     {
        for(LongInt j=0; j<n; j++)
         {
            s << "BlockIndex = (" << i << "," << j << ")" << endl;

            const AMatrixType theType = H(i,j)->type();

            if(theType==isHMatrixInterface)
             {
                s << *(HMatrixInterface*)(H(i,j)) << endl;
             }
            else if(theType==isRMatrix)
             {
                s << *(RMatrix*)(H(i,j)) << endl;
             }
            else if(theType==isRkRMatrix)
             {
                s << *(RkRMatrix*)(H(i,j)) << endl;
             }
            else if(theType==isCMatrix)
             {
                s << *(CMatrix*)(H(i,j)) << endl;
             }
            else if(theType==isRkCMatrix)
             {
                s << *(RkCMatrix*)(H(i,j)) << endl;
             }
            else
             {
                throw SimpleException(IString("Warning In HMatrixInterface::operator *= or /=, Bad BlockType !!!"));
             }            
         }
     }
    
   return s;
 }


LongReal HMatrixInterface::frobeniusNorm () const
 {
    const LongInt m = numberOfRowBlocks();
    const LongInt n = numberOfColumnBlocks();

    LongReal norm = 0.0, temp = 0.0;

    for(LongInt i=0; i<m; i++)
     {
        for(LongInt j=0; j<n; j++)
         {
            temp  = ((*this)(i,j))->frobeniusNorm();
            norm += temp*temp;
         }
     }
    
   return sqrt(norm);
 }


LongInt HMatrixInterface::rowIndex (const LongInt& i, const LongInt& j) const
 {
    LongInt index = 0;

    for(int l=0; l<=i; l++)
     {
        index += (*this)(l,j)->numberOfRows();
     }

   return index;
 }


LongInt HMatrixInterface::columnIndex (const LongInt& i, const LongInt& j) const
 {
    LongInt index = 0;

    for(int l=0; l<=j; l++)
     {
        index += (*this)(i,l)->numberOfColumns();
     }

   return index;
 }

