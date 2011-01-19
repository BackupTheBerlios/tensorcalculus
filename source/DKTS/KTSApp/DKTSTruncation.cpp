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

// DKTSTruncation.cpp: Implementierung der Klasse DKTSTruncation.
//
//////////////////////////////////////////////////////////////////////

#include "DKTSTruncation.hpp"


DKTSTruncation::DKTSTruncation(const LongInt& d, const LongInt& r, const LongInt& R, const LongInt& n)
 {
    (*this)
     .allocateDataSpace(d, r, R, n)
     .generateRandomExample(d, r, R, n)
					.setDefaultParameter()
    ;
 }
 
 
DKTSTruncation::DKTSTruncation(const IString& fileName, const LongInt& r)
 {
    DKTS A;        
    
    A.readDataFrom(fileName);        
                       
    const LongInt R = A.k();
               
    attr_indexSet = (LongIntPointer)  new LongInt[R];
    attr_values   = (LongRealPointer) new LongReal[R];
 
    (*this)
					.setInputFileName(fileName)
     .prepareInputData(A)
     .setInitialGuestBiggestSum(r)
					.setDefaultParameter()
    ; 
 }


DKTSTruncation::DKTSTruncation(const IString& fileName, const IString& coFileName, const IString& ortFileName, const LongInt& r)
 {        
    attr_Z
     .readDataFrom(ortFileName)
    ;
    
    attr_a
     .readDataFrom(coFileName)
    ;

    const LongInt R = attr_a.k();

    (*this)
     .setOriginalRank(R)
    ;

    attr_indexSet = (LongIntPointer)  new LongInt[R];
    attr_values   = (LongRealPointer) new LongReal[R];

    (*this)
     .setInputFileName(fileName)
     .sortIndexSet()
     //.computeMainPart(1.0e-16)
     .setInitialGuestBiggestSum(r)
					.setDefaultParameter()
    ;

 }


DKTSTruncation::DKTSTruncation(const DKTS& A, const LongInt& r)
 {
    const LongInt R = A.k();
    
    attr_indexSet = (LongIntPointer)  new LongInt[R];
    attr_values   = (LongRealPointer) new LongReal[R];
 
    (*this)
     .prepareInputData(A)
     .sortIndexSet()
     .setInitialGuestBiggestSum(r)
					.setDefaultParameter()
    ;  
 }


DKTSTruncation::DKTSTruncation(const DKTS& a, const DKTS& A, const DKTS& Z, const LongInt& r)
 {
    attr_a = a;
    attr_Z = Z;
    
    const LongInt R = attr_a.k();
    
    attr_indexSet = (LongIntPointer)  new LongInt[R];
    attr_values   = (LongRealPointer) new LongReal[R];
 
    (*this)
     .sortIndexSet()
     .setInitialGuestBiggestSum(r)
					.setDefaultParameter()
    ; 
 }

     
DKTSTruncation::~DKTSTruncation()
 {
   (*this)
    .deleteDataSpace()
   ;
 }


DKTSTruncation& DKTSTruncation::setDefaultParameter()
 {				
			 ProtocolFormat pF;
				
			 pF
				 .setTablePosition("H")
				 .setMaxTableSize(25)
			 ;
				
			 ProtocolProperties pP(pF);
				
			 pP
			  .setTopicString("R1 Updates")
				 .setTeXFileName("lastRunR1.tex")
			 ;	
	
	   attr_truncationLogR1.setProperties(pP);
				
		 return (*this);
	}		


DKTSDIterationInfo DKTSTruncation::bestR1Truncation(const LongInt& r, DKTS& a, DKTS& x, const LongReal& normD)
 {    
    Random rGenerator;    
      
    decomposer()    
     .setPrintCout(false)
    ;

    LongReal error    = 1.0e20;
    LongReal errorMin = error;
    
    Timer time;        

    const LongInt d = a.d();
    const LongInt n = a.n();

    x.resize(d, 1, n);

    DKTS x0(x);

    DKTSDIterationInfo infoEntry;
    
    time.startTiming();
    
    const LongInt R  = attr_a.k();
    const LongInt mR = MIN(7+1, R);
    
    LongInt index = 0;
                
    LongInt max     = 0;
    LongInt maxR    = 23;
    LongInt quality = -1;    
    
    bool noOne = true;
    
    DKTSDDataBlock dataBlock;
    
    addBeginTrancationR1(r, attr_truncationLogR1);
    
    cout << "Initial guess : ";
    
    for(LongInt i=1; i<mR || max==maxR; i++)
     {
        cout << ".";
        
        addTrancationR1(i, r, attr_truncationLogR1);
        
        const LongInt rIndex = rGenerator.randomLongInt(R);
                
        index = rIndex;

        x0.copyFromTo(attr_a, -index, index, index);
    
        const LongReal normX0 = frobeniusNorm(x0);
        
        x0 *= (normD/normX0);           
    
        dataBlock.removeAllEntrys();
        infoEntry = attr_decomposer.decompose(a, x0, normD, dataBlock);                
            
        attr_truncationLogR1.add(dataBlock);    
            
        error = infoEntry.error();    
        
        if(1.0 < error && noOne)
         {
            i--;
         }
        else
         {                
            if(error < errorMin)
             {
                noOne    = false;
                x = x0;
                errorMin = error;
                quality  = attr_decomposer.quality();
                
                if(errorMin < 0.92)
                 {
                    i = mR;
                 }
             }
         }                        
     }        
         
    LongReal sec = time.elapsedTimeSec();
    cout << errorMin << ", (" << quality << "), " << sec << "." << endl << endl;
                               
    addEndTrancation(infoEntry, r, attr_truncationLogR1);
    
    infoEntry.setCalculationTime(sec);

    decomposer()
     .setPrintCout(true)
    ;

   return infoEntry;
 }


/*
DKTSDIterationInfo DKTSTruncation::bestR1Truncation(const LongInt& r, DKTS& a, DKTS& x, const LongReal& normD)
 {    
    decomposer()    
     .setPrintCout(false)
    ;

    LongReal error    = 1.0e20;
    
    Timer time;        

    const LongInt d = a.d();
    const LongInt R = a.k();
    const LongInt n = a.n();
    
    LongInt quality = -1;    

    x.resize(d, 1, n);


    for(LongInt mu=0; mu<d; mu++)
     {
        DKTVector& v = x(0,mu);
        
        for(LongInt i=0; i<R; i++)
         {
            v += a(i,mu);
         }        
        v.normalized();
     }

    //x *= normD;    

    DKTSDIterationInfo infoEntry;
    
    time.startTiming();
    
    DKTSDDataBlock dataBlock;
    
    addBeginTrancationR1(r, attr_truncationLogR1);
    
    cout << "Initial guess : ";    
        
    addTrancationR1(1, r, attr_truncationLogR1);
                          
    dataBlock.removeAllEntrys();
    
    infoEntry = attr_decomposer.decompose(a, x, normD, dataBlock);                
            
    attr_truncationLogR1.add(dataBlock);    
            
    error   = infoEntry.error();
    quality = attr_decomposer.quality();    
                 
    LongReal sec = time.elapsedTimeSec();
    cout << error << ", (" << quality << "), " << sec << "." << endl << endl;
                               
    addEndTrancation(infoEntry, r, attr_truncationLogR1);
    
    infoEntry.setCalculationTime(sec);

    decomposer()
     .setPrintCout(true)
    ;

   return infoEntry;
 }
*/

DKTSDIterationInfo DKTSTruncation::truncate2(const LongInt& rT, DKTS& a, DKTS& x, const LongReal& normA)//, Protocol& protocol)
 { 
    const LongInt d = a.d();
    const LongInt R = a.k();
    const LongInt n = a.n();
    
    DKTSDIterationInfo infoEntry;            
    DKTSDDataBlock dataBlock;
    
    infoEntry = bestR1Truncation(1, a, x, normA);
    
    LongReal error = infoEntry.error();
    
    DKTS x0(d, 1, n), residuum(d, R, n);
    
    for(LongInt r=1; r<rT; r++)
     {     
        residuum.setSumOf(1.0, a, -1.0, x);        
	
        bestR1Truncation(r, residuum, x0, normA*error);
		
       	DKTS xt(x);
        
	       x.setSumOf(1.0, xt, 1.0, x0);
        x.reScaled();
   
        decomposer()
         .setPrintCout(false)
        ;
              
        infoEntry = attr_decomposer.decompose(a, x, normA, dataBlock);    
        
        error = infoEntry.error();
     }

    decomposer()
     .setPrintCout(true)
    ;
    
   return infoEntry; 
 }

DKTSDIterationInfo DKTSTruncation::startTruncation(DKTS& a, DKTS& x, const LongReal& normA)
 {    
    DKTSDIterationInfo infoEntry;
    Timer time;
    DKTSDDataBlock dataBlock;    

    time.startTiming();    
    infoEntry = attr_decomposer.decompose(a, x, normA, dataBlock);    
    LongReal sec = time.elapsedTimeSec();    
     
    if(attr_decomposer.printCout())
     {
        cout << "Time [sec] = " << sec << endl;
     }

    infoEntry.setCalculationTime(sec);    
   
    addDataBlock(dataBlock, infoEntry, x.k(), attr_truncationLog);

   return infoEntry;
 }


DKTSDIterationInfo DKTSTruncation::truncate()
 {        
    if(attr_decomposer.printCout())
     {
        writeParameter(cout);
     }
                 
   return startTruncation(attr_a, attr_x, normA());
 }
 
 
DKTSTruncation& DKTSTruncation::writeSolutionTo(const IString& fileName)
 {
    DKTS X;
    
    X
     .regeneratedBy(attr_x, attr_Z)
     .reScaled()
     .writeDataTo(fileName)
    ;
        
   return (*this);
 }
 
 
DKTSTruncation& DKTSTruncation::allocateDataSpace(const LongInt& d, const LongInt& r, const LongInt& R, const LongInt& n)
 {
    attr_indexSet = (LongIntPointer)  new LongInt [R];
    attr_values   = (LongRealPointer) new LongReal[R];
     
    attr_Z.resize(d, R, n);    
    
    attr_a.resize(d, R, R);
    attr_x.resize(d, r, R);
    
   return (*this);
 }


DKTSTruncation& DKTSTruncation::deleteDataSpace()
 {
    delete [] attr_indexSet;
    delete [] attr_values;
    
   return (*this);
 }


DKTSTruncation& DKTSTruncation::resize(const LongInt& d, const LongInt& r, const LongInt& R, const LongInt& n)
 {
    const LongReal d1 = attr_a.d();
    const LongReal R1 = attr_a.k();
    const LongReal r1 = attr_x.k();
    const LongReal n1 = attr_Z.n();
    
    if(d1==d && r1==r && R1==R && n1==n)
     {
        attr_a.setNull();
        attr_x.setNull();
        attr_Z.setNull();        
     }
    else
     {
        (*this)
         .deleteDataSpace()
         .allocateDataSpace(d, r, R, n)
        ;
     }
     
   return (*this);
 }


DKTSTruncation& DKTSTruncation::prepareInputData(const DKTS& A)
 {
    (*this)
     .prepareInputData(A, 1)
    ;
         
   return (*this);
 }

 
DKTSTruncation& DKTSTruncation::prepareInputData(const DKTS& A, const LongInt& r)
 {
    const LongInt d = A.d();
    const LongInt R = A.k();
    const LongInt n = A.n();
    
    (*this)
     .resize(d, r, R, n)
     .setOriginalRank(R)
    ;
    
    Timer time;
    time.startTiming();            
    
    IString date(time.date());
    
    attr_Z.setOrthogonal2(A);
    attr_a.setCoefficientsSystemOf(A, attr_Z);        
    
    const LongReal sec = time.elapsedTimeSec();
    
    (*this)
     .setPreCalculationTime(sec)
     .sortIndexSet()
     //.computeMainPart(1.0e-16);
    ;
   
   return (*this);
 }


DKTSTruncation& DKTSTruncation::computeMainPart(const LongReal& eps)
 {
    const LongInt R = attr_a.k();

    // compute the new Rank
    
    LongReal normA = DKTSTruncation::normA();       
    LongReal epsA  = eps*normA;    
    LongReal rest  = attr_values[attr_indexSet[R-1]];

    LongInt diff = 0;

    LongInt i = R-2;
            
    while(rest<epsA && 0<i)
     {
       rest += attr_values[attr_indexSet[i]];

       diff++;
       i--;
     }    
    
    const LongInt rNew = R - diff;
    
    if(0<rNew)
     {
        setOriginalRank(R);
    
        LongIntPointer  indexSet = (LongIntPointer)  new LongInt [rNew];
        LongRealPointer values   = (LongRealPointer) new LongReal[rNew];    
    
        for(i=0; i<rNew; i++)
         {
            indexSet[i] = attr_indexSet[i];
            values[i]   = attr_values[i];
         }
    
        delete [] attr_indexSet;
        delete [] attr_values;

        attr_indexSet = (LongIntPointer)  new LongInt [rNew];
        attr_values   = (LongRealPointer) new LongReal[rNew];
    
        for(i=0; i<rNew; i++)
         {
            attr_indexSet[i] = i;
            attr_values[i]   = values[i];
         }
     
        DKTS a(attr_a);    
    
        const LongInt d = a.d();
    
        attr_a.resize(d, rNew, a.n());
        
        for(LongInt mu=0; mu<d; mu++)
         {
            for(i=0; i<rNew; i++)
             {            
                attr_a(i, mu) = a(indexSet[i], mu);
             }
         }    
    
        delete [] indexSet;
        delete [] values;
     }
    
   return (*this);
 }

 
DKTSTruncation& DKTSTruncation::prepareInputData(const DKTS& A, const DKTS& X)
 {
    const LongInt d = A.d();
    const LongInt R = A.k();
    const LongInt r = X.k();
    const LongInt n = A.n();
    
    (*this)
     .resize(d, r, R, n)
    ;
 
    Timer time;
    time.startTiming();
    
    IString date(time.date());
    
    attr_Z.setOrthogonal2(A);
    attr_a.setCoefficientsSystemOf(A, attr_Z);
    attr_x.setCoefficientsSystemOf(X, attr_Z);
    attr_x.reScaled();
    
    const LongReal sec = time.elapsedTimeSec();
        
    (*this)
     .sortIndexSet()
     .setPreCalculationTime(sec)
     //.computeMainPart(1.0e-16);     
    ;
      
   return (*this);
 }


DKTSTruncation& DKTSTruncation::generateRandomExample(const LongInt& d, const LongInt& r, const LongInt& R, const LongInt& n, const LongReal& eps)
 {
    DKTS A(d, R, n);
    
    A
     .setRand()
     .scale(eps)
    ;        
    
    (*this)
     .prepareInputData(A, r)
     .setInitialGuestBiggestSum(r)
    ;

   return (*this); 
 }


DKTSTruncation& DKTSTruncation::writePreComputationDataTo (const IString& fileNameCt, const IString& fileNameOb)
 {
    (*this)
     .writeCoefficientSystemTo(fileNameCt)
     .writeOrthogonalBasisTo(fileNameOb)
    ;
    
   return (*this);
 }


DKTSTruncation& DKTSTruncation::resizeInitialGuest(const LongInt& d, const LongInt& r, const LongInt& R, const LongInt& n)
 {
    const LongInt rank = MIN(r, R);
    
    attr_x.resize(d, rank, n);
    
   return (*this);
 }


DKTSTruncation& DKTSTruncation::writeCoefficientSystemTo(const IString& fileName)
 {
    attr_a.writeDataTo(fileName);
    
   return (*this);
 }


DKTSTruncation& DKTSTruncation::writeOrthogonalBasisTo(const IString& fileName)
 {
    attr_Z.writeDataTo(fileName);
    
   return (*this);
 }
 

DKTSTruncation& DKTSTruncation::setInitialGuestBiggestSum(const LongInt& r)
 {
    resizeInitialGuest(attr_a.d(), r, attr_a.k(), attr_a.n());    
 
    const LongInt k = MIN(attr_x.k(), attr_a.k());
    const LongInt d = attr_x.d();

    for(LongInt j=0; j<k; j++)
     {
        const LongInt index = attr_indexSet[j];
        
        for(LongInt mu=0; mu<d; mu++)
         {                          
            attr_x(j,mu) = attr_a(index,mu);
         }
     }
    
    attr_x.reScaled();

   return (*this); 
 }
 
 
DKTSTruncation& DKTSTruncation::readInitialGuestFrom(const IString& fileName)
 {     
    DKTS X;
    X.readDataFrom(fileName);
    X.writeIndexSet(cout);
    
    attr_x.setCoefficientsSystemOf(X, attr_Z);   
    
    attr_x.reScaled();
    
   return (*this);
 }
 

DKTSTruncation& DKTSTruncation::sortIndexSet()
 {
    const LongInt R = attr_a.k();
    
    const LongReal normA = frobeniusNorm(attr_a);
    
    (*this)
     .setNormA(normA)
    ;    
    
    for(LongInt i=0; i<R; i++)
     {
        attr_values[i]   = attr_a.frobeniusNormOfSummand(i);        
        attr_indexSet[i] = i;
     }

    quickSort(0, R-1, attr_values);        

   return (*this);
 }

 
LongInt DKTSTruncation::partition(LongInt low, LongInt high, LongRealPointer f) 
 {
    LongInt i,j;
    LongReal pivot = f[attr_indexSet[low]];
    i = low;

    for(j=low+1; j<=high; j++) 
     {
       if (pivot<=f[attr_indexSet[j]])
        {
           i++;
           swapIndex(i, j);
        }
     }

    swapIndex(i, low);

   return i;
 }


DKTSTruncation& DKTSTruncation::quickSort(LongInt low, LongInt high, LongRealPointer f)
 {
    LongInt m = 0;

    if (low < high) 
     {
        m = partition(low, high, f);
        quickSort(low, m-1, f);
        quickSort(m+1, high, f);
     }

   return (*this);
 }


DKTSTruncation& DKTSTruncation::swapIndex (const LongInt& i, const LongInt& j)
 {
    LongInt& a_i = attr_indexSet[i];
    LongInt& a_j = attr_indexSet[j];

    const LongInt temp = a_i;

    a_i = a_j;
    a_j = temp;

   return (*this);
 }


bool DKTSTruncation::writeParameter(ostream& s) const
 {
    bool value = true;

    const LongInt d = attr_x.d();
    const LongInt k = attr_x.k();
    const LongInt m = attr_x.n();
    const LongInt l = originalRank();
    const LongInt n = attr_Z.n();

    s << "d = " << d << endl;
    s << "R = " << l << endl;
    s << "r = " << k << endl;
    s << "N = " << n << endl;
    s << "m = " << m << endl;    

    s << setprecision(2);
    s << "Data Storage Memory = " << setw(4) << (LongReal)( (2*d*l*n + d*k*n + d*k*l + d*l*l)*sizeof(LongReal))/(LongReal)1048576 << " MByte" << endl;

   return value;
 }


DescriptionData DKTSTruncation::parameterDescription() const
 {
    DescriptionData dD;

    const LongInt d = attr_x.d();
    const LongInt k = attr_x.k();
    const LongInt m = attr_x.n();
    const LongInt l = attr_a.k();
    const LongInt n = attr_Z.n();

    dD.addString (IString("d ") + IString("$ = ") + IString(d) + IString("$,") + IString("\\hspace{15pt}")
                + IString("R ") + IString("$ = ") + IString(l) + IString("$,") + IString("\\hspace{15pt}")
	               + IString("r ") + IString("$ = ") + IString(k) + IString("$,") + IString("\\hspace{15pt}")
	               + IString("N ") + IString("$ = ") + IString(n) + IString("$,") + IString("\\hspace{15pt}")
	               + IString("m ") + IString("$ = ") + IString(m) + IString("$"));
    dD.addString ("");
    
   return dD;
 }

 
bool DKTSTruncation::writeNormsOfA(DescriptionData& dDIn) const
 {
    const LongInt R = attr_a.k();

    for(LongInt i=0; i<R; i++)
     {
        dDIn.addString (IString(attr_indexSet[i]) + IString(" : ") 
	                     + IString(attr_a.frobeniusNormOfSummand(attr_indexSet[i]))
		                    + IString("\\\\"));
     }

   return true;
 }


bool DKTSTruncation::writeNormsOfA(ostream& s) const
 {
    const LongInt R = attr_a.k();
  
    for(LongInt i=0; i<R; i++)
     {
        s << setw(5) << attr_indexSet[i] << " : " << attr_a.frobeniusNormOfSummand(attr_indexSet[i]) << endl;
     }
  
   return true;
 }


DKTSTruncation& DKTSTruncation::addInputTensorInformation(const IString& date, Protocol& truncationLog)
 {
    const LongInt d    = attr_a.d();
    const LongInt R    = attr_a.k();
    const LongInt l    = attr_a.n();
    const LongInt k    = attr_x.k();
    const LongInt n    = attr_Z.n();
    const LongInt oR   = originalRank();


    DescriptionData dD;    
    dD.addString(IString("\\chapter{Input Tensor Sum Information}"));
    dD.addString(IString("Date : ")                             + IString(date)       + IString("\\\\"));
    // toDo kein _ im FileName
    dD.addString(IString("Reading Initial-Tensor from : ")      + IString("$")        + IString("tensor\\_input.ten") + IString("$"));
    dD.addString(IString(""));
    dD.addString(IString("Time for Data Precalculation [sec.]") + IString("$\\ =\\ ") + IString(preCalculationTime()) + IString("$"));
    dD.addString(IString(""));

    dD.addString( IString("$\\|A\\| = $") + IString(normA()));
    
    dD.addString("");
    dD.addString (IString("d  ") + IString("$ = ") + IString(d) + IString(",$") + IString("\\hspace{15pt}")
                + IString("R  ") + IString("$ = ") + IString(R) + IString(",$") + IString("\\hspace{15pt}")
                + IString("oR ") + IString("$ = ") + IString(oR) + IString(",$") + IString("\\hspace{15pt}")
                + IString("n  ") + IString("$ = ") + IString(n) + IString("$"));
    dD.addString("");

    dD.addString(  IString("Data\\ Storage\\ Memory\\ ") + IString("$\\ =\\ ")
                 + IString((LongReal)( (2*d*l*n + d*k*n + d*k*l + d*l*l)*sizeof(LongReal))/(LongReal)1048576)
       	         + IString("$") + IString("\\ MByte"));

    dD.addString("");

    dD.addString(IString("$i:\\|A_i\\|$") + IString("\\\\"));
    
    (*this)
     .writeNormsOfA(dD)
    ;
    
    dD.addString(IString(""));
    dD.addString(IString("$i:\\|A_i\\|$") + IString("\\\\"));
    
    attr_a.writeIndexSet(dD);
    
    truncationLog.add(dD);    

   return (*this);
 }


DKTSTruncation& DKTSTruncation::addBeginTrancation (const LongInt& r, Protocol& truncationLog)
 {    
    DescriptionData dD;
	
    dD.addString(IString("\\chapter{Best Rank ") + IString(r) + IString(" Truncation}"));
    truncationLog.add(dD);    

   return (*this);
 }


DKTSTruncation& DKTSTruncation::addBeginTrancationR1(const LongInt& r, Protocol& truncationLog)
 {    
    DescriptionData dD;
	
    dD.addString(IString("\\chapter{Compute Initial Guest for Best ") + IString(r) + IString(" Truncation}"));
    truncationLog.add(dD);    

   return (*this);
 }


DKTSTruncation& DKTSTruncation::addTrancationR1(const LongInt& run ,const LongInt& rank, Protocol& truncationLog)
 {    
    DescriptionData dD;
	
    dD.addString (IString("\\section{$ ") + IString(run)  + IString(" ^{st}$ Test ") 
				            + IString("Best ")        + IString(rank) +	IString(" Truncation}"));
    truncationLog.add(dD);    

   return (*this);
 }
 
 
DKTSTruncation& DKTSTruncation::addDataBlock(const DKTSDDataBlock& dataBlock, const DKTSDIterationInfo& infoBlock, 
                                             const LongInt& k, Protocol& truncationLog)
 {
   addBeginTrancation(k, truncationLog);
   truncationLog.add(dataBlock);
   addEndTrancation(infoBlock, k, truncationLog);

  return (*this);
 }


DKTSTruncation& DKTSTruncation::addEndTrancation(const DKTSDIterationInfo& infoBlock, const LongInt& r, Protocol& truncationLog)
 {    
    DescriptionData dD;

    const LongInt  nSte = infoBlock.numberOfNewtonSteps();
    const LongReal dOld = infoBlock.startError();
    const LongReal dNew = infoBlock.error();
    const LongReal diff = dOld-dNew;
    const LongReal time = infoBlock.calculationTime();


    dD.addString(IString("Working Memory\\ \\ =\\ ") + IString("$\\ \\ ") + IString((LongReal)(attr_decomposer.memory()*sizeof(LongReal))/(LongReal)1048576) + IString("$ MByte"));    
    dD.addString("\\\\");
    dD.addString(IString("Time\\ [sec]\\      =\\ ") + IString("$\\ \\ ") + IString(time) + IString("$"));
    dD.addString("");

    dD.addString(  IString("$") + IString("\\frac{\\|A-X_0\\|}{\\|A\\|} = ")   + IString(dOld) + IString("$") + IString("\\hspace{15pt}")
                 + IString("$") + IString("\\frac{\\|A-X_{")    + IString(nSte)
                 + IString("}\\|}{\\|A\\|} = ") + IString(dNew) + IString("$") + IString("\\hspace{15pt}") 
                 + IString("\\\\"));

    dD.addString(  IString("diff ")      + IString("$\\ =\\ ") + IString(diff)          + IString("$") + IString("\\hspace{15pt}")
                 + IString("diff[\\%] ") + IString("$\\ =\\ ") + IString(diff/dOld*100) + IString("$") + IString("\\%"));

    dD.addString("");
    dD.addMathString("i:\\|X_i\\|");
    dD.addString("");

    attr_x.writeIndexSet(dD);
        	           			
    //dD.addString(IString("Write Solution to\\ File\\ :\\ ") + IString("tensor\\_best\\_r=") + IString(r) + IString("\\_.ten"));

    truncationLog.add(dD);    

   return (*this);
 }


DKTSTruncation& DKTSTruncation::addInfoBlock(const DKTSDIterationInfoBlock& infoBlock, const LongReal& totalTime, 
                                             const LongReal& eps, const LongReal& epsN, const LongReal& preC)
 {
    const LongInt d    = attr_a.d();
    const LongInt R    = originalRank();
    const LongInt l    = attr_a.n();
    const LongInt k    = attr_x.k();
    const LongInt n    = attr_Z.n();        

    DescriptionData dD, dE;
    
    dD.addString(IString("\\chapter{Summary}"));
    
    dE.addString (IString("d ") + IString("$ = ") + IString(d) + IString(",$") + IString("\\hspace{15pt}")
                + IString("R ") + IString("$ = ") + IString(R) + IString(",$") + IString("\\hspace{15pt}")
                + IString("n ") + IString("$ = ") + IString(n) + IString("$"));
    dE.addString ("");
    
    dE.addString (IString("eps ")          + IString("$ = ") + IString(eps)  + IString("$") + IString("\\\\"));
    dE.addString (IString("minPrecision ") + IString("$ = ") + IString(preC) + IString("$") + IString("\\\\"));
	   dE.addString (IString("epsNewton ")    + IString("$ = ") + IString(epsN) + IString("$") + IString("\\\\"));
    
    dE.addString (IString("Total Time\\ [sec]\\ \\ ")  + IString("$\\ =\\ ") + IString(totalTime) + IString("$"));
    
    attr_truncationLog.add(dD);
    attr_truncationLog.add(infoBlock);
    attr_truncationLog.add(dE);
  
   return (*this);
 }


 DKTSTruncation& DKTSTruncation::addInfoBlockR1(const DKTSDIterationInfoBlock& infoBlock)
 {
    DescriptionData dD;

    dD.addString(IString("\\appendix"));
    dD.addString(IString("\\chapter{Best Rank 1  Updates}"));

    attr_truncationLog.add(dD);    
    attr_truncationLog.add(infoBlock);

   return (*this);
 }


