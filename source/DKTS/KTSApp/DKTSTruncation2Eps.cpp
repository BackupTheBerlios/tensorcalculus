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

// DKTSTruncation2Eps.cpp: Implementierung der Klasse DKTSTruncation2Eps.
//
//////////////////////////////////////////////////////////////////////

#include "DKTSTruncation2Eps.hpp"


DKTSTruncation2Eps::DKTSTruncation2Eps(const IString& fileName, const LongInt& rMin)
:DKTSTruncation(fileName, rMin)
 {
    
 }


DKTSTruncation2Eps::DKTSTruncation2Eps(const IString& fileName, const IString& coFileName, const IString& ortFileName, const LongInt& rMin)
:DKTSTruncation(fileName, coFileName, ortFileName, rMin)
 {
 
 }


DKTSTruncation2Eps::DKTSTruncation2Eps(const LongInt& d, const LongInt& R, const LongInt& n)
:DKTSTruncation(d, 1, R, n)
 {
 
 }
  

DKTSTruncation2Eps::~DKTSTruncation2Eps()
 {
    
 }


DKTSTruncation2Eps& DKTSTruncation2Eps::truncateF2(const LongReal& eps, const LongReal& epsN, const LongReal& preC)
 {
    const LongInt d    = attr_a.d();
    const LongInt oR   = originalRank();    
    const LongInt rMin = attr_x.k();
    const LongInt R    = attr_a.k();
    const LongInt m    = attr_a.n();
    const LongInt n    = attr_Z.n();        

    Protocol& truncationLog   = attr_truncationLog;
    Protocol& truncationLogR1 = attr_truncationLogR1;
    
    Timer time, fTimer;
    
    const IString logFile = IString("log_truncate2Eps_")  + IString(eps) + IString("_rMin=") + IString(rMin) + IString("_.log");
    ofstream log(logFile);
    
    LongReal sec=0.0, gTime=0.0, fTime=0.0;
    
    IString date(time.date());

    addInputTensorInformation(date, truncationLog);
    
    const LongReal normA = DKTSTruncation::normA();
    
    log << "Tensor Decompositions       : " << endl;
    log << "Date                        : " << date;
    log << endl;
    log << "Truncate2Eps "   << endl;
    log << endl;
				log << scientific << setprecision(5) ;
				log << "eps          = " << eps  << endl;
				log << "minPrecision = " << preC << endl;
				log << "epsNewton    = " << epsN << endl;
				log << resetiosflags( ::std::ios::scientific );
				log << endl;
    log << "d   = " << d << endl;
    log << "R   = " << oR << endl;
    log << "n   = " << n << endl;
    log << "|A| = " << normA << endl;
    log << endl;

    decomposer()
     .setEpsilon   (epsN)
     .setPrecision (preC)
     .setPrintCout (true)
    ;            
    
    attr_a.writeIndexSet(log);
         
    IString file, newFile, newFileG;
    
    time.startTiming();    

    DKTSDIterationInfoBlock infoBlock, infoBlockR1;

    infoBlock  .setCaptionString("Kronecker Rank Truncation");        
    infoBlockR1.setCaptionString("BestR1Truncation");

    LongReal error = 1.0e20;
    
    DKTSDIterationInfo infoEntry;
    
    infoEntry.setError(error);
    
    LongInt  r = 1;

    log << "kRank"    << '\t' << "|A-X_0|/|A|" << '\t' << "|A-X|/|A|"   << '\t' << "relDiff.[%]"	<< '\t';
				log << "Gradient" << '\t' << "Steps"      << '\t' << "time[sec.]" << endl;
								
				
				if(rMin==1)
     {    
        newFile  = IString("tensor_best_r=") + IString(rMin) + IString("_.ten");
        newFileG = IString("igTensor_best_r=") + IString(r) + IString("_.ten");

/*        
        (*this)
         .writeSolutionTo(newFileG)
        ;
*/

        infoEntry = truncate();                
        error     = infoEntry.error();
                
        infoEntry.setStep(rMin);       
	       log << infoEntry;
	
	       (*this)
         .writeSolutionTo(newFile)
        ;
       
        infoBlock.addEntry(infoEntry);       

        r = rMin + 1;        
     }
    else
     {
        newFile = IString("tensor_best_r=") + IString(rMin-1) + IString("_.ten");
	
        (*this)
         .readInitialGuestFrom(newFile)
        ;
               
        error = attr_x.distanceTo(attr_a)/normA;
        r     = rMin;
     }       
    
    DKTS residuum, x0(d, 1, m); 
 
    while(eps<error && r<=R)
     {
        newFile  = IString("tensor_best_r=") + IString(r) + IString("_.ten");

        residuum.setSumOf(1.0, attr_a, -1.0, attr_x);
        
        DKTSDIterationInfo infoEntryR1;
		
        infoEntryR1 = bestR1Truncation(r, residuum, x0, normA*error);
	       //infoEntryR1 = DKTSTruncation::truncate2(2, residuum, x0, normA*error);	       
	       //x0 *= 1.0/error;
	
       	DKTS xt(attr_x);
        
	       attr_x.setSumOf(1.0, xt, 1.0, x0);        

        r = attr_x.k();
        infoEntryR1.setStep(r);

        newFileG = IString("igTensor_best_r=") + IString(r) + IString("_.ten");
/*                
        (*this)        
         .writeSolutionTo(newFileG)
        ;
*/
        infoEntry = truncate();
	       error     = infoEntry.error();
                        
        cout << endl;
        cout << "Write Solution to File : " << newFile << " (";

        fTimer.startTiming();

        (*this)     
         .writeSolutionTo(newFile)
        ;                                        
        
        fTime = fTimer.elapsedTimeSec();
        
	       cout << "Time [sec] = " << fTime << ")"<< endl << endl;
               
        infoEntry.setStep(r);
	
	       log << infoEntry;
								
	       infoBlock.addEntry(infoEntry);
        infoBlockR1.addEntry(infoEntryR1);
       
        r++;
     }    
    
    gTime = time.elapsedTimeSec();
    

    cout << endl;
    cout << "Total Time [sec] = " << setw(8) << gTime << endl;    

    addInfoBlock  (infoBlock, gTime, eps, epsN, preC);
    addInfoBlockR1(infoBlockR1);               

    log << endl;
    log << "Total Time [sec] = " << setw(8) << gTime << endl;    

   return (*this);
 }


DKTSTruncation2Eps& DKTSTruncation2Eps::truncate2(const LongReal& eps, const DKTS& A, DKTS& X, const bool& plotInfo, const bool& useX)
 {
 
 
   return (*this);
 }
