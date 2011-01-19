/*
 * Copyright (C) Marcel Schuster
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

// DKTSDDataBlock.cpp: Implementierung der Klasse DKTSDDataBlock.
//
//////////////////////////////////////////////////////////////////////

#include "DKTSDDataBlock.hpp"


DKTSDDataBlock::DKTSDDataBlock()
 {
    (*this)
     .setDefaultParameter()
    ;
 }


DKTSDDataBlock::DKTSDDataBlock(const DKTSDDataBlock& data)
 {
    (*this) = data;
 }    
  
  
DKTSDDataBlock::~DKTSDDataBlock()
 {
 
 }


DKTSDDataBlock DKTSDDataBlock::setDefaultParameter()
 {
    (*this)
     .setStepString                ("Step")
     .setNormOfGradientString      ("relGrad")
     .setDeltaString               ("relDelta")
     .setStepSizeString            ("a_k")
     .setErrorString               ("relError")
					.setRelFunctionValueString    ("relFunctionValue")
     .setStepTeXString             ("Step")
     .setNormOfGradientTeXString   ("\\frac{\\|f'(x_k)|}{\\|A\\|}")
     .setDeltaTeXString            ("\\frac{\\|\\Delta_k\\|}{\\|A\\|}")
     .setStepSizeTeXString         ("a_k")
     .setErrorTeXString            ("\\frac{\\|A-X_{k}\\|}{\\|A\\|}")
					.setRelFunctionValueTeXString ("\\frac{f(x_k)}{\\|A\\|}")
    ;
    
   return (*this);
 }  
 
 
DKTSDDataBlock& DKTSDDataBlock::operator = (const DKTSDDataBlock& data) 
 {
    (*this)
     .setStepString                (data.stepString())
     .setNormOfGradientString      (data.normOfGradientString())
     .setDeltaString               (data.deltaString())
     .setStepSizeString            (data.stepSizeString())
     .setErrorString               (data.errorString())
					.setRelFunctionValueString    (data.relFunctionValueString())
     .setStepTeXString             (data.stepTeXString())
     .setNormOfGradientTeXString   (data.normOfGradientTeXString())
     .setDeltaTeXString            (data.deltaTeXString())
     .setStepSizeTeXString         (data.stepSizeTeXString())
     .setErrorTeXString            (data.errorTeXString())
					.setRelFunctionValueTeXString (data.relFunctionValueTeXString())
     .setCaptionString             (data.captionString())
    ;
    
    const LongInt size = data.numberOfEntrys();
    
    for(LongInt i=0; i<size; i++)
     {
        (*this).addEntry(data.getEntryAt(i));
     }
  
   return (*this);
 } 


ostream& operator << (ostream& s, const DKTSDDataBlock& dD)
 {
	   const IString stString  = dD.stepString();
				const IString erString  = dD.errorString();
    const IString nogString = dD.normOfGradientString();
    const IString delString = dD.deltaString();
    const IString sSString  = dD.stepSizeString();
				const IString rFVString = dD.relFunctionValueString();
				
				s << stString << '\t' << erString << '\t' << nogString << '\t' << delString << '\t' << sSString << '\t' << rFVString << endl;
				
    const LongInt size = dD.numberOfEntrys();
    
				for(LongInt i=0; i<size; i++)
				 {
        const DKTSDData data = dD.getEntryAt(i);
								
								s << data;
					}
					
			return s;
	}	
				 
					
bool DKTSDDataBlock::plot(ostream& log, const ProtocolProperties& pP) const
 {
    bool value = false;

    plotEntry(log, pP);
    
   return value;
 }
 
 
bool DKTSDDataBlock::plotEntry(ostream& log, const ProtocolProperties& pP) const
 {
    bool value = false;
    
    //Bestimmung der Laenge der Strings
    
    const LongInt prec = pP.format().precision();
    
    const LongInt lengthStepString             = stepString()             .length();
    const LongInt lengthNormOfGradientString   = normOfGradientString()   .length();
    const LongInt lengthDeltaString            = deltaString()            .length();
    const LongInt lengthStepSizeString         = stepSizeString()         .length();
    const LongInt lengthErrorString            = errorString()            .length();
				const LongInt lengthRelFunctionValueString = relFunctionValueString() .length();
    
    const LongInt lengthOfNumber = (pP.format().precision() + 6);
    
    LongInt disStepString,             disStepNumber;
    LongInt disNormOfGradientString,   disNormOfGradientNumber;
    LongInt disDeltaString,            disDeltaNumber;
    LongInt disStepSizeString,         disStepSizeNumber;
    LongInt disErrorString,            disErrorNumber;
				LongInt disRelFunctionValueString, disRelFunctionValueNumber; 
    
    //Berechnung der Werte zur Zentrierung
    
    centeringEntrys(4,              lengthStepString,             disStepString,             disStepNumber);
    centeringEntrys(lengthOfNumber, lengthNormOfGradientString,   disNormOfGradientString,   disNormOfGradientNumber);
    centeringEntrys(lengthOfNumber, lengthDeltaString,            disDeltaString,            disDeltaNumber);
    centeringEntrys(lengthOfNumber, lengthStepSizeString,         disStepSizeString,         disStepSizeNumber);
    centeringEntrys(lengthOfNumber, lengthErrorString,            disErrorString,            disErrorNumber);
				centeringEntrys(lengthOfNumber, lengthRelFunctionValueString, disRelFunctionValueString, disRelFunctionValueNumber);
    
    const LongInt dis  = MAX(4,              lengthStepString);
    const LongInt dis1 = MAX(lengthOfNumber, lengthNormOfGradientString); 
    const LongInt dis2 = MAX(lengthOfNumber, lengthDeltaString);
    const LongInt dis3 = MAX(lengthOfNumber, lengthStepSizeString);
    const LongInt dis4 = MAX(lengthOfNumber, lengthErrorString);
				const LongInt dis5 = MAX(lengthOfNumber, lengthRelFunctionValueString);

    //Ausgabe der Strings als Tabellenkopf

    log << setw(                          dis                          - disStepString)             << ::std::right << stepString(); 
    log << setw(disStepString           + dis1 + pP.format().tabSize() - disNormOfGradientString)   << ::std::right << normOfGradientString();
    log << setw(disNormOfGradientString + dis2 + pP.format().tabSize() - disDeltaString)            << ::std::right << deltaString();
    log << setw(disDeltaString          + dis3 + pP.format().tabSize() - disStepSizeString)         << ::std::right << stepSizeString();
    log << setw(disStepSizeString       + dis4 + pP.format().tabSize() - disErrorString)            << ::std::right << errorString();
    log << setw(disErrorString          + dis5 + pP.format().tabSize() - disRelFunctionValueString) << ::std::right << relFunctionValueString();
    log << endl;

    const LongInt size = numberOfEntrys();
    
    //Ausgabe der Daten
    
    for(LongInt i=0; i<size; i++)
     {
        const DKTSDData data = (*this).getEntryAt(i);
	
	       const LongInt  st  = data.step();
        const LongReal er  = data.error();
        const LongReal nOG = data.normOfGradient();
        const LongReal del = data.delta();
        const LongReal sS  = data.stepSize();
								const LongReal rFV = data.relFunctionValue();
	
	       IString zeros;
	       zeros.setZeros(st,4);

	       log << setw(                          dis                          - disStepNumber)            << zeros;
	       log << setprecision(prec) << scientific << uppercase;   //  << Genauigkeit << Zehnerpotenz << Grossbuchstaben
        log << setw(disStepNumber           + dis1 + pP.format().tabSize() - disNormOfGradientNumber)   << nOG;   
        log << setw(disNormOfGradientNumber + dis2 + pP.format().tabSize() - disDeltaNumber)            << del;
        log << setw(disDeltaNumber          + dis3 + pP.format().tabSize() - disStepSizeNumber)         << sS;
        log << setw(disStepSizeNumber       + dis4 + pP.format().tabSize() - disErrorNumber)            << er;
								log << setw(disErrorNumber          + dis5 + pP.format().tabSize() - disRelFunctionValueNumber) << rFV; 
	       log << resetiosflags( ::std::ios::scientific );
	       log << endl;
     }
    
   return value;
 }    	


bool DKTSDDataBlock::centeringEntrys(const LongInt& lengthOfNumber, const LongInt& lengthString, LongInt& disString, LongInt& disNumber)const
 {
    bool value = false;
    
    if(lengthOfNumber > lengthString)
     {
        disString = (lengthOfNumber - lengthString) / 2;
        disNumber = 0;
     }
    else
     {
        disNumber = (lengthString - lengthOfNumber) / 2;
        disString = 0;
     } 
     
   return value;    
  }     
     
     
bool DKTSDDataBlock::generatedPlotFiles(const IString& blockIndex, const ProtocolProperties& pP) const
 {
    const bool value = (*this).plot2File(blockIndex, pP);
    
   return value;
 }
 

bool DKTSDDataBlock::plot2File(const IString& blockIndex, const ProtocolProperties& pP)const
 {
    bool value;
    
    //Erstellen der Dateinamen
    const IString erString  = errorString();
    const IString nogString = normOfGradientString();
    const IString delString = deltaString();
    const IString sSString  = stepSizeString();
				const IString rFVString = relFunctionValueString();
    
    const IString dateiEndung = ".dat";
    
    IString fileNameError            = fileName(erString,  blockIndex, dateiEndung);
    IString fileNameNormOfGradient   = fileName(nogString, blockIndex, dateiEndung);
    IString fileNameDelta            = fileName(delString, blockIndex, dateiEndung);
    IString fileNameStepSize         = fileName(sSString,  blockIndex, dateiEndung);
				IString fileNameRelFunctionValue = fileName(rFVString,  blockIndex, dateiEndung);
    
    //Schreiben der Daten ins File
    ofstream outer  (fileNameError);
    ofstream outnog (fileNameNormOfGradient);
    ofstream outdel (fileNameDelta);
    ofstream outsS  (fileNameStepSize);
				ofstream outrFV (fileNameRelFunctionValue);

    const LongInt size =    numberOfEntrys();
    const LongInt prec = pP.format().precision();

    for(LongInt i=0; i<size; i++)
     {
        const DKTSDData data =(*this).getEntryAt(i);
	
	       const LongInt  st  = data.step();
        const LongReal er  = data.error();
	       const LongReal nog = data.normOfGradient();
       	const LongReal del = data.delta();
	       const LongReal sS  = data.stepSize();
								const LongReal rFV = data.relFunctionValue();
    
	       //  << Genauigkeit << Zehnerpotenz << Grossbuchstaben
	       outer  << setprecision(prec) << scientific << uppercase << st << '\t' << er  << endl;
	       outnog << setprecision(prec) << scientific << uppercase << st << '\t' << nog << endl;        
	       outdel << setprecision(prec) << scientific << uppercase << st << '\t' << del << endl;        
	       outsS  << setprecision(prec) << scientific << uppercase << st << '\t' << sS  << endl;
								outrFV << setprecision(prec) << scientific << uppercase << st << '\t' << rFV << endl;
     }
    
    const IString name       = IString("main-") + IString(blockIndex) + IString(".plt");
    const IString nameMainPs = IString("ps-")   + IString("main-")    + IString(blockIndex) + IString(".plt");
	   const IString nameps     = IString("'")     + IString("main-")    + IString(blockIndex) + IString(".ps") + IString("'");
	
    MainFilesWrite mFW;
				MainFilesWrite mFWps;
				
    ofstream out(name);
				ofstream outps(nameMainPs);

    mFW.mainGnuPlotFileHead(out);
   
    out << " '" << fileNameError            << "', ";
    out << " '" << fileNameNormOfGradient   << "', ";
    out << " '" << fileNameDelta            << "', ";
    out << " '" << fileNameStepSize         << "', ";
				out << " '" << fileNameRelFunctionValue << "'  ";
   
    mFW.mainGnuPlotFileEnd(out);
				
				// gnuplot-main-file-fuer-PostScript
				mFWps.mainGnuPlotFileHead(outps, "postscript color", nameps);
				
				outps << " '" << fileNameError            << "', ";
    outps << " '" << fileNameNormOfGradient   << "', ";
    outps << " '" << fileNameDelta            << "', ";
    outps << " '" << fileNameStepSize         << "', ";
				outps << " '" << fileNameRelFunctionValue << "'  ";
				
    mFWps.mainGnuPlotFileEnd(outps);
				
    value = true;
   
   return value;
 }    	      


IString DKTSDDataBlock::timeFileName(const IString& name, const IString& blockIndex, const IString& dateiEndung) const
 { 
    Timer time;
    
    const IString zeitstempel(time.timeStamp());
       
    const IString value  = IString(zeitstempel) + IString("-") + IString(blockIndex) + IString("-") + IString(name) + IString(dateiEndung);
  
   return value;
 }  	    


IString DKTSDDataBlock::fileName(const IString& name, const IString& blockIndex, const IString& dateiEndung) const
 {
    const IString value  = IString(blockIndex) + IString("-") + IString(name) + IString(dateiEndung);
    
   return value;
 }
    

IString DKTSDDataBlock::generatedTeXFiles(const IString& blockIndex, const ProtocolProperties& pP) const
 {
    //Erstellen des Dateinamens
    const IString dateiEndung = ".tex";
    
    IString file = fileName("DKTSDDataBlock", blockIndex, dateiEndung);
    
    //Schreiben der Daten ins File
    ofstream out(file);
    
    const LongInt size      =    numberOfEntrys();
    const LongInt prec      = pP.format().precision();
    const LongInt tableSize = pP.format().maxTableSize();
    const IString position  = pP.format().tablePosition();
				
    LongInt rest = size%tableSize;
    LongInt numberOfTables;
    
    if(rest != 0)
     {
        numberOfTables = size/tableSize + 1;
     }
    else
     {
        numberOfTables = size/tableSize;
     }

    LongInt lineNumber     = 0;
    
    for(LongInt i=0; i<numberOfTables; i++)
     {
        //TeX - Tabellenanfang
    
        const IString hline = IString("\\hline");
    
        const IString erTeXString  = errorTeXString();
        const IString stTeXString  = stepTeXString();
        const IString nogTeXString = normOfGradientTeXString();
        const IString delTeXString = deltaTeXString();
        const IString sSTeXString  = stepSizeTeXString();
								const IString rFVTeXString = relFunctionValueTeXString(); 
    
        out << "\\begin{table}[" << position << "]" << endl;
        out << "\\begin{center}"                    << endl;
        out << "\\begin{tabular}{cccccc}"           << endl;
        
								out << " $"  << stTeXString  << "$ " << " & ";
        out << " $"  << nogTeXString << "$ " << " & ";
        out << " $"  << delTeXString << "$ " << " & ";
        out << " $"  << sSTeXString  << "$ " << " & ";
        out << " $"  << erTeXString  << "$ " << " & ";
								out << " $"  << rFVTeXString << "$ ";
        out << "\\"  << "\\"         << endl;
        out << hline << endl;
     
        //TeX - Tabellenelemente
        out << setprecision(prec) << uppercase;    //  << Genauigkeit << Grossbuchstaben
        
	       LongInt length = MIN(lineNumber + tableSize, size); 
	
        for(LongInt j=lineNumber; j<length; j++)
         {
            const DKTSDData data = (*this).getEntryAt(j);
        
	           const LongInt  st  = data.step();
            const LongReal er  = data.error();
	           const LongReal nog = data.normOfGradient();
            const LongReal del = data.delta();
            const LongReal sS  = data.stepSize();
												const LongReal rFV = data.relFunctionValue();
	
            out << st  << " & "; 
	           out << scientific;     // << Zehnerpotenz
												out << nog << " & ";
												out << del << " & ";
	           out << sS  << " & ";
	           out << er  << " & ";
												out << rFV << "\\\\";
												out << resetiosflags( ::std::ios::scientific );
												out << endl;
         }
     
        out << hline << endl;
        //TeX - Tabellenabschluss 	
    
        out << "\\end{tabular}" << "\\caption{" << captionString();
								
								if(numberOfTables > 1)
								 {
								    out << ",\\hspace{15pt}$" << i + 1 << "$" << "$\\backslash$" << "$" << numberOfTables << "$" << "}" << endl;
									}
								else
								 {
									   out << "}" << endl;
									}
																
        out << "\\end{center}"  << endl;
        out << "\\end{table}";
	 
        lineNumber += tableSize;
     }	 
    
   return file;
 }


bool DKTSDDataBlock::generatedAll(const IString& blockIndex, const ProtocolProperties& pP) const
 {
    bool value = false;
    
    generatedPlotFiles (blockIndex, pP);
    generatedTeXFiles  (blockIndex, pP);
    
   return value;
 }   
    
