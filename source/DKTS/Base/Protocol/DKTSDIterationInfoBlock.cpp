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

// DKTSDIterationInfoBlock.cpp: Implementierung der Klasse DKTSDIterationInfoBlock.
//
//////////////////////////////////////////////////////////////////////

#include "DKTSDIterationInfoBlock.hpp"


DKTSDIterationInfoBlock::DKTSDIterationInfoBlock ()
 {
    (*this)
     .setDefaultParameter()
    ;
 }
 
 
DKTSDIterationInfoBlock::DKTSDIterationInfoBlock(const DKTSDIterationInfoBlock& data)
 {
    (*this) = data;
 }
 
        
DKTSDIterationInfoBlock::~DKTSDIterationInfoBlock()
 {
 
 }


DKTSDIterationInfoBlock DKTSDIterationInfoBlock::setDefaultParameter()
 {
    (*this)
     .setStepString                   ("Step")
     .setStartErrorString             ("StartError")
     .setEndErrorString               ("EndError")
     .setRelativeDifferenzString      ("relDiff")
					.setGradientString               ("relGradient")
     .setNumberOfNewtonStepsString    ("Steps")
     .setCalculationTimeString        ("Time")
     .setStepTeXString                ("kRank")
     .setStartErrorTeXString          ("\\frac{\\|A-X_0^{kR}\\|}{\\|A\\|}")
     .setEndErrorTeXString            ("\\frac{\\|A-X^{kR}\\|}{\\|A\\|}")
     .setRelativeDifferenzTeXString   ("relDiff.[\\%]")
					.setGradientTeXString            ("\\frac{\\|f'(x^{KR})\\|}{\\|A\\|}")
     .setNumberOfNewtonStepsTeXString ("Steps")
     .setCalculationTimeTeXString     ("time[sec.]")
    ;
    
   return (*this);
 }


DKTSDIterationInfoBlock& DKTSDIterationInfoBlock::operator = (const DKTSDIterationInfoBlock& data)
 {
    (*this)
     .setStepString                   (data.stepString())
     .setStartErrorString             (data.startErrorString())
     .setEndErrorString               (data.endErrorString())
     .setNumberOfNewtonStepsString    (data.numberOfNewtonStepsString())
     .setCalculationTimeString        (data.calculationTimeString())
     .setRelativeDifferenzString      (data.relativeDifferenzString())
					.setGradientString               (data.gradientString())
     .setStepTeXString                (data.stepTeXString())
     .setStartErrorTeXString          (data.startErrorTeXString())
     .setEndErrorTeXString            (data.endErrorTeXString())
     .setNumberOfNewtonStepsTeXString (data.numberOfNewtonStepsTeXString())
     .setCalculationTimeTeXString     (data.calculationTimeTeXString())
     .setRelativeDifferenzTeXString   (data.relativeDifferenzTeXString())
					.setGradientTeXString            (data.gradientTeXString())
     .setCaptionString                (data.captionString())
    ;
    
    const LongInt size = data.numberOfEntrys();
    
    for(LongInt i=0; i<size; i++)
     {
        (*this).addEntry(data.getEntryAt(i));
     }
     
   return (*this);
 }


ostream& operator << (ostream& s, DKTSDIterationInfoBlock& data)
 {
	   const IString stepString = data.stepString();
				const IString serString  = data.startErrorString();
				const IString errString  = data.endErrorString(); 
				const IString raString   = data.relativeDifferenzString();
				const IString grString   = data.gradientString(); 
    const IString nonsString = data.numberOfNewtonStepsString(); 
    const IString ctString   = data.calculationTimeString();
    
 
	   s << stepString << '\t' << serString << '\t' << errString  << '\t';
				s << raString   << '\t' << grString  << '\t' << nonsString << '\t';
				s << ctString   << '\t' << endl;
				
    const LongInt size = data.numberOfEntrys();
	
    for(LongInt i=0; i<size; i++)
				 {
        const DKTSDIterationInfo infoData = data.getEntryAt(i);
								
								s << infoData;
					}
					
			return s;
	}						
	

bool DKTSDIterationInfoBlock::plot(ostream& log, const ProtocolProperties& pP)const
 {
    bool value = false;
    
    plotEntry(log, pP);
    
   return value;
 }
 
 
bool DKTSDIterationInfoBlock::plotEntry(ostream& log, const ProtocolProperties& pP)const
 {
    bool value = false;
    
    //Bestimmung der Laenge der Strings
    
    const LongInt prec = pP.format().precision();
    
    const LongInt lengthStepString                = stepString()                .length();
    const LongInt lengthStartErrorString          = startErrorString()          .length();
    const LongInt lengthEndErrorString            = endErrorString()            .length();
				const LongInt lengthRelativeDifferenzString   = relativeDifferenzString()   .length();
				const LongInt lengthGradientString            = gradientString()            .length();
    const LongInt lengthNumberOfNewtonStepsString = numberOfNewtonStepsString() .length();
    const LongInt lengthCalculationTimeString     = calculationTimeString()     .length();
    
    
    const LongInt lengthOfNumber = (pP.format().precision() + 6);

    LongInt disStepString,                disStepNumber;
    LongInt disStartErrorString,          disStartErrorNumber;
    LongInt disEndErrorString,            disEndErrorNumber;
				LongInt disRelativeDifferenzString,   disRelativeDifferenzNumber;
				LongInt disGradientString,            disGradientNumber;
    LongInt disNumberOfNewtonStepsString, disNumberOfNewtonStepsNumber;
    LongInt disCalculationTimeString,     disCalculationTimeNumber;
    
             	      
    //Berechnung der Werte zur Zentrierung

    centeringEntrys(4,              lengthStepString,                disStepString,                disStepNumber);
    centeringEntrys(lengthOfNumber, lengthStartErrorString,          disStartErrorString,          disStartErrorNumber);
    centeringEntrys(lengthOfNumber, lengthEndErrorString,            disEndErrorString,            disEndErrorNumber);
				centeringEntrys(lengthOfNumber, lengthRelativeDifferenzString,   disRelativeDifferenzString,   disRelativeDifferenzNumber);
				centeringEntrys(lengthOfNumber, lengthGradientString,            disGradientString,            disGradientNumber);
    centeringEntrys(lengthOfNumber, lengthNumberOfNewtonStepsString, disNumberOfNewtonStepsString, disNumberOfNewtonStepsNumber);
    centeringEntrys(lengthOfNumber, lengthCalculationTimeString,     disCalculationTimeString,     disCalculationTimeNumber);
    
    
    const LongInt dis  = MAX(4,              lengthStepString);
    const LongInt dis1 = MAX(lengthOfNumber, lengthStartErrorString); 
    const LongInt dis2 = MAX(lengthOfNumber, lengthEndErrorString);
				const LongInt dis3 = MAX(lengthOfNumber, lengthRelativeDifferenzString);
				const LongInt dis4 = MAX(lengthOfNumber, lengthGradientString);
    const LongInt dis5 = MAX(lengthOfNumber, lengthNumberOfNewtonStepsString);
    const LongInt dis6 = MAX(lengthOfNumber, lengthCalculationTimeString);
    
    
    //Ausgabe der Strings als Tabellenkopf
   
    log << setw(                               dis                          - disStepString)                << stepString(); 
    log << setw(disStepString                + dis1 + pP.format().tabSize() - disStartErrorString)          << startErrorString();
    log << setw(disStartErrorString          + dis2 + pP.format().tabSize() - disEndErrorString)            << endErrorString();
				log << setw(disEndErrorString            + dis3 + pP.format().tabSize() - disRelativeDifferenzString)   << relativeDifferenzString();
    log << setw(disRelativeDifferenzString   + dis4 + pP.format().tabSize() - disGradientString)            << gradientString();
				log << setw(disGradientString            + dis5 + pP.format().tabSize() - disNumberOfNewtonStepsString) << numberOfNewtonStepsString();
    log << setw(disNumberOfNewtonStepsString + dis6 + pP.format().tabSize() - disCalculationTimeString)     << calculationTimeString();
    log << endl;
    
    const LongInt size = numberOfEntrys();
    
    //Ausgabe der Daten
    
    for(LongInt i=0; i<size; i++)
     {
        const DKTSDIterationInfo data = (*this).getEntryAt(i);
	
	       const LongInt  st   = data.step();
	       const LongReal ser  = data.startError();
        const LongReal eer  = data.error();
        const LongReal ra   = data.relativeDifferenz();
								const LongReal gr   = data.gradient();
        const LongReal nons = data.numberOfNewtonSteps();
        const LongReal ct   = data.calculationTime();
	
	       IString zeros;
	       zeros.setZeros(st,4);
	
	       log << setw(                               dis                          - disStepNumber)                << ::std::right << zeros; 
	       log << setprecision(prec) << scientific;
	       log << setw(disStepNumber                + dis1 + pP.format().tabSize() - disStartErrorNumber)          << ::std::right << ser;
        log << setw(disStartErrorNumber          + dis2 + pP.format().tabSize() - disEndErrorNumber)            << ::std::right << eer;
        log << setw(disEndErrorNumber            + dis3 + pP.format().tabSize() - disRelativeDifferenzNumber)   << ::std::right << ra;
								log << setw(disRelativeDifferenzNumber   + dis4 + pP.format().tabSize() - disGradientNumber)            << ::std::right << gr;
	       log << resetiosflags( ::std::ios::scientific );
								log << setw(disGradientNumber            + dis3 + pP.format().tabSize() - disNumberOfNewtonStepsNumber) << ::std::right << nons;
	       log << setw(disNumberOfNewtonStepsNumber + dis4 + pP.format().tabSize() - disCalculationTimeNumber)     << ::std::right << ct;
        log << endl;
     }
      
   return value;
 }
 
      
bool DKTSDIterationInfoBlock::centeringEntrys(const LongInt& lengthOfNumber, const LongInt& lengthString, LongInt& disString, LongInt& disNumber)const
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
	
	
bool DKTSDIterationInfoBlock::generatedPlotFiles(const IString& blockIndex, const ProtocolProperties& pP) const
 {
    const bool value = (*this).plot2File(blockIndex, pP);
    
   return value;
 }
 
 
bool DKTSDIterationInfoBlock::plot2File(const IString& blockIndex, const ProtocolProperties& pP)const
 {
    bool value = false;
    
    //Erstellen der Dateinamen
    const IString serString  = startErrorString();
    const IString errString  = endErrorString();
    const IString raString   = relativeDifferenzString();
				const IString grString   = gradientString();
    const IString nonsString = numberOfNewtonStepsString();
    const IString ctString   = calculationTimeString();
    
    const IString dateiEndung = ".dat";
    
    IString fileNameStartError = fileName(serString,  blockIndex, dateiEndung);    
    IString fileNameEndError   = fileName(errString,  blockIndex, dateiEndung);
    IString fileNameRA         = fileName(raString,   blockIndex, dateiEndung);
				IString fileNameGr         = fileName(grString,   blockIndex, dateiEndung);
    IString fileNameNONS       = fileName(nonsString, blockIndex, dateiEndung);
    IString fileNameCT         = fileName(ctString,   blockIndex, dateiEndung);
    
    //Schreiben der Daten ins File
    ofstream outser  (fileNameStartError);
    ofstream outeer  (fileNameEndError);
    ofstream outra   (fileNameRA);
				ofstream outgr   (fileNameGr);
    ofstream outnons (fileNameNONS);
    ofstream outct   (fileNameCT);

    const LongInt size        =    numberOfEntrys();
    const LongInt prec        = pP.format().precision();
				const IString projektName = pP.projektName();

    for(LongInt i=0; i<size; i++)
     {
        const DKTSDIterationInfo data =(*this).getEntryAt(i);
	
	       const LongInt  st   = data.step();
	       const LongReal ser  = data.startError();
        const LongReal eer  = data.error();
        const LongReal ra   = data.relativeDifferenz();
								const LongReal gr   = data.gradient();
        const LongReal nons = data.numberOfNewtonSteps();
        const LongReal ct   = data.calculationTime();
    
		      //  << Genauigkeit << Zehnerpotenz << Grossbuchstaben
		      outser  << setprecision(prec) <<                            st << '\t' << ser  << endl;
		      outeer  << setprecision(prec) << scientific << uppercase << st << '\t' << eer  << endl;
								outra   << setprecision(prec) << scientific << uppercase << st << '\t' << ra   << endl;
								outgr   << setprecision(prec) << scientific << uppercase << st << '\t' << gr   << endl;
		      outnons << setprecision(prec) <<                            st << '\t' << nons << endl;
		      outct   << setprecision(prec) <<                            st << '\t' << ct   << endl;
     }
     
    const IString name       = IString("main-") + IString(blockIndex) + IString(".plt");
	   const IString nameMainPs = IString("ps-")   + IString("main-")    + IString(blockIndex) + IString(".plt");
	   const IString nameps     = IString("'")     + IString("main-")    + IString(blockIndex) + IString(".ps")  + IString("'");
    
				ofstream out(name);
				ofstream outps(nameMainPs);

    MainFilesWrite mFW;
				MainFilesWrite mFWps;
    
				mFW.mainGnuPlotFileHead(out);
    
    out << " '" << fileNameStartError << "',";
    out << " '" << fileNameEndError   << "',";
				out << " '" << fileNameRA         << "',";
				out << " '" << fileNameGr         << "',";
    out << " '" << fileNameNONS       << "',";
    out << " '" << fileNameCT         << "' ";
      
    
    mFW.mainGnuPlotFileEnd(out);
    
				// gnuplot-main-file-fuer-PostScript
				mFWps.mainGnuPlotFileHead(outps, "postscript color", nameps);
				
				out << " '" << fileNameStartError << "',";
    out << " '" << fileNameEndError   << "',";
				out << " '" << fileNameRA         << "',";
				out << " '" << fileNameGr         << "',";
    out << " '" << fileNameNONS       << "',";
    out << " '" << fileNameCT         << "' ";
				
    mFWps.mainGnuPlotFileEnd(outps);
				
    value = true;
    
   return value;
 }    	      


IString DKTSDIterationInfoBlock::timeFileName(const IString& name, const IString& blockIndex, const IString& dateiEndung) const
 { 
    Timer time;
    
    const IString zeitstempel(time.timeStamp());
       
    const IString value = IString(zeitstempel) + IString("-") + IString(blockIndex) + IString("-") + IString(name) + IString(dateiEndung);
  
   return value;
 }  	    


IString DKTSDIterationInfoBlock::fileName(const IString& name, const IString& blockIndex, const IString& dateiEndung) const
 {
    const IString value = IString(blockIndex) + IString("-") + IString(name) + IString(dateiEndung);
    
   return value;
 }
 

IString DKTSDIterationInfoBlock::generatedTeXFiles(const IString& blockIndex, const ProtocolProperties& pP) const
 {
    //Erstellen des Dateinamens
    const IString dateiEndung = ".tex";
    
    IString file = fileName("DKTSDIterationInfoBlock", blockIndex, dateiEndung);
    
    //Schreiben der Daten ins File
    ofstream out(file);
    
    const LongInt size      = numberOfEntrys();
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
     	
    LongInt lineNumber = 0;

    for(LongInt i=0; i<numberOfTables; i++)
     {
        //TeX - Tabellenanfang
    
        const IString hline = IString("\\hline");
    
        const IString stTeXString   = stepTeXString();
        const IString serTeXString  = startErrorTeXString();
        const IString errTeXString  = endErrorTeXString();
        const IString raTeXString   = relativeDifferenzTeXString();
								const IString grTeXString   = gradientTeXString(); 
        const IString nonsTeXString = numberOfNewtonStepsTeXString();
        const IString ctTeXString   = calculationTimeTeXString();
    
        out << "\\begin{table}[" << position << "]" << endl; 
        out << "\\begin{center}"                    << endl;
        out << "\\begin{tabular}{ccccccc}"          << endl;
        //out << hline << endl;
        out << " $"  << stTeXString   << "$ " << " & ";
        out << " $"  << serTeXString  << "$ " << " & ";
        out << " $"  << errTeXString  << "$ " << " & ";
								out << " $"  << raTeXString   << "$ " << " & ";
								out << " $"  << grTeXString   << "$ " << " & ";
        out << " $"  << nonsTeXString << "$ " << " & ";
        out << " $"  << ctTeXString   << "$ ";
        out << "\\"  << "\\"          << endl;
        out << hline << endl;
     
        //TeX - Tabellenelemente

        //  << Genauigkeit << Zehnerpotenz << Grossbuchstaben
        out << setprecision(prec) << uppercase;
    
	       LongInt length = MIN(lineNumber + tableSize, size);
    
        for(LongInt j=lineNumber; j<length; j++)
         {
            const DKTSDIterationInfo data = (*this).getEntryAt(j);
        
            const LongInt  st   = data.step();
            const LongReal ser  = data.startError();
            const LongReal eer  = data.error();
            const LongReal ra   = data.relativeDifferenz();
												const LongReal gr   = data.gradient();
            const LongReal nons = data.numberOfNewtonSteps();
            const LongReal ct   = data.calculationTime();
	
            out <<               st   << " & ";
	           out << scientific << ser  << " & ";
	           out <<               eer  << " & ";
												out <<               ra   << " & ";
												out <<               gr   << " & ";
	           out << resetiosflags(::std::ios::scientific);
	           out <<               nons << " & ";
	           out <<               ct;
												out << "\\\\";
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
        out << "\\end{table}"   << endl;
	
        lineNumber += tableSize;
     }	
    
   return file;
 }


bool DKTSDIterationInfoBlock::generatedAll(const IString& blockIndex, const ProtocolProperties& pP) const
 {
    bool value = false;
    
    generatedPlotFiles (blockIndex, pP);
    generatedTeXFiles  (blockIndex, pP);
    
   return value;
 }   
    
      
