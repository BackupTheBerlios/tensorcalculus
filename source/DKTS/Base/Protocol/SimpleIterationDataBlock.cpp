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

// SimpleIterationDataBlock.cpp: Implementierung der Klasse SimpleIterationDataBlock.
//
//////////////////////////////////////////////////////////////////////

#include "SimpleIterationDataBlock.hpp"


SimpleIterationDataBlock::SimpleIterationDataBlock ()
 {
    (*this)
     .setDefaultParameter()
    ; 
 }


SimpleIterationDataBlock::SimpleIterationDataBlock(const SimpleIterationDataBlock& data)
 {
    (*this) = data;
 }   

     
SimpleIterationDataBlock::~SimpleIterationDataBlock()
 {
 
 }


SimpleIterationDataBlock SimpleIterationDataBlock::setDefaultParameter()
 {
    (*this)
     .setStepString     ("Step")
     .setErrorString    ("EndError")
     .setStepTeXString  ("Step")
     .setErrorTeXString ("\\frac{\\|A-X_{itr}\\|}{\\|A\\|}")
    ;
    
   return (*this); 
 }
 

SimpleIterationDataBlock& SimpleIterationDataBlock::operator = (const SimpleIterationDataBlock& data)
 {
    (*this)
     .setStepString     (data.stepString())
     .setErrorString    (data.errorString())
     .setStepTeXString  (data.stepTeXString())
     .setErrorTeXString (data.errorTeXString())
    ;
  
    const LongInt size = data.numberOfEntrys();
  
    for(LongInt i=0; i<size; i++)
     {
	       (*this).addEntry(data.getEntryAt(i));
     }
   
   return (*this);    
 }       


ostream& operator << (ostream& s, const SimpleIterationDataBlock& sD)
 {
	   const IString stepString  = sD.stepString();
				const IString errorString = sD.errorString();
				
				s << stepString << '\t' << errorString << endl;
				
				const LongInt size = sD.numberOfEntrys();

    for(LongInt i=0; i<size; i++)
     {
        const SimpleIterationData data = sD.getEntryAt(i);
								
								s << data;
					}
					
			return s;
	}		
								

bool SimpleIterationDataBlock::plot(ostream& log, const ProtocolProperties& pP) const
{
   bool value = true;
   
   plotEntry(log, pP);
   
  return value;
}      


bool SimpleIterationDataBlock::plotEntry(ostream& log, const ProtocolProperties& pP) const
 {
    bool value = true;
     
    const LongInt lengthStepString  = stepString()  .length();
    const LongInt lengthErrorString = errorString() .length();
   
    const LongInt lengthOfNumber = (pP.format().precision() + 6);
   
    LongInt disStepString;
    LongInt disStepNumber;
    LongInt disErrorString;
    LongInt disErrorNumber; 
   
    centeringEntrys(4,              lengthStepString,  disStepString,  disStepNumber);
    centeringEntrys(lengthOfNumber, lengthErrorString, disErrorString, disErrorNumber);
     
    const LongInt disStep  = MAX(4,              lengthStepString);
    const LongInt disError = MAX(lengthOfNumber, lengthErrorString);  	
   
    log << setw(disStep                                          - disStepString)  << stepString  ();
    log << setw(disStepString + disError + pP.format().tabSize() - disErrorString) << errorString ();
    log << endl; 
   
    const LongInt size = numberOfEntrys();

    for(LongInt i=0; i<size; i++)
     {
        const SimpleIterationData data = (*this).getEntryAt(i);
       
        const LongInt  st = data.step();
        const LongReal er = data.error();
        
	       IString zerosStep;
	       zerosStep.setZeros(st,4);
       
        log << setw(disStep                                          - disStepNumber)  << zerosStep;
        log << setw(disStepNumber + disError + pP.format().tabSize() - disErrorNumber) << er;
        log << endl;
     }
   return value;
 }   
  
  
bool SimpleIterationDataBlock::centeringEntrys(const LongInt& lengthOfNumber, const LongInt& lengthString, LongInt& disString, LongInt& disNumber)const
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
   

bool SimpleIterationDataBlock::generatedPlotFiles(const IString& blockIndex, const ProtocolProperties& pP) const
 {
    const bool value = plot2File(blockIndex, pP);
      
   return value;
 }
 

bool SimpleIterationDataBlock::plot2File(const IString& blockIndex, const ProtocolProperties& pP) const
 {
    bool value = false;
    
    const LongInt prec = pP.format().precision();
				
    //Erstellen des Dateinamens
    const IString erString    = errorString();
    const IString dateiEndung = ".dat"; 
    
    IString fileNameError = fileName(erString, blockIndex, dateiEndung);	

    //Schreiben der Daten ins File	     
    ofstream outA(fileNameError);
    
    outA << setprecision(prec) << uppercase;       // << Genauigkeit << Grossbuchstaben				
    
				const LongInt size =    numberOfEntrys();
				
    for(LongInt i=0; i < size; i++)
     {
        const SimpleIterationData data = (*this).getEntryAt(i);
       
        const LongInt  st  = data.step();
        const LongReal er  = data.error(); 
	       
       
	       outA << st         << '\t';
								outA << scientific;                        // << Zehnerpotenz 
								outA << er         << endl;
	       outA << resetiosflags( ::std::ios::scientific );
     }
    
    MainFilesWrite mfw;
   
    IString name  = IString("main-") + IString(blockIndex) + IString(".plt");
	
    ofstream out(name);

    mfw.mainGnuPlotFileHead(out);
   
    out << " '" << fileNameError << "'";
    
    mfw.mainGnuPlotFileEnd(out);
    
				// gnuplot-main-file-fuer-PostScript
				MainFilesWrite mfwps;
   
    const IString nameMainPs = IString("ps")    + IString("main-")    + IString(blockIndex) + IString(".plt");
	   const IString nameps     = IString("'")     + IString("main-")    + IString(blockIndex) + IString(".ps") + IString("'");	
    
				ofstream outps(nameMainPs);
				
    mfwps.mainGnuPlotFileHead(outps, "postscript color", nameps);
   
    outps << " '" << fileNameError << "'";
    
    mfwps.mainGnuPlotFileEnd(outps);
				
				
    value = true;

   return value; 
 }


IString SimpleIterationDataBlock::fileName(const IString& name, const IString& blockIndex, const IString& dateiEndung) const
 {
    const IString value  = IString(blockIndex) + IString("-") + IString(name) + IString(dateiEndung);
    
   return value;
 }
 

IString SimpleIterationDataBlock::timeFileName(const IString& name, const IString& blockIndex, const IString& dateiEndung) const
 {
    Timer time;
    
    const IString zeitstempel(time.timeStamp());
       
    const IString value = IString(zeitstempel) + IString("-") + IString(blockIndex) + IString("-") + IString(name) + IString(dateiEndung);
  
   return value;
 }  	    


IString SimpleIterationDataBlock::generatedTeXFiles(const IString& blockIndex, const ProtocolProperties& pP) const
 {
    //Erstellen des Dateinamens
    const IString dateiEndung = ".tex";
    
    IString file = fileName("SimpleIterationDataBlock", blockIndex, dateiEndung);
    
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

    LongInt lineNumber = 0;
    
    for(LongInt i=0; i<numberOfTables; i++)
     {
        //TeX - Tabellenanfang
        const IString hline = IString("\\hline");
								
        const IString erTeXString = errorTeXString();
        const IString stTeXString = stepTeXString();
    
        out << "\\begin{table}[" << position << "]" << endl;
        out << "\\begin{center}"                    << endl;
        out << "\\begin{tabular}{|c|c|}"            << endl;
        out << hline << endl;
        out << " $"  << stTeXString << "$ "  << " & ";
        out << " $"  << erTeXString << "$ "; 
        out << "\\"  << "\\"        << endl;
        out << hline << endl;
     
        //TeX - Tabellenelemente
        out << setprecision(prec) << uppercase;   //  << Genauigkeit << Grossbuchstaben
	
	       LongInt length = MIN(lineNumber + tableSize, size); 
        
	       for(LongInt j=0; j<size; j++)
         {
            const SimpleIterationData data = (*this).getEntryAt(j);
        
	           const LongInt  st  = data.step();
            const LongReal er  = data.error();
	
	           out << st << " & ";
												out << scientific;                  // << Zehnerpotenz 
												out << er << "\\" << "\\" << endl;
	           out << resetiosflags( ::std::ios::scientific );
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
    
    
bool SimpleIterationDataBlock::generatedAll(const IString& blockIndex, const ProtocolProperties& pP) const
 {
    bool value = false;
    
    generatedPlotFiles (blockIndex, pP);
    generatedTeXFiles  (blockIndex, pP);
    
   return value;
 }
    
     

    
