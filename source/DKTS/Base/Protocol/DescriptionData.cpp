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

// DescriptionData.cpp: Implementierung der Klasse DescriptionData.
//
//////////////////////////////////////////////////////////////////////

#include "DescriptionData.hpp"


DescriptionData::DescriptionData()
 {
    
 }
 
     
DescriptionData::~DescriptionData()
 {
 
 }


ostream& operator << (ostream& s, const DescriptionData& data)
 {
	   const LongInt size = data.numberOfStrings();
        
    for(LongInt i=0; i<size; i++)
     { 
        const IString& string = data.getStringAt(i);
        
        s << string << endl;
     }
    
   return s;
 }


bool DescriptionData::plot(ostream& log, const ProtocolProperties& pP) const
 { 
    bool value = true;
    
    const LongInt size = numberOfStrings();
        
    for(LongInt i=0; i<size; i++)
     { 
        const IString& string = (*this).getStringAt(i);
        
        log << string << endl;
     }
    
   return value;
   
 }
 
    
IString DescriptionData::generatedTeXFiles(const IString& blockIndex, const ProtocolProperties& pP) const
 {
    //Erstellen des Dateinamens
    const IString dateiEndung = ".tex";
    
    IString file = fileName("DescriptionData", blockIndex, dateiEndung);
    
    //Schreiben der Daten ins File
    ofstream out(file);
    
    const LongInt size =    numberOfStrings();
    const LongInt prec = pP.format().precision();
    
    out << setprecision(prec) << scientific << uppercase;
    
    for(LongInt i=0; i<size; i++)
     {
	       const IString& string = (*this).getStringAt(i);
	
	       if((string != "") && (string != " "))
         {
            out << string << endl;
	        }
       	else
	        {
	           out << "\\bigskip\\" << "\\" << endl;  //Leerzeile im LaTeX
	        }      
     }
    
   return file;
 }


IString DescriptionData::fileName(const IString& name, const IString& blockIndex, const IString& dateiEndung) const
 {
    const IString value = IString(blockIndex) + IString("-") + IString(name) + IString(dateiEndung);
    
   return value;
 }
  

IString DescriptionData::timeFileName(const IString& name, const IString& blockIndex, const IString& dateiEndung) const
 {
    Timer time;
    
    const IString zeitstempel(time.timeStamp());
       
    const IString value = IString(zeitstempel) + IString("-") + IString(blockIndex) + IString("-") + IString(name) +	IString(dateiEndung);
  
   return value;
 } 


bool DescriptionData::generatedAll(const IString& blockIndex, const ProtocolProperties& pP) const
 {
    bool value = false;
    
    generatedPlotFiles (blockIndex, pP);
    generatedTeXFiles  (blockIndex, pP);
    
   return value;
 }


IString DescriptionData::addMathString(const IString& s)
 {
    const IString st = "$";
				
    IString mathString = st + s + st;
    
    (*this).addString(mathString);
    
   return mathString;
 }
    
    
