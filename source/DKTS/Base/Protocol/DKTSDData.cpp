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

// DKTSDData.cpp: Implementierung der Klasse DKTSDData.
//
//////////////////////////////////////////////////////////////////////

#include "DKTSDData.hpp"
#include <typeinfo>

DKTSDData::DKTSDData (const LongInt&  st,  const LongReal& er, 
                      const LongReal& nOG, const LongReal& del, const LongReal& sS, const LongReal& rFV)
 :SimpleIterationData(st, er)
 {
    (*this)
     .setNormOfGradient(nOG)
     .setDelta(del)
     .setStepSize(sS)
					,setRelFunctionValue(rFV)
     .setDefaultParameter()
    ;
 }
 

DKTSDData::DKTSDData(const DKTSDData& dD)
 :SimpleIterationData(dD.step(), dD.error())
 {     
    (*this) = dD;     
 }   
     
DKTSDData::~DKTSDData()
 {
 
 }

  
DKTSDData DKTSDData::setDefaultParameter()
 {
    (*this)
    ;
    
   return (*this);  
 }
 

DKTSDData& DKTSDData::operator = (const DKTSDData& data)
 {
    (*this)
     
     .setNormOfGradient   (data.normOfGradient())
     .setDelta            (data.delta())
     .setStepSize         (data.stepSize())
     .setRelFunctionValue (data.relFunctionValue())
					.setStep             (data.step())
     .setError            (data.error())
    ;
   
   return (*this); 
 }

 
bool DKTSDData::operator == (const ProtocolData& D)
 {
    bool value = false;
    
    if(isInstanceOf(&D, DKTSDData) == true)
     {
        DKTSDData& data = (DKTSDData&) D;
        
        if(step() == data.step() && fabs(error()-data.error())<EPS_NULL)
         {
            value = true;
         } 
     }
    
   return value;
 }        


ostream& operator << (ostream& s, const DKTSDData& data)
 {
    const LongInt  st  = data.step();
    const LongReal er  = data.error();
    const LongReal nOG = data.normOfGradient();
    const LongReal del = data.delta();
    const LongReal sS  = data.stepSize();
				const LongReal rFV = data.relFunctionValue();
    
				IString zeros;
    
    zeros.setZeros(st, 4);
				
    s << zeros << '\t';
				s << scientific << uppercase;  // << Zehnerpotenz << Grossbuchstaben
				s << er << '\t' << nOG << '\t' << del << '\t' << sS << '\t' << rFV << endl;
				s << resetiosflags( ::std::ios::scientific );
    
   return s;
 }   


istream& operator >> (istream& s, DKTSDData& data)
 {
	   LongInt  st  = 0;
    LongReal er  = 0;
    LongReal nOG = 0;
    LongReal del = 0;
    LongReal sS  = 0;
				LongReal rFV = 0;
				
				s >> st;
				s >> er;
				s >> nOG;
				s >> del;
				s >> sS;
				s >> rFV;
				
				data.setStep(st);
				data.setError(er);
				data.setNormOfGradient(nOG);
				data.setDelta(del);
				data.setStepSize(sS);
				data.setRelFunctionValue(rFV);
				
			return s;	
	}			
				 

bool DKTSDData::plot(ostream& log, const ProtocolProperties& pP) const
 {
    bool value = false;
    
    const LongInt  st  = step();
    const LongReal er  = error();
    const LongReal nOG = normOfGradient();
    const LongReal del = delta();
    const LongReal sS  = stepSize();
				const LongReal rFV = relFunctionValue();
    
    IString zeros;
    
    zeros.setZeros(st, 4);
    
    log << setw(4)                                                  << zeros;
    log << setw(pP.format().precision()+ pP.format().tabSize() + 6) << nOG;   
    log << setw(pP.format().precision()+ pP.format().tabSize() + 6) << del;
    log << setw(pP.format().precision()+ pP.format().tabSize() + 6) << sS;
    log << setw(pP.format().precision()+ pP.format().tabSize() + 6) << er;
				log << setw(pP.format().precision()+ pP.format().tabSize() + 6) << rFV;
    log << endl;
    
   return value;
 }
