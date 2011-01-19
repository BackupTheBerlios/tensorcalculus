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

// DKTSDIterationInfo.cpp: Implementierung der Klasse DKTSDIterationInfo.
//
//////////////////////////////////////////////////////////////////////

#include "DKTSDIterationInfo.hpp"
#include <typeinfo>

DKTSDIterationInfo::DKTSDIterationInfo(const LongInt& st,  const LongReal& er,   const LongReal& ser, const LongReal& ra,
                                       const LongReal& gr, const LongReal& nons, const LongReal& ct)
 :SimpleIterationData(st, er)
 {
    (*this)
     .setStartError(ser)
					.setRelativeDifferenz(ra)
					.setGradient(gr)
     .setNumberOfNewtonSteps(nons)
     .setCalculationTime(ct)
     .setDefaultParameter()
    ; 
 }
    

DKTSDIterationInfo::DKTSDIterationInfo(const DKTSDIterationInfo& data)
 :SimpleIterationData(data.step(), data.error())
 {
    (*this) = data;
 }
 
         
DKTSDIterationInfo::~DKTSDIterationInfo()
 {
 
 }


DKTSDIterationInfo DKTSDIterationInfo::setDefaultParameter()
 {
    (*this)
    ;
    
   return (*this);
 }
 
 
DKTSDIterationInfo& DKTSDIterationInfo::operator = (const DKTSDIterationInfo& data)
 {
    (*this)
     .setStartError          (data.startError())
     .setNumberOfNewtonSteps (data.numberOfNewtonSteps())
     .setCalculationTime     (data.calculationTime())
     .setRelativeDifferenz   (data.relativeDifferenz())
					.setGradient            (data.gradient())
     .setStep                (data.step())
     .setError               (data.error())
    ;
    
   return (*this);
 }


bool DKTSDIterationInfo::operator == (const ProtocolData& data)
 {
    bool value = false;
    
    if(isInstanceOf(&data, DKTSDIterationInfo) == true)
     {
        DKTSDIterationInfo& dii = (DKTSDIterationInfo&) data;
	
	if(step() == dii.step() && fabs(error()-dii.error())<EPS_NULL)
         {
            value = true;
         } 
     }
    
   return value;
 }    


ostream& operator << (ostream& s, const DKTSDIterationInfo& data)
 {
    const LongInt  st   = data.step();
    const LongReal er   = data.error();
    const LongReal ser  = data.startError();
    const LongReal nons = data.numberOfNewtonSteps();
    const LongReal ct   = data.calculationTime();
    const LongReal ra   = data.relativeDifferenz();
				const LongReal gr   = data.gradient();
				
				IString zeros;
    
    zeros.setZeros(st, 4);
				 
				s << zeros         << '\t';
				s << scientific    << uppercase;
				s << ser           << '\t' << er   << '\t';
				s	<< ra            << '\t' << gr   << '\t';
				s << resetiosflags( ::std::ios::scientific );
				s << setw(7)       << nons << '\t';
				s << setw(10)      << ct;
    
   return s;
 }   


istream& operator >> (istream& s, DKTSDIterationInfo& data)
 {
	   LongInt  st   = 0;
    LongReal er   = 0;
    LongReal ser  = 0;
				LongReal ra   = 0;
				LongReal gr   = 0;
    LongReal nons = 0;
    LongReal ct   = 0;
				
				
				s >> st;
				s >> ser;
				s >> er;
				s >> ra;
				s >> gr;
				s >> nons;
				s >> ct;
				
				
				data.setStep(st);
				data.setError(er);
				data.setStartError(ser);
				data.setRelativeDifferenz(ra);
				data.setGradient(gr);
				data.setNumberOfNewtonSteps(nons);
				data.setCalculationTime(ct);
				
			return s;	
	}
	

bool DKTSDIterationInfo::plot(ostream& log, const ProtocolProperties& pP)const
 {
    bool value = false;
    
    const LongInt  prec = pP.format().precision();
    const LongInt  st   = step();
    const LongReal er   = error();
    const LongReal ser  = startError();
    const LongReal ra   = relativeDifferenz();
				const LongReal gr   = gradient();
    const LongReal nons = numberOfNewtonSteps();
    const LongReal ct   = calculationTime();
    
    IString zeros;
    
    zeros.setZeros(st, 4);
    
    log << setprecision(prec) << uppercase;
    log << setw(4)                                                  << zeros;
    log << scientific;
    log << setw(pP.format().precision()+ pP.format().tabSize() + 6) << ser;   
    log << setw(pP.format().precision()+ pP.format().tabSize() + 6) << er;
    log << setw(pP.format().precision()+ pP.format().tabSize() + 6) << ra;
				log << setw(pP.format().precision()+ pP.format().tabSize() + 6) << gr;
    log << resetiosflags( ::std::ios::scientific );
    log << setw(pP.format().precision()+ pP.format().tabSize() + 6) << nons;
    log << setw(pP.format().precision()+ pP.format().tabSize() + 6) << ct;
    log << endl;
    
   return value;
 }
