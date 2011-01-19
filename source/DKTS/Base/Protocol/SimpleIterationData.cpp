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

// SimpleIterationData.cpp: Implementierung der Klasse SimpleIterationData.
//
//////////////////////////////////////////////////////////////////////

#include "SimpleIterationData.hpp"
#include <typeinfo>


SimpleIterationData::SimpleIterationData (const LongInt& st, const LongReal& er)
 {
    (*this)
     .setStep             (st)
     .setError            (er)     
     .setDefaultParameter ()
    ; 
 }
 

SimpleIterationData::SimpleIterationData(const SimpleIterationData& sID)
 {
    (*this) = sID;     
 } 
 
     
SimpleIterationData::~SimpleIterationData()
 {
 
 }


SimpleIterationData SimpleIterationData::setDefaultParameter()
 {
    (*this)
    ;
   
   return (*this); 
 }    


SimpleIterationData& SimpleIterationData::operator = (const SimpleIterationData& sD)
 {
    (*this)
     .setStep  (sD.step())
     .setError (sD.error())
    ; 
    
   return (*this);
 }  
   
   
bool SimpleIterationData::operator == (const ProtocolData& D)
 { 
    bool value = false;
    
    if(isInstanceOf(&D, SimpleIterationData) == true)
     {
        SimpleIterationData& data = (SimpleIterationData&) D;

        if(step()==data.step() && fabs(error()-data.error())<EPS_NULL)
         {
            value = true;
         }   
     }

   return value;   
 }


ostream& operator << (ostream& s, const SimpleIterationData& sD)
 {
	   const LongInt  st  = sD.step();
    const LongReal er  = sD.error();
				
				s << st << '\t';
				s << scientific << uppercase;  // << Zehnerpotenz << Grossbuchstaben 
				s << er << endl;
				s << resetiosflags( ::std::ios::scientific );
			
			return s;
	}
				

istream& operator >> (istream& s, SimpleIterationData& data)
 {
				LongInt  st = 0;
				LongReal er = 0;
				
				s >> st;
				s >> er;
				
				data.setStep(st);
				data.setError(er);
				
			return s;
	}
				

bool SimpleIterationData::plot(ostream& log, const ProtocolProperties& pP) const
 {
    bool value = true;
    
    const LongInt  st = step();    
    const LongReal er = error();    
    
    const LongInt prec     = pP.format().precision();
    const LongInt freesize = pP.format().tabSize();  
    
    IString zeros;
    
    zeros.setZeros(st, 4);
    
    log << setw(4) << zeros << setw(prec + freesize + 6) << er << endl;

   return value;    	
    
 }
