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

// Timer.cpp: Implementierung der Klasse Timer.
//
//////////////////////////////////////////////////////////////////////

#include "Timer.hpp"


Timer& Timer::startTiming()
 {
    startTimeS  = clock();
   return *this;
 }


LongReal Timer::elapsedTimeSec()
 {
    LongReal endTimeS = clock();

   return ((endTimeS - startTimeS)/CLOCKS_PER_SEC);
 }


IString Timer::date() const
 {
    time_t t;
    time(&t);
       
   return IString(ctime(&t));
 }


IString Timer::timeStamp() const
 {
    time_t zeitstempel;
    
    tm *nun;
    
    zeitstempel = time(0);
    
    nun = localtime(&zeitstempel);
    
    IString timeStamp  = IString(setZeros(nun->tm_mday))+IString(".")+IString(setZeros(nun->tm_mon+1))+IString(".")+IString(nun->tm_year+1900)+IString("_");
            timeStamp += IString(setZeros(nun->tm_hour))+IString("-")+IString(setZeros(nun->tm_min))  +IString("-")+IString(setZeros(nun->tm_sec)); 

   return timeStamp;
 }	    
	    

IString Timer::setZeros(const LongInt& number)const
 {
    IString value = "";
    
    if(number < 10)
     {
        value += IString("0");
        value += IString(number);
     }
    else
     {
        value += IString(number);
     }	
    
   return value;
 }  
     
