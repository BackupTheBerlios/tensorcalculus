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

// ProtocolFormat.cpp: Implementierung der Klasse ProtocolFormat.
//
//////////////////////////////////////////////////////////////////////

#include "ProtocolFormat.hpp"

ProtocolFormat::ProtocolFormat()
 {
    (*this)
				 .setDefaultParameter()
    ;
 }
	

ProtocolFormat::ProtocolFormat(const ProtocolFormat& pF)
 {
    (*this) = pF;
 }


ProtocolFormat::ProtocolFormat(const LongInt& prec, const LongInt& tabSize, const IString& position, const LongInt& tableSize)
 {
    (*this)
     .setPrecision     (prec)
     .setTabSize       (tabSize)
     .setTablePosition (position)
     .setMaxTableSize  (tableSize)
    ;
 }
 

ProtocolFormat::~ProtocolFormat()
 {
 
 }
 
 
ProtocolFormat ProtocolFormat::setDefaultParameter()
 {
    (*this)
     .setPrecision     (4)
     .setTabSize       (3)
     .setTablePosition ("ht")
     .setMaxTableSize  (35)
    ;
   
   return (*this); 
 }             
 
 
ProtocolFormat& ProtocolFormat::operator = (const ProtocolFormat& pF)
 {
     (*this)
      .setPrecision     (pF.precision())
      .setTabSize       (pF.tabSize())
      .setTablePosition (pF.tablePosition())
      .setMaxTableSize  (pF.maxTableSize())
     ;
     
    return (*this);
 }
     

ostream& operator << (ostream& s, const ProtocolFormat& pF)
 {
    s << "Precision     = " << pF.precision()     << endl;
    s << "tabSize       = " << pF.tabSize()       << endl;
    s << "tablePosition = " << pF.tablePosition() << endl;
    s << "maxTableSize  = " << pF.maxTableSize()  << endl;
    
   return s;
 }   
