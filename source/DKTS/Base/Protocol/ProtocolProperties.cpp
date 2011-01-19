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

// ProtocolProperties.cpp: Implementierung der Klasse ProtocolProperties.
//
//////////////////////////////////////////////////////////////////////

#include "ProtocolProperties.hpp"


ProtocolProperties::ProtocolProperties(const LongInt& prec, const LongInt& tabSize, const IString& position, const LongInt& tableSize)
 :attr_format(prec, tabSize, position, tableSize)
 {
    (*this)
     .setDefaultParameter()
    ; 
 }
 
 
ProtocolProperties::ProtocolProperties(const ProtocolFormat& pF)
 {
    (*this)
     .setFormat(pF)
     .setDefaultParameter()
    ; 
 }
 

ProtocolProperties::ProtocolProperties(const ProtocolProperties& pP)
 {
    (*this) = pP;
 }   
  
     
ProtocolProperties::~ProtocolProperties()
 {
 
 }


ProtocolProperties ProtocolProperties::setDefaultParameter()
 {
    Timer timer;
    
    const IString time = timer.timeStamp(); 
    
    (*this)
     .setTime(time)
     .setInstitutNamePart1 ("Max-Planck-Institute for Mathematics")
     .setInstitutNamePart2 ("in the Sciences")
     .setTopicString       ("Kronecker Sum Truncation")
     .setAuthorString      ("Mike Espig, Marcel Schuster")
     .setCityString        ("Leipzig")
					.setTeXFileName       ("lastRun.tex")
					.setProjektName       ("test")
    ;
    
   return (*this);
 }
 
 
ProtocolProperties& ProtocolProperties::operator = (const ProtocolProperties& pP)
 {
    (*this)
     .setFormat            (pP.format())
     .setTime              (pP.time())
     .setInstitutNamePart1 (pP.institutNamePart1())
     .setInstitutNamePart2 (pP.institutNamePart2())
     .setTopicString       (pP.topicString())
     .setAuthorString      (pP.authorString())
     .setCityString        (pP.cityString())
					.setTeXFileName       (pP.TeXFileName())
					.setProjektName       (pP.projektName())
    ;
    
   return (*this);
 }   
    
     
   
ostream& operator << (ostream& s, const ProtocolProperties& pP)
 {
    s <<                     pP.format()            << endl;
    s << "date: "         << pP.time()              << endl;
    s << "name1       = " << pP.institutNamePart1() << endl;
    s << "name2       = " << pP.institutNamePart2() << endl;
    s << "topic       = " << pP.topicString()       << endl;
    s << "author      = " << pP.authorString()      << endl;
    s << "city        = " << pP.cityString()        << endl;
				s << "TeXFileName = " << pP.TeXFileName()       << endl;
				s << "ProjektName = " << pP.projektName()       << endl;
    
   return s;
 }   
    
