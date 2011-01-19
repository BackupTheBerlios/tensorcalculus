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

#include "IString.hpp"
#include "SimpleException.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


IString::IString()
 {
  _buffer     = 0;
  _bufferSize = 0;
  _length     = 0;

  *this = "";
 }



IString::IString(const char* string)
 {
  _buffer     = 0;
  _bufferSize = 0;
  _length     = 0;

  *this = string;
 }



IString::IString(const unsigned char* string)
 {
  _buffer     = 0;
  _bufferSize = 0;
  _length     = 0;
  *this = (char*)string;
 }


IString::IString(const IString& string)
 {
  _buffer     = 0;
  _bufferSize = 0;
  _length     = 0;
  if (&string != 0) *this = string;
 }






IString::IString(signed char c)
 {
  _buffer     = 0;
  _bufferSize = 0;
  _length     = 0;
  char t[20];
  sprintf(t, "%c", c);
  *this = t;
 }




IString::IString(unsigned char c)
 {
  _buffer     = 0;
  _bufferSize = 0;
  _length     = 0;
  char t[20];
  sprintf(t, "%c", c);
  *this = t;
 }





IString::IString(float f)
 {
  _buffer     = 0;
  _bufferSize = 0;
  _length     = 0;
  char t[20];
  sprintf(t, "%f", double(f));
  *this = t;
  stripTrailing("0");
  stripTrailing(".");
 }






IString::IString(const LongReal& d, const LongInt& prec)
 {
  _buffer     = 0;
  _bufferSize = 0;
  _length     = 0;
  
  char t[20];  
  IString s;
  
  s = IString("%.") + IString(prec) + IString("E");    
        
  sprintf(t, s, d);
  
  *this = t;
  //stripTrailing("0");
  stripTrailing(".");
 }




IString::IString(long double d)
 {
  _buffer     = 0;
  _bufferSize = 0;
  _length     = 0;
  char t[20];
  sprintf(t, "%Lf", d);
  *this = t;
  stripTrailing("0");
  stripTrailing(".");
 }





IString::IString(signed short s)
 {
  _buffer     = 0;
  _bufferSize = 0;
  _length     = 0;
  char t[20];
  sprintf(t, "%ld", (signed long)s);
  *this = t;
 }







IString::IString(unsigned short s)
 {
  _buffer     = 0;
  _bufferSize = 0;
  _length     = 0;
  char t[20];
  sprintf(t, "%lu", (unsigned long)s);
  *this = t;
 }







IString::IString(signed long l)
 {
  _buffer     = 0;
  _bufferSize = 0;
  _length     = 0;
  char t[20];
  sprintf(t, "%ld", l);
  *this = t;
 }






IString::IString(unsigned long l)
 {
  _buffer     = 0;
  _bufferSize = 0;
  _length     = 0;
  char t[20];
  sprintf(t, "%lu", l);
  *this = t;
 }








IString::IString(signed int i)
 {
  _buffer     = 0;
  _bufferSize = 0;
  _length     = 0;
  char t[20];
  sprintf(t, "%d", i);
  *this = t;
 }






IString::IString(unsigned int i)
 {
  _buffer     = 0;
  _bufferSize = 0;
  _length     = 0;
  char t[20];
  sprintf(t, "%u", i);
  *this = t;
 }



IString::~IString()
 {
  flush();
 }


IString& IString::operator = (const char* string)
 {
  if (string == 0) string = "";

  int  newLen = (int) strlen(string);
  if (newLen+1 > _bufferSize)
   {
    flush();
    _bufferSize = newLen+1;
    _buffer = new char[_bufferSize];
   }

  strcpy(_buffer, string);
  _length = newLen;

  return *this;
 }


IString& IString::operator = (const IString& string)
 {
  *this = (char*)string;
  return *this;
 }



double IString::asDouble() const
 {
  double d;

  if (_buffer != 0) d = atof(_buffer);
  else d = 0;
  return d;
 }




signed char IString::asChar() const
 {
  if (_buffer != 0) return _buffer[0];
  else return '\0';
 }




signed short IString::asShort() const
 {
  return (signed short)asLong();
 }




signed int IString::asInt() const
 {
   LongReal temp = asDouble();
  return (signed int)temp;
 }




signed long IString::asLong() const
 {
  signed long l;

  if (_buffer != 0) sscanf(_buffer, "%ld", &l);
  else l = 0;
  return l;
 }




unsigned char IString::asUChar() const
 {
  if (_buffer != 0) return _buffer[0];
  else return '\0';
 }




unsigned short IString::asUShort() const
 {
  return (unsigned short)asLong();
 }




unsigned int IString::asUInt() const
 {
  return (unsigned int)asLong();
 }




unsigned long IString::asULong() const
 {
  unsigned long l;

  if (_buffer != 0) sscanf(_buffer, "%lu", &l);
  else l = 0;
  return l;
 }





// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



IString::operator char* () const
 {
  if (_buffer == 0) throw(SimpleException( IString("[IString::operator char*]: access to NULL buffer")));

  return _buffer;
 }



// ************************************************************
// * �berladene Stream-Operatoren
// ************************************************************

ostream& operator << (ostream& os, const IString& s)
 {
  os << (char*)s;
  return os;
 }




// ************************************************************
// * gesch�tzte Klassendienste 
// ************************************************************

void IString::flush()
 {
  if (_buffer) delete _buffer;
  _buffer     = 0;
  _bufferSize = 0;
  _length     = 0;
 }



int IString::length() const
 {
  return _length;
 }






IString& IString::operator += (const IString& string)
 {  
  int newLen = length() + string.length();
  if (newLen+1 > _bufferSize)
   { // der Puffer reicht nicht f�r beide Teilstrings
    char* p = new char[newLen+1];

    strcpy(p, _buffer);
    strcpy(p+_length, string._buffer);

    flush();
    _bufferSize = newLen+1;
    _buffer = p;
    _length = newLen;
   }
  else
   { // der Puffer reicht f�r beide Teilstrings
    strcpy(_buffer+_length, string._buffer);
    _length = newLen;
   }

  return *this;
 }





bool IString::operator == (const char* string) const
 {
  if ((int)length() != (int)strlen(string)) return false;
  else return (strcmp(_buffer, string) == 0);
 }





bool IString::operator != (const char* string) const
 {
  return !(*this == string);
 }





bool IString::operator < (const char* string) const
 {
  return (strcmp(_buffer, string) < 0);
 }





bool IString::operator > (const char* string) const
 {
  return (strcmp(_buffer, string) > 0);
 }




bool IString::operator == (const IString& string) const
 {
  if (length() != string.length()) return false;
  else return (strcmp(_buffer, string._buffer) == 0);
 }





bool IString::operator != (const IString& string) const
 {
  return !(*this == string);
 }





bool IString::operator < (const IString& string) const
 {
  return (strcmp(_buffer, string._buffer) < 0);
 }





bool IString::operator > (const IString& string) const
 {
  return (strcmp(_buffer, string._buffer) > 0);
 }







IString operator + (const IString& s1, const IString& s2)
 {
  IString s(s1);
  s += s2;
  return s;
 }


IString operator + (const IString& s1, const char* s2)
 {
  return IString(s1) + IString(s2);
 }





IString IString::subString(int from, int to) const
 {
  IString s;
  if (from >= 1 && to <= length())
   {
    char* t = new char[to-from+2];
    memcpy(t, _buffer+from-1, to-from+2);
    t[to-from+1] = 0;
    s = t;
    delete[] t;
   }
  else 
   {
    IString exc = IString("[IString::subString]: illegal substring range [");
    exc += IString(from) + IString(";") + IString(to) + IString("]");
    SimpleException e(exc);
    throw(e);
   }

  return s;
 }



IString& IString::stripLeading(const IString& t)
 {
  int l = t.length();

  while (strncmp(t._buffer, _buffer, l) == 0)
   {
    memcpy(_buffer, _buffer+l, _length);
    _length -= l;
   }

  return *this;
 }



IString& IString::stripTrailing(const IString& t)
 {
  int l = t.length();

  while (strncmp(t._buffer, _buffer+_length-l, l) == 0)
   {
    _length -= l;
    _buffer[_length] = 0;
   }

  return *this;
 }



IString& IString::strip(const IString& t)
 {
  stripLeading(t);
  stripTrailing(t);

  return *this;
 }



char& IString::operator [] (int index) const
 {
  if (index < 1 || index > length()) 
   {
    IString exctext = "[IString::operator []]: illegal index: ";
    exctext += IString(index);
    SimpleException e(exctext);
    throw(e);
   }

  return _buffer[index-1];
 }





IString convertNumberFormat(const IString& t, int sourceBase, int targetBase)
 {
  // ************************************************************************
  // * zun�chst den Wert des Source-Strings im sourceBase-System ermitteln
  // ************************************************************************
  char* p = t;

   // Ende finden
  int index = 0;
  while (p[index]) index++;  
  index--;

  // Wert ermitteln
  int number = 1;
  int allValue = 0;
  while (index >= 0)
   {
    if (p[index] == '-') allValue = -allValue;
    else
     {
      char digit = p[index];
      int value = 0;
      if (sourceBase == 256) // Zeichen-Index
        value = digit;
      else
       {
        switch (digit)
         {
          case '0': case '1': case '2': case '3': case '4':
          case '5': case '6': case '7': case '8': case '9':
            value = digit - '0';
            break;

          case 'a': case 'A':  value = 10; break;
          case 'b': case 'B':  value = 11; break;
          case 'c': case 'C':  value = 12; break;
          case 'd': case 'D':  value = 13; break;
          case 'e': case 'E':  value = 14; break;
          case 'f': case 'F':  value = 15; break;
          default: allValue = 0; index = -1; break; // Fehler im Zahlenformat!
         }
       }

      if (value < sourceBase) // g�ltige Ziffer?
       {
        allValue += number*value;
        number *= sourceBase;
       }
      else
       {
        index = -1;
        allValue = 0;
       }
     }
    index--;
   }


  // ************************************************************************
  // * jetzt den ermittelten Wert ins Ziel-Zahlensystem umwandeln
  // ************************************************************************

  char T[20];
  index = sizeof (T)-1;  // rueckw�rts schreiben

  T[index--] = '\0';

  if (allValue < 0)
    { T[index--] = '-'; allValue = -allValue; }  // Vorzeichen beachten

  const char digit[] = "0123456789abcdef";

  if (allValue == 0) return IString("0");

  while (allValue != 0)
   {
    int value = allValue % targetBase;
    T[index] = digit[value];
    index--;
    allValue /= targetBase;
   }
  index++;


  return IString(&T[index]);
 }


IString convertString(const IString& t, const LongInt& targetBase)
 {
    const LongInt len = t.length();

    IString s = "";

    for(LongInt i=1; i<=len; i++)
     {
        IString c = IString(char(t[i]));
       
        s += convertNumberFormat(c, 256, targetBase);
     }

   return s;
 }


IString& IString::b2o()
 {
    (*this) = convertNumberFormat((*this), 2, 8);
   return (*this);
 }


IString& IString::b2d()
 {
    (*this) = convertNumberFormat((*this), 2, 10);
   return (*this);
 }


IString& IString::b2x()
 {
    (*this) = convertNumberFormat((*this), 2, 16);
   return (*this);
 }


IString& IString::o2b()
 {
    (*this) = convertNumberFormat((*this), 8, 2);
   return (*this);
 }


IString& IString::o2d()
 {
    (*this) = convertNumberFormat((*this), 8, 10);
   return (*this);
 }


IString& IString::o2x()
 {
    (*this) = convertNumberFormat((*this), 8, 16);
   return (*this);
 }


IString& IString::d2b()
 {
    (*this) = convertNumberFormat((*this), 10, 2);
   return (*this);
 }


IString& IString::d2o()
 {
    (*this) = convertNumberFormat((*this), 10, 8);
   return (*this);
 }


IString& IString::d2x()
 {
    (*this) = convertNumberFormat((*this), 10, 16);
   return (*this);
 }


IString& IString::x2b()
 {
    (*this) = convertNumberFormat((*this), 16, 2);
   return (*this);
 }


IString& IString::x2o()
 {
    (*this) = convertNumberFormat((*this), 16, 8);
   return (*this);
 }


IString& IString::x2d()
 {
    (*this) = convertNumberFormat((*this), 16, 10);
   return (*this);
 }


IString& IString::c2b()
 {
    (*this) = convertString((*this), 2);
   return (*this);
 }


IString& IString::c2o()
 {
    (*this) = convertString((*this), 8);
   return (*this);
 }


IString& IString::c2d()
 {
    (*this) = convertString((*this), 10);
   return (*this);
 }


IString& IString::c2x()
 {
    (*this) = convertString((*this), 16);
   return (*this);
 }


IString IString::b2o(const IString& binString)
 {
   return IString(binString).b2o();
 }


IString IString::b2d(const IString& binString)
 {
   return IString(binString).b2d();
 }


IString IString::b2x(const IString& binString)
 {
   return IString(binString).b2x();
 }


IString IString::o2b(const IString& octString)
 {
   return IString(octString).o2b();
 }


IString IString::o2d(const IString& octString)
 {
   return IString(octString).o2d();
 }


IString IString::o2x(const IString& octString)
 {
   return IString(octString).o2x();
 }


IString IString::d2b(const IString& decString)
 {
   return IString(decString).d2b();
 }


IString IString::d2o(const IString& decString)
 {
   return IString(decString).d2o();
 }


IString IString::d2x(const IString& decString)
 {
   return IString(decString).d2x();
 }


IString IString::x2b(const IString& hexString)
 {
   return IString(hexString).x2b();
 }


IString IString::x2o(const IString& hexString)
 {
   return IString(hexString).x2o();
 }


IString IString::x2d(const IString& hexString)
 {
   return IString(hexString).x2d();
 }


IString IString::c2b(const IString& string)
 {
   return IString(string).c2b();
 }


IString IString::c2o(const IString& string)
 {
   return IString(string).c2o();
 }


IString IString::c2d(const IString& string)
 {
   return IString(string).c2d();
 }


IString IString::c2x(const IString& string)
 {
   return IString(string).c2x();
 }


bool IString::isBinaryDigits() const
 {
    char* p = _buffer;

    bool b = true;
  
    if(p != 0)
     {
        while(*p)
         {
            if(!(*p >= '0' && *p <= '1'))
             {
                b = false;
                break;
             }
            p++;
         }
     }
     
   return b;
 }


bool IString::isDigits() const
 {
    char* p = _buffer;

    bool b = true;
  
    if(p != 0)
     {
        while(*p)
         {
            if(!(*p >= '0' && *p <= '9'))
             {
                b = false;
                break;
             }
            p++;
         }
     }
     
   return b;
 }


bool IString::isAlphabetic() const
 {
    char* p = _buffer;

    bool b = true;
  
    if(p != 0)
     {
        while(*p)
         {
            if(!((*p >= 'a' && *p <= 'z') || (*p >= 'A' && *p <= 'Z')))
             {
                b = false;
                break;
             }
            p++;
         }
     }
     
   return b;
 }


bool IString::isAlphanumeric() const
 {
    char* p = _buffer;

    bool b = true;
  
    if(p != 0)
     {
        while(*p)
         {
            if (!((*p >= 'a' && *p <= 'z') || (*p >= 'A' && *p <= 'Z') || (*p >= '0' && *p <= '9')))
             {
                b = false;
                break;
             }
            p++;
         }
     }
     
   return b;
 }


bool IString::isHexDigits() const
 {
    char* p = _buffer;

    bool b = true;
  
    if(p != 0)
     {
        while(*p)
         {
            if(!((*p >= 'a' && *p <= 'f') || (*p >= 'A' && *p <= 'F') || (*p >= '0' && *p <= '9')))
             {
                b = false;
                break;
             }
            p++;
         }
     }
     
   return b;
 }


bool IString::isControl() const
 {
    char* p = _buffer;

    bool b = true;
  
    if(p != 0)
     {
        while(*p)
         {
            if(!(*p >= 0x00 && *p <= 0x1f) && *p != 0x7f)
             {
                b = false;
                break;
             }
            p++;
         }
     }
     
   return b;
 }


bool IString::isLowerCase() const
 {
    char* p = _buffer;

    bool b = true;
  
    if(p != 0)
     {
        while(*p)
         {
            char c1 = *p;
            char c2 = c1-35;//tolower(c1);
            
	    if(c1 != c2)
             {
                b = false;
                break;
             }
            p++;
         }
     }
     
   return b;
 }


bool IString::isUpperCase() const
 {
    char* p = _buffer;

    bool b = true;
  
    if(p != 0)
     {
        while(*p)
         {
            char c1 = *p;
            char c2 = c1+35; //toupper(c1);
            
	    if(c1 != c2)
             {
                b = false;
                break;
             }
            p++;
         }
     }
     
   return b;
 }


IString& IString::leftJustify(const LongInt& width, const IString& padString)
 {
    IString t = (*this);

    while(t.length() < width) 
     {
        t += padString;
     }
     
    (*this) = t.subString(1, width);

   return (*this);
 }


IString& IString::rightJustify(const LongInt& width, const IString& padString)
 {
    IString t = padString;

    while(t.length() < width) 
     {
        t += padString;
     }
     
    t += (*this);

    (*this) = t.subString(t.length()-width+1, t.length());

   return (*this);
 }


IString& IString::center(const LongInt& width, const IString& padString)
 {
    IString t = (*this);

    while(t.length() < width)
     {
        t = padString + t + padString;
     }
  
    LongInt middle = t.length() / 2;
    LongInt len    = width/2;

    (*this) = t.subString(middle-len+1, middle+len);

   return (*this);
 }


IString& IString::change(const IString& oldWord, const IString& newWord)
 {
    LongInt maxLen = length() * newWord.length() / oldWord.length() + 10;
    LongInt index = 0;
    
    char* buffer = new char[maxLen];
    char* src    = _buffer;
    char* dst    = buffer;

    while(*src)
     {
        // jetzt kommt das Wort
	if(strncmp(src, oldWord, oldWord.length()) == 0) 
         {
            strcpy(dst, newWord);
            
	    src += oldWord.length();
            dst += newWord.length();
         }
        else
         {
            *dst = *src;
      
            src++;
            dst++;
         }
     }
    *dst = 0;

    (*this) = IString(buffer);

    delete[] buffer;

   return (*this);
 }


IString& IString::copy(const LongInt& count)
 {
    IString t = "";

    for(LongInt i=0; i<count; i++) 
     {
        t += (*this);
     }	

    (*this) = t;

   return (*this);
 }


LongInt IString::occurrencesOf(const IString& string) const
 {
    LongInt index = 0;
    LongInt count = 0;
  
    while(index < length())
     {
        if(strncmp(&_buffer[index], string, string.length()) == 0) 
         {
	    count++;
	 }   
        index++;
     }
   
   return count;
 }


LongInt IString::indexOf(const IString& string) const
 {
    LongInt index = 0;
  
    while(index < length())
     {
        if(strncmp(&_buffer[index], string, string.length()) == 0) 
         {
	    return index+1;
	 }
        else 
         {
	    index++;
	 }   
     }
   
   // nicht gefunden
   return 0; 
 }


LongInt IString::indexOfAnyOf(const IString& string) const
 {
    LongInt index = 0;
    
    while(index < length())
     {
       for(LongInt i=1; i<=string.length(); i++)
        {
           if(_buffer[index] == string[i]) 
	    {
	       return index+1;
	    }   
        }
       index++;
     }
   
   // nicht gefunden
   return 0; 
 }


LongInt IString::indexOfAnyBut(const IString& string) const
 {
    LongInt index = 0;
  
    while(index < length())
     {
        bool found = false;
        
	for(LongInt i=1; i<=string.length(); i++)
         {
            if(_buffer[index] == string[i]) 
	     {
	        found = true;
	     }	
         }
        if(!found) 
	 {
	    return index+1;
	 }   
        else 
	 {
	    index++;
	 }   
     }
   
   // nicht gefunden
   return 0; 
 }


LongInt IString::indexOfWord(const LongInt& wordNumber, char delimeter) const
 {
    LongInt index = 0;
    LongInt  word = 0;
    
    while(index < length())
     {
        // Vorruecken zum n�chsten Wort
        if(index == 0 || _buffer[index] == delimeter)
         {
            while(_buffer[index] == delimeter) 
	     {
	        index++;
	     }	
            word++;
         }

        // ist es dieses Wort?
        if(word == wordNumber) 
	 {
	   return index+1;
	 }  

        // aufessen des Wortes
        while (_buffer[index] != delimeter && index < length()) 
	 {
	    index++;
	 }   
      }

   return 0;
 }


LongInt IString::numWords(char delimeter) const
 {
    LongInt index = 0;
    LongInt words = 0;
  
    while(index < length())
     {
        // Vorruecken zum n�chsten Wort
        if(index == 0 || _buffer[index] == delimeter)
         {
            while(_buffer[index] == delimeter) 
	     {
	        index++;
	     }	
            words++;
         }

        // aufessen des Wortes
        while(_buffer[index] != delimeter && index < length()) 
	 {
	    index++;
	 }   
     }

   return words;
 }


LongInt IString::lengthOfWord(const LongInt& wordNumber, char delimeter) const
 {
    LongInt index = 0;
    LongInt  word = 0;
  
    while(index < length())
     {
        // Vorruecken zum n�chsten Wort
        if(index == 0 || _buffer[index] == delimeter)
         {
            while(_buffer[index] == delimeter) 
	     {
	        index++;
	     }	 
            word++;
          }

        // ist es dieses Wort?
        if(word == wordNumber)
         {
            LongInt count = 0;
             
            while(_buffer[index] != delimeter && index < length()) 
             {
                index++;
                count++;
             }
           return count;
         }

        // aufessen des Wortes
        while (_buffer[index] != delimeter && index < length())
         {
	    index++;
	 }   
     }

   return 0;
 }


IString& IString::upperCase()
 {
    for(LongInt i=0; i<length(); i++)
     {
        _buffer[i] = _buffer[i]+35;
     }
     	
   return *this;
 }


IString& IString::lowerCase()
 {
    for(LongInt i=0; i<length(); i++)
     {
	_buffer[i] = _buffer[i]-35;
     }
     	
   return *this;
 }


IString& IString::setZeros(const LongInt& a, const LongInt& fieldsize)
 {
    LongInt numberOfZeros = 1;
    
    (*this) = IString (""); 
    
    if(a<10)
     {
        numberOfZeros = fieldsize - 1;
     }
    else if(a<100)
     {
        numberOfZeros = fieldsize -2;
     }
    else if(a<1000)
     {
        numberOfZeros = fieldsize - 3;
     }
    else
     {
        numberOfZeros = fieldsize -4;
     }	  

    for(LongInt i=0; i < numberOfZeros; i++)
     {
        (*this) += IString ("0");
     }
    
    (*this) += IString (a);
   
   return (*this);
 }    	


IString IString::word(const LongInt& wordNumber, const char delimeter) const
 {
    LongInt index     = 0;
    LongInt wordIndex = 0;
  
    while (index < length())
     {
        // Vorruecken zum n�chsten Wort
        if(index == 0 || _buffer[index] == delimeter)
         {
            while (_buffer[index] == delimeter)
	     {
	        index++;
             }   
	    wordIndex++;
         }

        // ist es dieses Wort?
        if (wordIndex == wordNumber)
         {
            // L�nge ermitteln
            LongInt count = 0;
      
            while (_buffer[index+count] != delimeter && index+count < length()) 
             {
                count++;
             }
            if (_buffer[index+count] == delimeter)
	     {
	        count--;
	     }

            char* buf = new char[count+1];
      
            memset(buf, 0x00, 255);
            strncpy(buf, &_buffer[index], count+1);
      
            IString t(buf);
      
            delete[] buf;
      
           return t;
          }

        // aufessen des Wortes
        while (_buffer[index] != delimeter && index < length()) 
         {
            index++;
         }
     }

    return IString("");
 }


IString& IString::reverse()
 {
    IString t = *this;
    LongInt l = length();
    
    for(LongInt i=1; i<=l; i++) 
     {
        (*this)[i] = t[l-i+1];
     }
   
   return *this;
 }
     
