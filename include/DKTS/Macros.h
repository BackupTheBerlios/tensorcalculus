/*
 * Copyright (C) Andreas Killaitis
 *               2010 Philipp Waehnert
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
/*!
 ******************************************************************************

 * \author
 *
 */


#ifndef MACROS__H
#define MACROS__H

#include <vector>

#include <math.h>
// #include <string>

#pragma warning(disable:4786)

typedef float       Real;
typedef double      LongReal;
typedef long double BigReal;

typedef LongReal* LongRealPointer;

typedef unsigned int ULongInt;
typedef ULongInt     ULongIntPointer;
typedef int          LongInt;
typedef LongInt*     LongIntPointer;
typedef int          Int;
typedef int          integer;

#define EPS_NULL 1.0e-20
#define NUM_INFINITE 1.0e80

#ifndef Pi
#define Pi 3.141592653589793238462643
#endif


#ifndef MAX
#define MAX(a,b)	(((a)<=(b))? (b):(a))
#endif

#ifndef SIGN
#define SIGN(a)	 (((a)<(0.0))? (-1.0):(1.0))
#endif


#ifndef MIN
#define MIN(a,b)	(((a)<=(b))? (a):(b))
#endif

#ifndef POT2
#define POT2(l)	    (MAX(1,(2<<(l-1))))
#endif


double erf(double x);

// read/write - Attribut
#define ATTRIBUTE(type, getAttr, setAttr)                               \
             private:                                                   \
               type attr_##getAttr;                                     \
             public:                                                    \
               type getAttr() const { return attr_##getAttr; }          \
               This& setAttr( type __x ) { attr_##getAttr = __x; return *this; }

#define PATTRIBUTE(type, getAttr, setAttr)                               \
             public:                                                   \
               type attr_##getAttr;                                     \
             public:                                                    \
               type getAttr() const { return attr_##getAttr; }          \
               This& setAttr( type __x ) { attr_##getAttr = __x; return *this; }


// read only - Attribut
#define CRITICAL_ATTRIBUTE(type, getAttr, setAttr)                      \
             private:                                                   \
               type attr_##getAttr;                                     \
             public:                                                    \
               type getAttr() const { return attr_##getAttr; }          \
             protected:                                                 \
               This& setAttr( type __x ) { attr_##getAttr = __x; return *this; }

// read/write - Klassen-Attribut
#define CLASS_ATTRIBUTE(type, getAttr, setAttr)                         \
             private:                                                   \
               static type attr_##getAttr;                              \
             public:                                                    \
               static type getAttr() const { return attr_##getAttr; }   \
               static This& setAttr(const type& getAttr)                \
                 { attr_##getAttr = getAttr; return *this; }


// read only - Klassen-Attribut
#define CRITICAL_CLASS_ATTRIBUTE(type, getAttr, setAttr)              \
             private:                                                 \
               static type attr_##getAttr;                            \
             public:                                                  \
               static type getAttr() const { return attr_##getAttr; } \
             protected:                                               \
               static type setAttr(const type& getAttr)               \
                 { attr_##getAttr = getAttr; return getAttr(); }


// Flag-Attribut
#define FLAG(fl)                                                                \
           private:                                                             \
             bool flag_##fl;                                                    \
           public:                                                              \
             bool is##fl() const { return flag_##fl; }                          \
             This& set##fl(bool b = true) { flag_##fl = b; return *this; }      \
             This& clear##fl() { flag_##fl = false; return *this; }             \
             This& enable##fl() { set##fl(true); return *this; }                       \
             This& disable##fl() { set##fl(false); return *this; }


// readonly - Flag-Attribut
#define CRITICAL_FLAG(fl)                                                       \
           private:                                                             \
             bool flag_##fl;                                                    \
           public:                                                              \
             bool is##fl() const { return flag_##fl; }                          \
           protected:                                                           \
             This& set##fl(bool b = true) { flag_##fl = b; return *this; }      \
             This& clear##fl() { flag_##fl = false; return *this; }             \
             This& enable##fl() { return set##fl(true); }                       \
             This& disable##fl() { return set##fl(false); }



#define BITFLAG32(name)             \
    private:                          \
        unsigned int bflag_##name;        \
    public:                               \
        This& add##name(const int index) { bflag_##name |= (1 << index); return *this; } \
        This& disableAll##name##s() { bflag_##name = 0; return *this; } \
This& enableAll##name##s() { bflag_##name = 0xFFFFFFFF; return *this; }  \
        bool  is##name(const int index) const { return ((bflag_##name & (1 << index)) > 0); } 
            


// Vorbereitung abgeleiteter Klassen
#define DECLARE(classname)                   \
           private:                          \
             typedef classname This;         \
           public:                           \
             static unsigned long classId()    \
              {                                \
               unsigned long id;               \
	       const char name[] = #classname; \
               int i=0;                        \
               id = 0;                         \
               while (name[i])                    \
                {                              \
                 unsigned long c = ((name[i]-32)<<i);  \
                 id += c;                      \
                 i++;                          \
                }                              \
               return id;                      \
              }


//
//           class ElementName##s : public SizeableArray<ElementType> {}; \
//             ElementName##s cont_##ElementName;			\

// Unterstützung für Container-Beziehungen
#define CONTAINS(ElementType, ElementName)                                      \
           private:                                                              \
   	     typedef std::vector<ElementType> ElementName##s; \
             ElementName##s cont_##ElementName;\
           public:                                                               \
             ElementName##s& get##ElementName##Array() const   \
               { return (ElementName##s&)cont_##ElementName; }   \
             This& add##ElementName(const ElementType& element)                  \
               { cont_##ElementName.push_back(element); \
                 return *this; }          \
             This& removeFirst##ElementName()                                  \
               { cont_##ElementName.erase(cont_##ElementName.begin());                               \
                 return *this; }                                                 \
             This& removeLast##ElementName()                                   \
               { cont_##ElementName.erase(cont_##ElementName.end());                               \
                 return *this; }                                                 \
             This& setSpaceFor##ElementName##sTo(int num)     \
               { cont_##ElementName.reserve(num); \
                 return *this; }           \
             int numberOf##ElementName##s() const                                \
               { return cont_##ElementName.size(); }                 \
             int capacityOf##ElementName##s() const                                \
               { return cont_##ElementName.capacity(); }                 \
             This& removeAll##ElementName##s()                                   \
               { cont_##ElementName.clear(); return *this; }                 \
             This& remove##ElementName##AtPosition(const int index)              \
               { cont_##ElementName.erase(cont_##ElementName.begin()+index);                     \
                 return *this; }                                              \
             ElementType& get##ElementName##At(int index) const            \
               { return (ElementType&)cont_##ElementName.at(index); }           \
             ElementType& first##ElementName() const              \
               { return (ElementType&)cont_##ElementName.front(); }                     \
             ElementType& last##ElementName()                \
               { return (ElementType&)cont_##ElementName.back(); }  \
             bool contains##ElementName(const ElementType& element)       \
              { \
              std::vector<ElementType >::iterator theIterator;                  \
               for (theIterator = cont_##ElementName.begin();               \
                         theIterator != cont_##ElementName.end();           \
                         theIterator++)                                     \
                {                                                           \
                 if(*theIterator == element) return true;                    \
                }                                                           \
               return false;                                                \
              }                                                             \
             This& remove##ElementName(const ElementType& element)      \
              { \
              std::vector<ElementType >::iterator theIterator;                  \
               for (theIterator = cont_##ElementName.begin();               \
                         theIterator != cont_##ElementName.end();           \
                         theIterator++)                                     \
                {                                                           \
                 if(*theIterator == element) cont_##ElementName.erase(theIterator);  \
                }                                                           \
               return *this;                                \
              }                                             \
              This& append##ElementName##sWith(const std::vector<ElementType>& array)      \
              {                                                             \
               cont_##ElementName.reserve(cont_##ElementName.size() + array.size()); \
               std::vector<ElementType >::const_iterator theIterator;       \
               for (theIterator = array.begin();               \
                         theIterator != array.end();           \
                         theIterator++)                                     \
                {                                                           \
                 cont_##ElementName.push_back(*theIterator);  \
                }                                                           \
               return *this;                                \
              }          \
              This& pad##ElementName##ArrayTo(const unsigned int size, const ElementType& el) \
               {                                                    \
                if(cont_##ElementName.size() < size)                \
                  cont_##ElementName.resize(size, el);              \
                return *this;                                       \
               }


class Object
 {
  DECLARE (Object)
 };



#pragma pack (1)

struct __Heritage;

struct __BaseClassEntry
 {
  const char*         className;
  const __Heritage*   heritage;
 };



struct __Heritage
 {
  const char*        thisClass;
  __BaseClassEntry   entry[32];

  bool inheritsFrom(const char* className) const;
 };

#pragma pack ()




#define SUPPORT_INHERITAGE()                                \
  public: static __Heritage __heritage;                     \
  public:                                                   \
    virtual bool inheritsFrom(const char* className) const  \
     {                                                      \
      return __heritage.inheritsFrom(className);            \
     }




#define BEGIN_HERITAGE(className)        \
__Heritage className::__heritage =       \
 {                                       \
  #className,                            \
   {



#define INHERIT_FROM(className)                   \
    { #className, &##className::__heritage },

#define END_HERITAGE()                   \
    { 0, 0 }                             \
   }                                     \
 };


#define isInstanceOf(OBJECT_PTR, CLASS) \
  (&typeid (*OBJECT_PTR) == &typeid (CLASS))


#endif // not defined MACROS__H

