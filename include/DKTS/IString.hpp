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

#ifndef __IString__
#define __IString__


#include <iostream>
#include "Macros.h"

using namespace std;

class SimpleException;

class IString
 {
  // ****************************************************************************
  // * gesch�tzte Implementierungsvariablen
  // ****************************************************************************

  private:
    char*   _buffer;      // Zeiger auf den Zeichenpuffer, der den String enth�lt
    int     _bufferSize;  // Gr��e des Puffers
    int     _length;      // aktuelle L�nge des Strings (_length < _bufferSize)


  // ****************************************************************************
  // * gesch�tzte Methoden zur Implementierung
  // ****************************************************************************

  protected:
    void flush(); // l�scht den Puffer und setzt die Implementierungsvaiablen zur�ck


  // ****************************************************************************
  // * Konstruktion / Destruktion
  // ****************************************************************************

  public:
    // Standardkonstruktor, initialisiert den String auf den Leerstring
    IString();

    // die folgenden Konstruktoren erwarten als Argument Variablen verschiedener
    // Typen, deren Wert vom String repr�sentiert wird, etwa IString(3.1) = "3.1"
    IString(const char* string);
    IString(const unsigned char* string);
    IString(const IString& string);
    IString(signed char c);
    IString(unsigned char c);
    IString(float f);
    IString(const LongReal& d, const LongInt& prec = 4);
    IString(long double d);
    IString(signed short s);
    IString(unsigned short s);
    IString(signed long l);
    IString(unsigned long l);
    IString(signed int i);
    IString(unsigned int i);

    // liest einen String aus der userResourceLibrary
    // IString(const IResourceId& resourceId);

    // gibt den verwendeten Speicher frei und r�umt die Instanz auf
    ~IString();


  // ****************************************************************************
  // * �berladene Operatoren
  // ****************************************************************************

  public:
    // Typ-Konversion: IString -> (char*)
    operator char* () const;

    // liefert das Zeichen an der mit index bezeichneten Stelle.
    // Der Index ist 1-basiert, die Elemente heissen [1][2][3]..[length()]
    char& operator [] (int index) const;

    // kopiert die angegebene Zeichenkette in this
    IString& operator = (const char* string);

    // kopiert den angegebenen String in this
    IString& operator = (const IString& string);

    // kopiert den angegebenen String in this
    IString& operator += (const IString& string);

    // Die folgenden Operatoren vergleichen diesen String mit einem angegebenen.
    // Es existieren Versionen f�r Zeichenzeiger und Srings
    bool operator == (const char* string) const;
    bool operator != (const char* string) const;
    bool operator < (const char* string) const;
    bool operator > (const char* string) const;

    bool operator == (const IString& string) const;
    bool operator != (const IString& string) const;
    bool operator < (const IString& string) const;
    bool operator > (const IString& string) const;


  // ****************************************************************************
  // * Konvertierungen
  // ****************************************************************************

  public:
    // interpretiert den String als Ausdruck eines der folgenden Typen und liefert
    // dessen Wert zur�ck

    double asDouble() const;

    signed char  asChar() const;
    signed short asShort() const;
    signed int   asInt() const;
    signed long  asLong() const;

    unsigned char  asUChar() const;
    unsigned short asUShort() const;
    unsigned int   asUInt() const;
    unsigned long  asULong() const;

    // die folgenden Operatoren fassen den String als numerischen Wert in
    // einem bestimmten Zahlensystem auf, wandeln den Wert in ein neues
    // Zahlensystem um und setzen den String auf den neuen Wert, also etwa
    //   IString("255").d2x() = "ff"         ,
    //   IString("1111111").b2d() = "127"    usw.

    IString& b2o(); // bin�r -> oktal
    IString& b2d(); // bin�r -> dezimal
    IString& b2x(); // bin�r -> hexadezimal

    IString& o2b(); // oktal -> bin�r
    IString& o2d(); // oktal -> dezimal
    IString& o2x(); // oktal -> hexadezimal

    IString& d2b(); // dezimal -> bin�r
    IString& d2o(); // dezimal -> oktal
    IString& d2x(); // dezimal -> hexadezimal

    IString& x2b(); // hexadezimal -> bin�r
    IString& x2o(); // hexadezimal -> oktal
    IString& x2d(); // hexadezimal -> dezimal

    IString& c2b(); // character -> bin�r
    IString& c2o(); // character -> oktal
    IString& c2d(); // character -> dezimal
    IString& c2x(); // character -> hexadezimal


    // die folgenden Operatoren haben dieselbe Aufgabe wie ihre �quivalente oben,
    // sind jedoch statisch und liefern nur eine konvertierte Kopie des Strings

    static IString b2o(const IString& binString); // bin�r -> oktal
    static IString b2d(const IString& binString); // bin�r -> dezimal
    static IString b2x(const IString& binString); // bin�r -> hexadezimal

    static IString o2b(const IString& octString); // oktal -> bin�r
    static IString o2d(const IString& octString); // oktal -> dezimal
    static IString o2x(const IString& octString); // oktal -> hexadezimal

    static IString d2b(const IString& decString); // dezimal -> bin�r
    static IString d2o(const IString& decString); // dezimal -> oktal
    static IString d2x(const IString& decString); // dezimal -> hexadezimal

    static IString x2b(const IString& hexString); // hexadezimal -> bin�r
    static IString x2o(const IString& hexString); // hexadezimal -> oktal
    static IString x2d(const IString& hexString); // hexadezimal -> dezimal

    static IString c2b(const IString& string);    // character -> bin�r
    static IString c2o(const IString& string);    // character -> oktal
    static IString c2d(const IString& string);    // character -> dezimal
    static IString c2x(const IString& string);    // character -> hexadezimal


  // ****************************************************************************
  // * String-Tests und Informationen
  // ****************************************************************************

  public:
    bool isBinaryDigits() const;
    bool isDigits() const;
    bool isAlphabetic() const;
    bool isAlphanumeric() const;
    bool isHexDigits() const;
    bool isControl() const; // Steuercodes 0x00-0x1f, 0x7f
    bool isLowerCase() const;
    bool isUpperCase() const;

    // liefert die L�nge des Strings in Zeichen
    int length() const;


  // ****************************************************************************
  // * Wort-Methoden
  // ****************************************************************************

  public:
    // liefert das Wort an der angegebenen Stelle
    IString word(const LongInt& wordNumber, const char delimeter = ' ') const;

    // liefert den Index des ersten Zeichen des angegebenen Wortes
    int indexOfWord(const LongInt& wordNumber, char delimeter = ' ') const;

    // liefert die Anzahl der W�oerter in diesem String
    LongInt numWords(char delimeter = ' ') const;

    // liefert die L�nge der Wortes mit dem angegebenen Index
    LongInt lengthOfWord(const LongInt& wordNumber, char delimeter = ' ') const;


  // ****************************************************************************
  // * String-Formatierungen
  // ****************************************************************************

  public:
    // liefert den Teilbereich [from][from+1]...[to-1][to] des Strings
    // liegen die Bereiche from oder to au�erhalb der Stringgrenzen, 
    // so gelten die Stringgrenzen als Stringgrenzen des Substrings
    // Die Indizes from und to sind 1-basiert
    IString subString(int from, int to) const;

    // f�gt den String in einen String der L�nge width linksb�ndig ein und 
    // f�llt rechts, falls Platz �brig ist, mit dem F�llstring padString auf, 
    // etwa
    //   IString("abc").leftJustify(10) = "abc       "
    // Ist der vorhandene String l�nger als width Zeichen, wird der Rest
    // abgeschnitten
    IString& leftJustify(const LongInt& width, const IString& padString = IString(" "));

    // f�gt den String in einen String der L�nge width rechtsb�ndig ein und 
    // f�llt links, falls Platz �brig ist, mit dem F�llstring padString auf, 
    // etwa
    //   IString("abc").leftJustify(10) = "       abc"
    // Ist der vorhandene String l�nger als width Zeichen, wird der Rest
    // abgeschnitten
    IString& rightJustify(const LongInt& width, const IString& padString = IString(" "));

    // f�gt den String in einen String der L�nge width zentriert ein und 
    // f�llt links und rechts, falls Platz �brig ist, mit dem F�llstring 
    // padString auf, etwa
    //   IString("abc").center(10) = "    abc    "
    // Ist der vorhandene String l�nger als width Zeichen, wird der Rest
    // abgeschnitten
    IString& center(const LongInt& width, const IString& padString = IString(" "));

    // entfernt den angegebenen String am Anfang des Strings, auch mehrfach
    // z.B.: IString("abcabcabcXXX").stripLeading("abc") = "XXX";
    IString& stripLeading(const IString& stripString = IString(" "));

    // entfernt den angegebenen String am Ende des Strings, auch mehrfach
    // z.B.: IString("XXXabcabcabc").stripTrailing("abc") = "XXX";
    IString& stripTrailing(const IString& stripString = IString(" "));

    // entfernt den angegebenen String am Anfang und Ende des Strings
    // z.B.: IString("   XXX     ").strip(" ") = "XXX";
    IString& strip(const IString& t = IString(" "));

    // ersetzt alle Vorkommen des angegebenen Suchwortes oldWord durch 
    // das Wort newWord
    IString& change(const IString& oldWord, const IString& newWord);

    // die Methode copy erzeugt eine Mehrfachkopie von sich selbst
    //   IString("abc").copy(3) = "abcabcabc";
    IString& copy(const LongInt& count = 1);

    // liefert die Anzahl der Vorkommen des angegebenen Strings
    LongInt occurrencesOf(const IString& searchString) const;

    // liefert den Index des ersten Vorkommens des angegebenen Strings
    LongInt indexOf(const IString& string) const;

    // liefert den Index des ersten Vorkommens irgendeines Zeichens
    // des angegebenen Strings
    LongInt indexOfAnyOf(const IString& string) const;

    // liefert den Index des ersten Vorkommens irgendeines Zeichens,
    // das NICHT im angegebenen String enthalten ist
    LongInt indexOfAnyBut(const IString& string) const;

    // wandelt alle Kleinbuchstaben in Gro�buchstaben um
    IString& upperCase();

    // wandelt alle Gro�buchstaben in Kleinbuchstaben um
    IString& lowerCase();

    //Einfuegen von fuerenden nullen
    IString& setZeros(const LongInt& a, const LongInt& fieldsize);

    // ****************************************************************
    /*
    IString &overlayWith( const IString &aString,
                      unsigned       index        = 1,
                      char           padCharacter = ' '),
 IString& remove            ( unsigned startPos );
 IString& remove            ( unsigned startPos,
                      unsigned numChars );

 */

 IString& reverse           ( );
 };



// ****************************************************************************
// * �berladene Operatoren
// ****************************************************************************

IString  operator +  (const IString& s1, const IString& s2);
ostream& operator << (ostream& os, const IString& s);


#endif // not defined __IString__
