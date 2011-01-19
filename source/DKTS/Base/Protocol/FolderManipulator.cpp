/*
 * Copyright (C) Marcel Schuster
 *               2011 Philipp Wähnert
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

// FolderManipulator.cpp: Implementierung der Klasse FolderManipulator.
//
//////////////////////////////////////////////////////////////////////

#include "FolderManipulator.hpp"

#if defined(_MSC_VER) && (defined(_WIN32) || defined(_WIN64))
#include <direct.h>
#define TENSORCALCULUS_WINDOWS
#elif defined(__GNUC__) && defined(__linux__)
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#define TENSORCALCULUS_LINUX
#else
#error "Unsupported platform"
#endif

#include "IString.hpp"
#include "SimpleException.hpp"

bool FolderManipulator::createFolder(const IString& name)const
 {
#ifdef TENSORCALCULUS_LINUX
    if(mkdir(name, -1) == -1) // mode = -1 yields mode && ~umask == ~umask
#elif defined(TENSORCALCULUS_WINDOWS)
    if(_mkdir(name) == -1)
#endif
     {
       return false;
     }
   return true;
 }
 

bool FolderManipulator::changeFolder(const IString& name)const
 {
#ifdef TENSORCALCULUS_LINUX
    if(chdir(name) == -1) 
#elif defined(TENSORCALCULUS_WINDOWS)
    if(_chdir(name) == -1)
#endif
     {
        throw SimpleException(IString("Konnte nicht in das Verzeichnis wechseln" ));
      
       return false;
     }
   return true;
 }    

bool FolderManipulator::folderUp()const
 {
   return changeFolder("..");
 }
    

IString FolderManipulator::currentFolder()const
 {
    char puffer[200];
#ifdef TENSORCALCULUS_LINUX
    if(getcwd(puffer,sizeof(puffer)) == NULL)
#elif defined(TENSORCALCULUS_WINDOWS)
    if(_getcwd(puffer,sizeof(puffer)) == NULL)
#endif
     {
        throw SimpleException(IString("Fehler bei getcwd ..."));

       return IString();
     }
   
  return puffer;
 }
 
