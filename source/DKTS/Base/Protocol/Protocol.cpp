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

// Protocol.cpp: Implementierung der Klasse Protocol.
//
//////////////////////////////////////////////////////////////////////

#include "Protocol.hpp"
#include "FolderManipulator.hpp"


Protocol::Protocol()
 {
    (*this)
     .setDefaultParameter()
    ;
 }
 
Protocol::Protocol(const ProtocolProperties& pP)
:attr_properties(pP)
 {
 
 }


Protocol::Protocol(const Protocol& p)
 {
    (*this) = p;
 }   
     
     
Protocol::~Protocol()
 {
 
 }


Protocol Protocol::setDefaultParameter()
 {
    (*this)
    ;
    
   return (*this);
 }    
 

Protocol& Protocol::operator = (const Protocol& p)
 {
    (*this)
     .setProperties(p.properties())
    ;
    
    const LongInt size = p.numberOfEntrys();
    
    for(LongInt i=0; i<size; i++)
     {
        const ProtocolDataPointer element = p.getEntryAt(i);
        
        include(element);
     }   
    
   return (*this);
 }


ostream& operator << (ostream& s, const Protocol& p)
 {
	   p.plot2(s);
	
	  return s;
	}			


Protocol& Protocol::include(const ProtocolDataPointer element)
 { 
    if(isInstanceOf(element, SimpleIterationData) == true)
     {
        const SimpleIterationData& data = (SimpleIterationData&)*element;
            
        add(data);
     }
    else if(isInstanceOf(element, DKTSDData) == true)
     {
        const DKTSDData& data = (DKTSDData&)*element;
        
        add(data);
     }
    else if(isInstanceOf(element, DKTSDIterationInfo) == true)
     {
        const DKTSDIterationInfo& data = (DKTSDIterationInfo&)*element;
	
	       add(data);
     }	 
    else if(isInstanceOf(element, DescriptionData) == true) 
     {
        const DescriptionData& data = (DescriptionData&)*element;
        
        add(data);
     }
    else if(isInstanceOf(element, SimpleIterationDataBlock) == true)
     {
        const SimpleIterationDataBlock& data = (SimpleIterationDataBlock&)*element;
	
	       add(data);
     }
    else if(isInstanceOf(element, DKTSDDataBlock) == true)
     {
        const DKTSDDataBlock& data = (DKTSDDataBlock&)*element;
	
	       add(data);
     }
    else if(isInstanceOf(element, DKTSDIterationInfoBlock) == true)
     {
        const DKTSDIterationInfoBlock& data = (DKTSDIterationInfoBlock&)*element;
	
	       add(data);
     }	 	
    else    
     {
        throw SimpleException(IString("Warning In Protocol::include (const ProtocolDataPointer element), Bad ProtocolDataType !"));
     }
         
   return(*this);     
 }  


bool Protocol::plot2(ostream& log) const
 {
    bool value = false;
    
    const ProtocolProperties& pP = properties();
        
    const LongInt size = numberOfEntrys ();
    
    for(LongInt i=0; i<size; i++)
     {
        const ProtocolDataPointer element = getEntryAt(i);
	
        element->plot(log, pP);
     }
   
   return value; 
 }      


bool Protocol::plot2PlotFile(const IString& mainFoldername)const
 {
    bool value = false;
    
    const ProtocolProperties& pP = properties();
    
    const LongInt size = numberOfEntrys();
          LongInt number = 1;
    
    //Erzeugen der main-dateien fuer gnu-Plot
    const IString folderName = "plot";
    const IString homeFolder = createFolderTree(mainFoldername, pP, folderName);
   
    //Erstellen der einzelnen Datendatein
    for(LongInt i=0; i<size; i++)
     {
	       const ProtocolDataPointer element = getEntryAt(i);
	
	       const IString blockIndex = setZeros(number);
	      
	       const bool fit = element->generatedPlotFiles(blockIndex, pP);
	       
	       if(fit == true)
	        {
		          number++;
		       }   
     }	      
    
    //Wechseln ins urspruengliche Verzeichnis
    FolderManipulator fm;
    fm.changeFolder(homeFolder);
    
   return value; 
 }


bool Protocol::plot2TeXFile(const IString& mainFoldername, const IString& projektName)const
 {
    bool value = false;
   
    const ProtocolProperties& pP = properties();
				
				const IString TeXFileName = pP.TeXFileName();
//				pP.setProjektName(projektName);
    
    const LongInt size  = numberOfEntrys();
    
    //Erstellen der Verzeichnisstruktur
    FolderManipulator fm;
				
    fm.createFolder(mainFoldername);            
    fm.changeFolder(mainFoldername);
    
				ofstream lastOut(TeXFileName);    
    
				fm.folderUp();
    
    const IString folderName = "tex";
    const IString homeFolder = createFolderTree(mainFoldername, pP, folderName);
    
    ofstream out("main.tex");    
    
    //Erstellen des Kopfes der Projektdatei fuer LaTeX
    MainFilesWrite mfw;
    
    mfw.mainTeXFileHead(lastOut, pP);
    mfw.mainTeXFileHead(out, pP);    
    
    const IString data = "data";
    fm.createFolder(data);
    fm.changeFolder(data);
    
    const IString path = pP.time() + IString("/") + folderName + IString("/") + data;
    
    //Erstellen der einzelnen Datendatein
    for(LongInt i=0; i<size; i++)
     {
        const ProtocolDataPointer element = getEntryAt(i);
	
	       const IString blockIndex = setZeros(i+1);
	
	       const IString name = element->generatedTeXFiles(blockIndex, pP);
	
	       if(name != "0")
	        {	 
	           lastOut << "\\input{" << path << "/" << name << "}" << endl;          
	           out     << "\\input{" << data << "/" << name << "}" << endl;
	        }   
	
     }	      
    
    //Erstellen des Endes der Projektdatei fuer LaTeX
    mfw.mainTeXFileEnd(lastOut);    
    mfw.mainTeXFileEnd(out);
    
    //Wechseln ins urspruengliche Verzeichnis
    fm.changeFolder(homeFolder);
    
   return value; 
 }


bool Protocol::plotAll(const IString& mainFoldername, const IString& name)const
 {
    bool value = false;
    
    plot2TeXFile(mainFoldername, name);
    plot2PlotFile(mainFoldername);
    plotAll2File(mainFoldername, name);
     
   return value; 
 }


bool Protocol::plotAll2File(const IString& mainFoldername, const IString& name)const
 {
    bool value = false;
   
    const ProtocolProperties& pP = properties();
   
    //Erstellen der Verzeichnisstruktur
    const IString folderName("doc");
    const IString homeFolder = createFolderTree(mainFoldername, pP, folderName);  
    
    ofstream write2(name+".log");
    plot2(write2);

    FolderManipulator fm;
    fm.changeFolder(homeFolder);
    
   return value;
 }   
    
    
Protocol& Protocol::add(const SimpleIterationData& data)
 {
    ProtocolDataPointer pd = new SimpleIterationData(data);
    
    addEntry(pd);
    
   return (*this);
 }  


Protocol& Protocol::add(const DKTSDData& data)
 {
    ProtocolDataPointer pd = new DKTSDData(data);
    
    addEntry(pd);
    
   return (*this);
 }   


Protocol& Protocol::add(const DKTSDIterationInfo& data)
 {
    ProtocolDataPointer pd = new DKTSDIterationInfo(data);
    
    addEntry(pd);
   
   return (*this);
 }
    

Protocol& Protocol::add(const DescriptionData& data)
 {
    ProtocolDataPointer pd = new DescriptionData(data);
    
    addEntry(pd);
   
   return (*this);
 }
 
 
Protocol& Protocol::add(const SimpleIterationDataBlock& data)
 {
    ProtocolDataPointer pd = new SimpleIterationDataBlock(data);
    
    addEntry(pd);
    
   return (*this);
 }   


Protocol& Protocol::add(const DKTSDDataBlock& data)
 {
    ProtocolDataPointer pd = new DKTSDDataBlock(data);
    
    addEntry(pd);
   
   return (*this);
 }   


Protocol& Protocol::add(const DKTSDIterationInfoBlock& data)
 {
    ProtocolDataPointer pd = new DKTSDIterationInfoBlock(data);
    
    addEntry(pd);
   
   return (*this);
 } 


IString Protocol::createFolderTree(const IString& mainFoldername, const ProtocolProperties& pP, const IString& name) const
 {
   //Erstellen der Verzeichnisstruktur
    const IString zeitstempel(pP.time());
    
    FolderManipulator fm;
    
    const IString homeFolder = fm.currentFolder();
    
    fm.createFolder(mainFoldername);            
    fm.changeFolder(mainFoldername);
    fm.createFolder(zeitstempel);
    fm.changeFolder(zeitstempel);
    fm.createFolder(name);
    fm.changeFolder(name);
   
   return homeFolder;
 }  


IString Protocol::setZeros(const LongInt& number)const
 {
    IString value = "";
    
    if(number < 10)
     {
        value += IString("00");
        value += IString(number);
     }
    else if(number < 100)
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
