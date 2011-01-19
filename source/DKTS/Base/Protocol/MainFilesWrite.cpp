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

// MainFilesWrite.cpp: Implementierung der Klasse MainFilesWrite.
//
//////////////////////////////////////////////////////////////////////

#include "MainFilesWrite.hpp"


MainFilesWrite::MainFilesWrite()
 {
    	
 }
 
     
MainFilesWrite::~MainFilesWrite()
 {
 
 }


bool MainFilesWrite::mainGnuPlotFileHead(ofstream& log, const IString& terminal, const IString& outputFileName)const
 {
    bool value = false;
				
    log << "#!/gnuplot"											                                     << endl;
    log << "#"												                                             << endl;
    log << "#" 												                                            << endl;    
    log << "#    	G N U P L O T"									                              << endl;
    log << "#    	Version 4.1 patchlevel 0" 							                    << endl;
    log << "#    	last modified Sat Jul  3 00:04:32 CEST 2004" 					   << endl;
    log << "#    	System: MS-Windows 32 bit" 							                   << endl; 
    log << "#" 												                                            << endl;    
    log << "#    	Copyright (C) 1986 - 1993, 1998, 2004" 						        << endl;
    log << "#    	Thomas Williams, Colin Kelley and many others" 				 	<< endl;
    log << "#" 												         										         										      << endl;    
    log << "#    	Type `help` to access the on-line reference manual." << endl;
    log << "#    	The gnuplot FAQ is available from" 						            << endl;
    log << "#    		http://www.gnuplot.info/faq/"																       << endl;
    log << "#"    																				         										         	    << endl;
    log << "#    	Send comments and requests for help to"					        	<< endl;
    log << "#    		<gnuplot-beta@lists.sourceforge.net>"					          << endl;
    log << "#    	Send bugs, suggestions and mods to"						            << endl;
    log << "#    		<gnuplot-beta@lists.sourceforge.net>"				          	<< endl;
    log << "#"    																					         										             << endl;
    log << "set terminal " << terminal     						         				         << endl;
    log << "set output "	  << outputFileName						         										  << endl;
    log << "unset clip points"														         						            << endl;
    log << "set clip one"																								         					        << endl;
    log << "unset clip two"																		         		               << endl;
    log << "set bar 1.000000"																         				             << endl;
    log << "set border 31 front linetype -1 linewidth 1.000"						     << endl;
    log << "set xdata"																		         										         			<< endl;
    log << "set ydata"																		         										         			<< endl;
    log << "set zdata"																		         										         			<< endl;
    log << "set x2data"																	         										         			<< endl;
    log << "set y2data"																		         										         		<< endl;
    log << "set timefmt x \"%d/%m/%y,%H:%M\" "  											         		 << endl;
    log << "set timefmt y \"%d/%m/%y,%H:%M\" "												         		  << endl;
    log << "set timefmt z \"%d/%m/%y,%H:%M\" "											         		  	<< endl;
    log << "set timefmt x2 \"%d/%m/%y,%H:%M\" "										         				 << endl;
    log << "set timefmt y2 \"%d/%m/%y,%H:%M\" "				         					    		<< endl;
    log << "set timefmt cb \"%d/%m/%y,%H:%M\" "											          			<< endl;
    log << "set boxwidth"														         						         		      << endl;
    log << "set style fill  empty border"										         				       << endl;
    log << "set dummy x,y"														         						         		     << endl;
    log << "set format x \"% g\" "									         						         				<< endl;
    log << "set format y \"% g\" "									         						         				<< endl;
    log << "set format x2 \"% g\" "									         						      						<< endl;
    log << "set format y2 \"% g\" "										         						      					<< endl;
    log << "set format z \"% g\" "											         					        				<< endl;
    log << "set format cb \"% g\" "										         			         					<< endl;
    log << "set angles radians"													         						         			<< endl;
    log << "set grid nopolar"												         						           				<< endl;
    log << "set grid xtics nomxtics ytics nomytics noztics nomztics \\"					                          << endl;
    log << " nox2tics nomx2tics noy2tics nomy2tics nocbtics nomcbtics"					                           << endl;
    log << "set grid layerdefault   linetype 0 linewidth 1.000,  linetype 0 linewidth 1.000" 		       << endl;
    log << "set key title \"\" " 													         						         						         						         		<< endl;
    log << "set key inside right top vertical Right noreverse noinvert enhanced samplen 4 spacing 1";
    log << "width 0 height 0 autotitles box linetype -2 linewidth 1.000"                             	<< endl;
    log << "unset label"													         						         						   	<< endl;
    log << "unset arrow"														         						         						   << endl;
    log << "unset style line"									         						         							  << endl;
    log << "unset style arrow"									         						         							 << endl;
    log << "set style histogram clustered gap 2 title  offset 0, 0, 0"	<< endl;
    log << "unset logscale"													         						         			    << endl;
    log << "set logscale y 10"											         						         					 << endl;
    log << "set offsets 0, 0, 0, 0"								         						         				<< endl;
    log << "set pointsize 1"														         						         		   << endl;
    log << "set encoding default"								         						         						<< endl;
    log << "unset polar"													         						         		       	<< endl;
    log << "unset parametric"												         						         			  	<< endl;
    log << "unset decimalsign"											         						         					 << endl;
    log << "set view 60, 30, 1, 1"									         						         				<< endl;
    log << "set samples 100, 100"											         						         			<< endl;
    log << "set isosamples 10, 10"											         						         		<< endl;
    log << "set surface"														         						         		       << endl;
    log << "unset contour"												         						         				     << endl;
    log << "set clabel '%8.3g'" 													         						         		<< endl;
    log << "set mapping cartesian" 										         						        			<< endl;
    log << "set datafile separator whitespace"										         		  		<< endl;
    log << "unset hidden3d"													         						         		    	<< endl;
    log << "set cntrparam order 4"												         						         	<< endl;
    log << "set cntrparam linear"													         						         	<< endl;
    log << "set cntrparam levels auto 5"								         						        << endl;
    log << "set cntrparam points 5"												         			            << endl;
    log << "set size ratio 0 1,1"													         						         	<< endl;
    log << "set origin 0,0"													         						         			    << endl;
    log << "set style data linespoints"										         					        << endl;
    log << "set style function lines"											         				          << endl;
    log << "set xzeroaxis linetype -2 linewidth 1.000"							          << endl;
    log << "set yzeroaxis linetype -2 linewidth 1.000"							          << endl;
    log << "set zzeroaxis linetype -2 linewidth 1.000"							          << endl;
    log << "set x2zeroaxis linetype -2 linewidth 1.000"							         << endl;
    log << "set y2zeroaxis linetype -2 linewidth 1.000"							         << endl;
    log << "set ticslevel 0.5"														         						         	 	<< endl;
    log << "set mxtics default"													         						         			<< endl;
    log << "set mytics default"													         						         			<< endl;
    log << "set mztics default"													         						         			<< endl;
    log << "set mx2tics default"												         						         			<< endl;
    log << "set my2tics default"													         						         		<< endl;
    log << "set mcbtics default"												         						         			<< endl;
    log << "set xtics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 2"		        << endl;
    log << "set ytics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 autofreq"		 << endl;
    log << "set ztics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autofreq"	<< endl;
    log << "set nox2tics"														         						         						         						         		   << endl;
    log << "set noy2tics"														         						         						         						         		   << endl;
    log << "set cbtics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 autofreq" 	<< endl;
    log << "set title \"\"  offset character 0, 0, 0 font \"\" norotate"    				                  << endl; //setzt einen Plottitel
    log << "set timestamp bottom "								        				         						         						         		   << endl; //fuegt datum und zeit in plot ein
    log << "set timestamp \"\" offset character 0, 0, 0 font \"\" norotate"								         		    << endl;
    log << "set rrange [ * : * ] noreverse nowriteback  # (currently [0.000000:10.0000] )"						  << endl;
    log << "set trange [ * : * ] noreverse nowriteback  # (currently [-5.00000:5.00000] )"						  << endl;
    log << "set urange [ * : * ] noreverse nowriteback  # (currently [-5.00000:5.00000] )"						  << endl;
    log << "set vrange [ * : * ] noreverse nowriteback  # (currently [-5.00000:5.00000] )"					   << endl;
    log << "set xlabel \"\"  offset character 0, 0, 0 font \"\" textcolor lt -1 norotate"						   << endl;
    log << "set x2label \"\" offset character 0, 0, 0 font \"\" textcolor lt -1 norotate"						   << endl;
    log << "set xrange [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )"						  << endl;
    log << "set x2range [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )"						 << endl;
    log << "set ylabel \"\"  offset character 0, 0, 0 font \"\" textcolor lt -1 rotate by 90"					<< endl;
    log << "set y2label \"\"  offset character 0, 0, 0 font \"\" textcolor lt -1 rotate by 90"				<< endl;
    log << "set yrange [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )"						  << endl;
    log << "set y2range [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )"						 << endl;
    log << "set zlabel \"\"  offset character 0, 0, 0 font \"\" textcolor lt -1 norotate"						   << endl;
    log << "set zrange [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )"						  << endl;
    log << "set cblabel \"\"  offset character 0, 0, 0 font \"\" textcolor lt -1 norotate"						  << endl;
    log << "set cbrange [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )"						 << endl;
    log << "set zero 1e-008"														         						         						         						         		<< endl;
    log << "set lmargin -1"													         						         						         						         			 << endl;
    log << "set bmargin -1"												         						         						         						         				 << endl;
    log << "set rmargin -1"												         						         						         						         			 	<< endl;
    log << "set tmargin -1"														         						         						         						         	 	<< endl;
    log << "set locale \"C\""												         						         						         						         			<< endl;
    log << "set pm3d explicit at s"										         						         						         						        << endl;
    log << "set pm3d scansautomatic"										         						         						         						       << endl;
    log << "set pm3d interpolate 1,1 flush begin noftriangles nohidden3d corners2color mean"		    << endl;
    log << "set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB" 			            << endl;
    log << "set palette rgbformulae 7, 5, 15"												         						         						         		<< endl;
    log << "set colorbox default"													         						         						         		           << endl;
    log << "set colorbox vertical origin 0.9,0.2 size 0.1,0.63 bdefault"								         						   << endl;
    log << "set loadpath" 													         						         						         						         					<< endl;
    log << "set fontpath"													         						         						         						         			   << endl;
    log << "set fit noerrorvariables"												         						         						         						    << endl;
    log << "plot "; 											
 
   return value;
 }  


bool MainFilesWrite::mainGnuPlotFileEnd(ofstream& log)const
 {
    bool value = false;
    
    log << endl;	
    log << "#" << '\t' << "EOF" << endl;
    
   return value;
 }   
    

bool MainFilesWrite::mainTeXFileHead(ofstream& log, const ProtocolProperties& pP)const
 {
    bool value = true;
    
    const IString institutNamePart1 = pP.institutNamePart1();
    const IString institutNamePart2 = pP.institutNamePart2();
    const IString topicString       = pP.topicString();
    const IString authorString      = pP.authorString();
    const IString cityString        = pP.cityString();
				const IString projektName       = pP.projektName();
    
    log << "\\documentclass [12pt]{report}"                                              			<< endl;
//    log << "\\usepackage {german}"                                                       			<< endl;
    log << "\\usepackage {a4}"                                                           			<< endl;
    log << "\\usepackage {amsmath}"                                                      			<< endl;
    log << "\\usepackage {amssymb}"                                                      			<< endl;
    log << "\\usepackage {amsthm}"                                                       			<< endl;
    log << "\\usepackage {alltt}"                                                        			<< endl;
    log << "\\usepackage {rotating}"                                                     			<< endl;
    log << "\\usepackage {subfigure}"                                                    			<< endl;
    log << "\\usepackage {epsfig}"                                                       			<< endl;
    log << "\\usepackage [latin1]{inputenc}"                                             			<< endl;
    log << "\\usepackage [dvips]{color}"                                                 			<< endl;
    log << "\\usepackage{fancyhdr}"                                                      			<< endl;
    log << "\\usepackage{float}"														         						         						         		     << endl;
    log << endl;
    log << endl;
    log << "\\markboth{chapter}{section}"                                                			<< endl;
    log << endl;
    log << "\\def\\refS#1{\\ref{#1}, Seite \\pageref{#1}}"                   					          << endl;
    log << endl;
    log << "\\graphicspath{{images/}}"                                       					          << endl;
    log << "\\DeclareGraphicsRule{.eps.gz}{eps}{.eps.gz}{`cat #1 | gzip -d}" 					          << endl;
    log << endl;
    log << "\\ifx\\undefined\\pdfpageheight"                                             			<< endl;
    log << '\t' << "\\def\\image#1#2{\\epsfig{file=#1,width=#2}}"            					          << endl;
    log << "\\else"                                                                      			<< endl;
    log << '\t' << "\\def\\image#1#2{\\epsfig{file=#1.pdf,width=#2}}"        					          << endl;
    log << "\\fi"                                                                        			<< endl;
    log << endl;
    log << "\\def\\kapitel#1"                                                            			<< endl;
    log << '\t' << "{"                                                                   			<< endl;
    log << '\t' << '\t' << "\\refstepcounter{chapter}"                                      << endl;
    log << '\t' << '\t' << "\\chapter*{#1}\\addcontentsline{toc}{chapter}{{\\thechapter}\\hspace{2.1mm} #1}" 	<< endl;
    log << '\t' << "}"                                                                   			<< endl;
    log << endl;
    log << endl;
    log << "%\\setlength{\\parskip}{1.5ex plus0.5ex minus0.5ex}"                            << endl;
    log << "\\setlength{\\parindent}{0.0em}"                                            		 	<< endl;
    log << endl; 
    log << "\\setlength{\\unitlength}{1.0cm} \\setcounter {secnumdepth}{5}"                 << endl;
    log << "\\setlength{\\unitlength}{1.0cm} \\setcounter {tocdepth}{5}"                    << endl;
    log << endl;
    log << "% Fuzz -------------------------------------------------------------------"     << endl;
    log << "\\hfuzz2pt % Don't bother to report over-full boxes if over-edge is < 2pt"      << endl;
    log << endl;
    log << endl;
    log << endl;
    log << "%\\input{MATHS/HEADMATH.TEX}"    									    << endl;
    log << "%\\input{MATHS/ENVIMATHDE.TEX}"  									    << endl;
    log << endl;
    log << "\\begin{document}"               									    << endl;
    log << endl;
    log << "\\pagenumbering{Roman}"  << endl;
    log << "\\title{"                << endl;
    log << "	\\begin{Large}"         << endl;
    log << "		{\\sc "                << institutNamePart1 << "}\\\\[5ex]" << endl;
    log << "		"                      << institutNamePart2 << "\\\\[6ex]"  << endl;
    log << "	\\end{Large}"           << endl;
    log << "	\\begin{large}"         << endl;
    log << "		\\textbf{"             << topicString       << "}\\\\[5ex]" << endl;
				log << "  "                      << projektName                       << endl;
    log << "	\\end{large}"           << endl;
    log << "	\\\\[5ex]"              << endl;
    log << "	\\begin{small}"         << endl;
    log << "	\\author{"              << authorString      << "}"          << endl;
    log << "	\\end{small}"           << endl;
    log << "}"                       << endl;
    log << endl;
    log << "\\date{"                 << cityString << "\\\\ \\today}"     << endl;
    log << endl;
    log << endl;
    log << "\\maketitle"             << endl;
    log << endl;
    log << "\\tableofcontents"       << endl;
    log << endl;
    log << endl;
    log << "\\pagenumbering{arabic}" << endl;
    log << endl;
    
   return value;
 }  

    
bool MainFilesWrite::mainTeXFileEnd(ofstream& log)const
 {
    bool value = false;
    
    log << endl;
    log << "\\end{document}" << endl;
    
   return value;
 }
