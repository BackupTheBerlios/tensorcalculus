/*
 * Copyright (C) 2010, 2011 Stefan Handschuh
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

#include <iostream>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream> 

#include "Utilities/Index.hpp"

	std::string vector2string(const std::vector<int> vector) {
		std::stringstream result;
		
		result << "(";
		for (int n = 0, i = vector.size(); n < i; n++) {
			result << vector[n];
			if (n < i-1) {
				result << ", ";	
			}
		}	
		result << ")";
		return result.str();
		
	}

int main() {
	std::vector<int> summations(4);
	
	summations[0] = 4;
	summations[1] = 5;
	summations[2] = 3;
	summations[3] = 4;
	
	for (Index index(summations); !index.end(); ++index) {
		//std::cout << index << std::endl;
	}
	
	std::vector<int> fixedIndices(2);
	fixedIndices[0] = 0;
	fixedIndices[1] = 2;
	
	std::vector<int> fixedValues(2);
	fixedValues[0] = 0;
	fixedValues[1] = 0;
	

	for (Index index(summations, fixedIndices, fixedValues); !index.end(); ++index) {
		//std::cout << index << std::endl;
	}
	
	std::vector<int> singleSummation(1);
	
	singleSummation[0] = 34;
	
	Index index(singleSummation);
	
	for (index.begin(); !index.end(); ++index) {
		//std::cout << index << std::endl;
	}

	std::cout << std::endl << " FixedEmpty " << std::endl << "===========" << std::endl;
	std::vector<int> all(4);
	all[0] = 0;
	all[1] = 1;
	all[2] = 2;
	all[3] = 3;
	std::vector<int> fixedAll(4);
	fixedAll[0] = 1;
	fixedAll[1] = 1;
	fixedAll[2] = 2;
	fixedAll[3] = 0;

	for (Index index(summations, all, fixedAll); !index.end(); ++index) {
		std::cout << index << std::endl;// this should only be written once
	}
	
	return 0;	
}
