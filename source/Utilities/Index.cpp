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

#include "Utilities/Index.hpp"
#include <ostream>
#include <stdexcept> 
#include <stdlib.h>
#include <iostream>
#include <numeric>
#include <functional>


int indexOf(const std::vector<int> &vector, int element) {
	for (int n = 0, i = vector.size(); n < i; n++) {
		if (element == vector[n]) {
			return n;
		}	
	}
	return -1;
}

void Index::init(const std::vector<int> * bounds) {
	(*this).bounds = bounds;
	boundsCount = bounds == 0 ? 0 : (*bounds).size();
	current.resize(boundsCount);
	position = 0;
	changeCount = -1;
	increment = 0;
	int factor = 1;
	
	partialProducts.resize(boundsCount);
	if (bounds != 0) {
		for (int n = 0; n < boundsCount; n++) {
			partialProducts[n] = (*bounds)[n]*factor;
			factor = partialProducts[n];
		}
	}

}

Index::Index(const std::vector<int> * bounds) {
	init(bounds);
}

Index::Index(const std::vector<int> &bounds) {
	init(&bounds);
	fixedIndices = 0;
}

Index::Index(const std::vector<int> &bounds, const std::vector<int> &fixedIndices, const std::vector<int> &fixedValues) {
	init(&bounds);

	elementCount = size();

	int i = fixedIndices.size();

#ifdef ARGUMENT_CHECKS_ON
	if (i >= boundsCount) {
		throw std::invalid_argument("The number of fixed indices >= number of bounds.");
	}
#endif
	(*this).fixedIndices = &fixedIndices;		
		
#ifdef ARGUMENT_CHECKS_ON
	if (i != fixedValues.size()) {
		throw std::invalid_argument("The number of fixed indices differs from the corresponding value.");
	}
#endif

	int k;
	
	for (int n = 0; n < i; n++, k = fixedIndices[n]) {
		k = fixedIndices[n];
		elementCount /= bounds[k];
#ifdef ARGUMENT_CHECKS_ON
		if (k >= i) {
			throw std::invalid_argument("Fixing index exceeds index count.");
		}
#endif
		current[k] = fixedValues[n];
		if (k > 0) {
			position += fixedValues[n]*partialProducts[k-1];
		} else {
			position += fixedValues[n];	
		} 
	}
}

int Index::getIncrement() const {
	return increment;
}

const std::vector<int>& Index::getCurrent() const {
	return current;	
}
		
Index& Index::moveTo(long n) {
	for (int i = 0; i < boundsCount; i++) {
//		if (fixedIndices == 0 || indexOf((*fixedIndices), i) == -1) { 
			current[n] = n % (*bounds)[i];
			n /= (*bounds)[i];	
//		}
	}
	position = n;
	return *this;
}
		
int Index::getChangeCount() const {
	return changeCount;	
}
		
const std::vector<long> Index::getPartialProducts() const {
	return partialProducts;	
}

const long Index::size() const {
	return partialProducts[boundsCount-1];
}

const long Index::count() const {
	return elementCount;
}

Index& Index::begin() {
	for (int n = 0; n < boundsCount; n++) {
		if (fixedIndices == 0 || indexOf((*fixedIndices), n) == -1) {
			current[n] = 0;
		}	
	}	
	changeCount = -1;
	position = 0;
	return *this;
}
		
bool Index::end() const {
	return changeCount == boundsCount;	
}

const std::vector<int> * Index::getBounds() const {
	return bounds;
}

		
Index& Index::operator ++ () {
	if (fixedIndices != 0) {
		for(changeCount = 0; changeCount < boundsCount && (current[changeCount] == (*bounds)[changeCount]-1 || indexOf((*fixedIndices), changeCount) > -1); changeCount++) {
			if (indexOf((*fixedIndices), changeCount) == -1) { 
				current[changeCount] = 0;	
			} else {
				if (changeCount > 0) {
					position += partialProducts[changeCount-1]*((*bounds)[changeCount]-1);
				} else {
					position += (*bounds)[0]-1;
				}
			}
		}
	} else {
		for(changeCount = 0; changeCount < boundsCount && current[changeCount] == (*bounds)[changeCount]-1; changeCount++) {
			current[changeCount] = 0;
		}
	}
	if (boundsCount > changeCount) {
		++current[changeCount];	
		++position;
	}
	increment++;
	return *this;
}

const int Index::operator [] (int n) const {
	return current[n];
}
		
Index::operator long () const {
	return position;
}

Index::operator const std::vector<int>& () {
	return current;	
}
		
long Index::getPosition() const {
	return position;	
}

std::ostream& operator << (std::ostream &stream, const Index &index) {
	const std::vector<int> &current = index.getCurrent();
	
	int size = current.size();
	
	stream << (int) index << "/" << index.getPartialProducts()[size-1] << " = (";
	for (int n = 0; n < size; n++) {
		stream << current[n];
		if (n < size-1) {
			stream << ", ";	
		}
	}	
	stream << ") change = " << index.getChangeCount();
	return stream;	
}

bool operator == (const Index &index1, const Index &index2) {
	return index1.getCurrent() == index2.getCurrent();	
}

bool operator == (const Index &index1, const std::vector<int> &index2) {
	return index1.getCurrent() == index2;
}

bool operator == (const std::vector<int> &index1, const Index &index2) {
	return index2 == index1;
}

bool operator != (const Index &index1, const Index &index2) {
	return !(index1 == index2);
}

bool operator != (const Index &index1, const std::vector<int> &index2) {
	return !(index1 == index2);
}

bool operator != (const std::vector<int> &index1, const Index &index2) {
	return index2 != index1;
}


