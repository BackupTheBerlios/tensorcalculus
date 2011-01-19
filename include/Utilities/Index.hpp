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
#ifndef INDEX_HPP_
#define INDEX_HPP_

#include <vector>
#include <iosfwd>

class Index {
	private:
		const std::vector<int> * bounds;
		
		const std::vector<int> * fixedIndices;
		
		std::vector<int> current;
		
		int changeCount;
		
		int boundsCount;
		
		int position;
		
		long elementCount;

		std::vector<long> partialProducts;
		
		void init(const std::vector<int> * bounds);
		
		int increment;

	public:
		Index(const std::vector<int> * bounds = 0);

		Index(const std::vector<int> &bounds);
		
		Index(const std::vector<int> &bounds, const std::vector<int> &fixedIndices, const std::vector<int> &fixedValues);

		const std::vector<int>& getCurrent() const;
					
		Index& moveTo(long n);
		
		int getChangeCount() const;
		
		const std::vector<long> getPartialProducts() const;
		
		const long size() const;
		
		const long count() const;

		Index& begin();
		
		bool end() const;
		
		const std::vector<int> * getBounds() const;
		
		Index& operator ++ ();
		
		int getIncrement() const;

		operator long () const;
		
		operator const std::vector<int>& ();
		
		const int operator [] (int) const;
		
		long getPosition() const;
		
};

std::ostream& operator << (std::ostream &stream, const Index &index);

bool operator == (const Index &index1, const Index &index2);

bool operator == (const Index &index1, const std::vector<int> &index2);

bool operator == (const std::vector<int> &index1, const Index &index2);

bool operator != (const Index &index1, const Index &index2);

bool operator != (const Index &index1, const std::vector<int> &index2);

bool operator != (const std::vector<int> &index1, const Index &index2);


#endif /*INDEX_HPP_*/
