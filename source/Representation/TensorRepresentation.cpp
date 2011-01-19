/*
 * Copyright (C) 2011 Stefan Handschuh
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

#include "Representation/TensorRepresentation.hpp"

namespace TensorCalculus {
std::vector<int> index2vector(const std::vector<int> &summations, int index) {
		int size = summations.size();

		std::vector<int> result(size);

		for (int n = 0; n < size; n++) {
			result[n] = index % summations[n];
			index /= summations[n];
		}
		return result;
	}

	std::vector<int> index2partialVector(const std::vector<int> &summations, int index,
			                             const std::vector<int> &incidenceRow) {
		int size = incidenceRow.size();

		std::vector<int> result(size);

		for (int n = 0; n < size; n++) {
			result[n] = index % summations[incidenceRow[n]];
			index /= summations[incidenceRow[n]];
		}
		return result;
	}

	int vector2index(const std::vector<int> &summations, const std::vector<int> &vector) {
#ifdef ARGUMENT_CHECKS_ON
		if (summations.size() != vector.size()) {
			std::invalid_argument("Vector count and summation count differ.");
		}
#endif
		int position = 0;

		int offset = 1;

		for (int n = 0, i = summations.size(); n < i; n++) {
			position += vector[n]*offset;
			offset *= summations[n];
		}
		return position;
	}

	int partialVector2index(const std::vector<int> &summations, const std::vector<int> &vector,
			                const std::vector<int> incidenceRow) {
		int size = incidenceRow.size();

		int position = 0;

		int offset = 1;

		if (vector.size() != size) {
			for (int n = 0; n < size; n++) {
				position += vector[incidenceRow[n]]*offset;
				offset *= summations[incidenceRow[n]];
			}
		} else {
			for (int n = 0; n < size; n++) {
				position += vector[n]*offset;
				offset *= summations[incidenceRow[n]];
			}
		}
		return position;
	}

	int partialVector2index(const std::vector<int> &summations, const std::vector<int> &vector,
			                const std::vector< std::vector<int> > &incidenceMatrix, int mu) {
		return partialVector2index(summations, vector, incidenceMatrix[mu]);
	}

	int partialIndex(const std::vector<int> &summations,
			         const std::vector< std::vector<int> > &incidenceMatrix, int index, int mu) {
		return partialVector2index(summations, index2vector(summations, index),
				                   incidenceMatrix, mu);
	}

	std::vector<int> getPartialSummations(const std::vector<int> &summations2,
			                              const std::vector<int> &incidenceRow2) {
		int n = incidenceRow2.size();

		std::vector<int> result(n);

		for (int i = 0; i < n; i++) {
			result[i] = summations2[incidenceRow2[i]];
		}
		return result;
	}

	std::vector<int> getPartialSummations(const std::vector<int> &summations,
			                              const std::vector<int> &incidenceRow, int k) {
		// k is the ignore-edge

		int n = incidenceRow.size()-1;

		std::vector<int> result(n);

		for (int i = 0, l = 0; i < n; i++, l++) {
			if (incidenceRow[i] == k) {
				++l;
			}
			result[i] = summations[incidenceRow[l]];
		}
		return result;
	}

	long getComponentProduct(const std::vector<int> &vector) {
		long result = 1;

		for (int n = 0, i = vector.size(); n < i; n++) {
			result *= vector[n];
		}
		return result;
	}

	void fillGenericVector(std::vector<int> &dest, int length, int offset, bool resize) {
		if (resize) {
			dest.resize(length);
		}
		for (int n = 0; n < length; n++) {
			dest[n] = n + offset;
		}
	}

	std::vector<int> genereateGenericVector(int length, int offset) {
		std::vector<int> result(length);

		fillGenericVector(result, length, offset, false);
		return result;
	}

}
