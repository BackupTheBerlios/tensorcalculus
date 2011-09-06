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

#ifndef __TENSORTRAINREPRESENTATION_HPP
#define __TENSORTRAINREPRESENTATION_HPP

#include "TensorChainRepresentation.hpp"



namespace TensorCalculus {
	template <typename T> class TensorTrainRepresentation: public TensorChainRepresentation<T> {
		public:

			TensorTrainRepresentation(int r, const std::vector<int> &componentDimensions, const std::vector< std::vector<T> > &v) {
				int length = v.size();
#ifdef ARGUMENT_CHECKS_ON
				if (length < 3) {
					throw std::invalid_argument("The length has to be larger than 2.");
					/* actually, length=2 would also be possible, but this is a trivial case*/
				}
#endif
				std::vector<int> summations(length, r);

				summations[0] = 1;
				TensorChainRepresentation<T>::init(summations, componentDimensions, v);
			}
	};

	template<typename T>
	TensorTrainRepresentation<T> createRandomTensorTrain(int ttRank, const std::vector<int> &componentDimensions, T (*randomNumberGenerator)()) {
		int d = componentDimensions.size();

		std::vector< std::vector<T> > v_tt(d);

		v_tt[0].resize(componentDimensions[0]*ttRank);
		for (int n = 1; n < d-1; n++) {
			v_tt[n].resize(componentDimensions[n]*ttRank*ttRank);
		}
		v_tt[d-1].resize(componentDimensions[d-1]*ttRank);
		for (int n = 0; n < d; n++) {
			for (unsigned int k = 0; k < v_tt[n].size(); k++) {
				v_tt[n][k] = generator();
			}
		}

		TensorTrainRepresentation<T> tensorTrain(ttRank, componentDimensions, v_tt);

		return tensorTrain;
	}

	template<typename T>
	TensorTrainRepresentation<T> createRandomTensorTrain(int rank, int d, int componentDimension, T (*randomNumberGenerator)()) {
		std::vector<int> componentDimensions(d, componentDimension);

		return createRandomTensorTrain(rank, componentDimensions, randomNumberGenerator);
	}

}

#endif  //__TENSORTRAINREPRESENTATION_HPP
