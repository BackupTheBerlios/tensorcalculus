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

#ifndef __CPTENSORREPRESENTATION_HPP
#define __CPTENSORREPRESENTATION_HPP

#include <iostream>
#include "TensorRepresentation.hpp"

namespace TensorCalculus {
	template <typename T> class CPTensorRepresentation: public TensorRepresentation<T> {
		public:
			CPTensorRepresentation(int r, int dimension, const std::vector< std::vector<T> > &v) {
				int d = v.size();
				
				std::vector<int> componentDimensions(d);
				
				for (int n = 0; n < d; n++) {
					componentDimensions[n] = dimension;
				}
				init(r, v, componentDimensions, std::vector<int>(0));
			}
			
			CPTensorRepresentation(int r, const std::vector<int> &componentDimensions,
					               const std::vector< std::vector<T> > &v) {
				init(r, v, componentDimensions, std::vector<int>(0));
			}

			CPTensorRepresentation(int r, int d, int componentDimension) {
				std::vector< std::vector<T> > v(d);

				std::vector<int> componentDimensions(d);

				for (int n = 0; n < d; n++) {
					componentDimensions[n] = componentDimension;
					v[n].resize(r*componentDimension);
				}
				init(r, v, componentDimensions, std::vector<int>(0));
			}

			CPTensorRepresentation(int r, int d, std::vector<int> &componentDimensions) {
				std::vector< std::vector<T> > v(d);

				for (int n = 0; n < d; n++) {
					v[n].resize(r*componentDimensions[n]);
				}
				init(r, v, componentDimensions, std::vector<int>(0));
			}


			CPTensorRepresentation(const char * filename) {
				TensorRepresentation<T>::init(filename);
#ifdef ARGUMENT_CHECKS_ON
				// check if this is really a CP-Tensor
				if (this->getL() > 0) {
					throw std::invalid_argument("Coefficient tensors are of invalid structure.");
				}
				if (this->summations.size() != 1) {
					throw std::invalid_argument("There has to be exactly one summation");
				}
				for (int n = 0; n < this->d; n++) {
					if (this->getIncidenceRow(n).size() != 1 || this->getIncidenceRow(n)[0] != 0) {
						throw std::invalid_argument("The given network structure is no "
													"CP-Tensor.");
					}
				}
#endif
			}

			CPTensorRepresentation(const std::vector<int> &componentDimensions,
					               const std::vector<T> &values) {
				long rank = componentProduct(componentDimensions);
#ifdef ARGUMENT_CHECKS_ON
				if (rank != values.size()) {
					throw std::invalid_argument("The number of data points does not match "
												"the dimension-product");
				}
#endif
				int d = componentDimensions.size();

				std::vector< std::vector<T> > v(d);

				for (int n = 0; n < d; n++) {
					v[n].resize(componentDimensions[n]*rank);
				}

				for (Index index(componentDimensions); !index.end(); ++index) {
					v[0][componentDimensions[0]*index + index[0]] = values[index];
					for (int n = 1; n < d; n++) {
						v[n][componentDimensions[n]*index + index[n]] = 1.0;
					}
				}
				// we should ballance here
				init(rank, v, componentDimensions, std::vector<int>(0));
			}

			void setValue(int mu, int j, int i, T value) {
				(*this).v[mu][j*(*this).componentDimensions[mu]+i] = value;
			}

			void setCrossApproximation(const CPTensorRepresentation<T> &representation) {
#ifdef ARGUMENT_CHECKS_ON
				this->checkCompatibility(representation);
				if (this->summations[0] >= representation.getSummation(0)) {
					throw std::invalid_argument("The rank of this representation is larger or "
												"equal to the given representation.");
				}
#endif

			}



	};

	template<typename T>
	CPTensorRepresentation<T> createRandomCPTensor(int rank, int d,
			                                       const std::vector<int> &componentDimensions,
			                                       T (*randomNumberGenerator)()) {
		std::vector< std::vector<T> > v_cp(d);
		for (int n = 0; n < d; n++) {
			v_cp[n].resize(rank*componentDimensions[n]);
			for (int k = 0, j = v_cp[n].size(); k < j; k++) {
				v_cp[n][k] = randomNumberGenerator();
			}
		}

		CPTensorRepresentation<T> cpTensor(rank, componentDimensions, v_cp);

		return cpTensor;
	}

	template<typename T>
	CPTensorRepresentation<T> createRandomCPTensor(int rank, int d, int componentDimension,
			                                       T (*randomNumberGenerator)()) {
		std::vector<int> componentDimensions(d);
		for (int n = 0; n < d; n++) {
			componentDimensions[n] = componentDimension;
		}
		return createRandomCPTensor(rank, d, componentDimensions, randomNumberGenerator);
	}

	template<typename T>
	CPTensorRepresentation<T> toCPTensorRepresentation(const TensorRepresentation<T> &representation) {
#ifdef ARGUMENT_CHECKS_ON
		if (representation.getL()) {
			throw std::invalid_argument("We do not want to increase d such that L has to be 0.");
		}
#endif
		Index index(representation.getSummations());

		int d = representation.getD(), rank = index.size();

		int TSIZE = sizeof(T);

		std::vector< std::vector<T> > v(d);

		for (int n = 0; n < d; n++) {
			v[n].resize(rank*representation.getComponentDimension(n));
		}

		for (index.begin(); !index.end(); ++index) {
			for (int n = 0; n < d; n++) {
				std::memcpy(&v[n][index*representation.getComponentDimension(n)],
							&representation[n][partialVector2index(representation.getSummations(),
									index, representation.getIncidenceRow(n))],
						    representation.getComponentDimension(n)*TSIZE);
			}
		}
		return CPTensorRepresentation<T>(rank,
				                         std::vector<int>(representation.getComponentDimensions()),
				                         v);
	}

}

#endif /* __CPTENSORREPRESENTATION_HPP */
