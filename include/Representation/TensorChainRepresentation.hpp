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

#ifndef __TENSORCHAINREPRESENTATION_HPP
#define __TENSORCHAINREPRESENTATION_HPP

#include "TensorRepresentation.hpp"
#include "CPTensorRepresentation.hpp"
#include "Matrix/MatrixOperators.hpp"

namespace TensorCalculus {

	template<typename T>
	struct dmrg_memory_t {
		std::vector<T> dummy;

		std::vector<T> a_ell_ell_plus_1;

		std::vector<T> u_ell_ell_plus_1;

		std::vector<T> workspace;

		std::vector<T> S;

		std::vector<T> U;

		std::vector<T> VT;

		std::vector<T> dummy_u;

		/* <a,a> */
		T norm_sqr_a;

		/* <u,u> */
		T norm_sqr_u;

		/* <u,u> - 2<a,u> */
		T norm_sqr_aou;

		minvert_result_t<T> minvert_result;
	};

	template<typename T>
	struct als_memory_t {
		std::vector<T> dummy_B_copy;

		std::vector<T> dummy;

		std::vector<T> dummy_v;

		/* <a,a> */
		T norm_sqr_a;

		/* <u,u> */
		T norm_sqr_u;

		/* <u,u> - 2<a,u> */
		T norm_sqr_aou;

		minvert_result_t<T> minvert_result;

	};

	template <typename T> class TensorChainRepresentation: public TensorRepresentation<T> {

		public:
			void init(const std::vector<int> &summations, const std::vector<int> &componentDimensions,
					  const std::vector< std::vector<T> > &v) {
				int length = v.size();

#ifdef ARGUMENT_CHECKS_ON
				if (length < 3) {
					throw std::invalid_argument("The length has to be larger than 2.");
				}
#endif
				std::vector< std::vector<int> > incidenceMatrix(length);

				for (int n = 0; n < length; n++) {
					incidenceMatrix[n].resize(2);
					incidenceMatrix[n][0] = n;
					incidenceMatrix[n][1] = (n+1) % length;
				}
				TensorRepresentation<T>::init(summations, v, componentDimensions, incidenceMatrix);
			}

		public:
			TensorChainRepresentation() {

			}

			TensorChainRepresentation(int r, int dimension, const std::vector< std::vector<T> > &v) {
				int length = v.size();

				std::vector<int> summations(length);
				std::vector<int> componentDimensions(length);
				for (int n = 0; n < length; n++) {
					summations[n] = r;
					componentDimensions[n] = dimension;
				}
				init(summations, componentDimensions, v);
			}
			
			TensorChainRepresentation(int r, const std::vector<int> &componentDimensions,
					                  const std::vector< std::vector<T> > &v) {
				int length = v.size();

				std::vector<int> summations(length);
				for (int n = 0; n < length; n++) {
					summations[n] = r;
				}
				init(summations, componentDimensions, v);
			}

			TensorChainRepresentation(const std::vector<int> &summations,
					                  const std::vector<int> &componentDimensions,
					                  const std::vector< std::vector<T> > &v) {
				init(summations, componentDimensions, v);
			}
			
			TensorChainRepresentation(const CPTensorRepresentation<T> &cpTensor) {
				int TSIZE = sizeof(T);

				int length = cpTensor.getV().size();

				std::vector< std::vector<T> > v(length);
				std::vector<int> componentDimensions = cpTensor.getComponentDimensions();

				int summation = cpTensor.getSummation(0);

				std::vector<int> summations(length, summation);

				for (int n = 0; n < length-2; n++) {
					v[n].resize(componentDimensions[n]*summations[n]*summations[n]);
					for (int i = 0; i < summations[n]; i++) {
						std::memcpy(&v[n][componentDimensions[n]*(i*summations[n]+i)],
								    &cpTensor[n][componentDimensions[n]*i],
								    TSIZE*componentDimensions[n]);
					}
				}
				v[length-2].resize(componentDimensions[length-2]*summation);
				std::memcpy(&v[length-2][0], &cpTensor[length-2][0],
							TSIZE*componentDimensions[length-2]*summation);
				summations[length-1] = 1;
				v[length-1].resize(componentDimensions[length-1]*summation);
				std::memcpy(&v[length-1][0], &cpTensor[length-1][0],
						    TSIZE*componentDimensions[length-1]*summation);

				init(summations, std::vector<int>(componentDimensions), v);
			}

			TensorChainRepresentation(const char * filename) {
				TensorRepresentation<T>::init(filename);

#ifdef ARGUMENT_CHECKS_ON
				/* now, we check for TC-structure */
				if (this->d != this->getSummationsCount()) {
					throw std::invalid_argument("The file does not contain a TC.");
				}

				/* take the first node as reference */
				std::vector<int> swapSummations(2);

				swapSummations[0] = this->getIncidenceRow(0)[0];
				swapSummations[1] = this->getIncidenceRow(0)[1];
				while(swapSummations.size() < this->getSummationsCount()) {
					for (int n = 1; n < this->d; n++) {
						if (swapSummations[swapSummations.size() - 1] == this->getIncidenceRow(n)[0]) {
							swapSummations.push_back(this->getIncidenceRow(n)[1]);
							break;
						}
					}
				}
				//for (int n = 0; n < this->getSummationsCount())
				using namespace VectorOperators;
				std::cout << swapSummations << std::endl;
#endif
			}

			/**
			 * Returns the factor-matrix in the given direction for H
			 */
			std::vector<T> generate_h_matrix(int direction)
			{
				using namespace VectorOperators;

				int dimension = this->componentDimensions[direction];

				std::vector<int> summations = this->getSummations(direction);

				std::vector<T> result(sqr(componentProduct(summations)));

				for (Index index1(summations); !index1.end(); ++index1) {
					for (Index index2(summations); !index2.end(); ++index2) {
						const int column = index2.getCurrent()[0]*summations[1] + index2.getCurrent()[1];

						const int row = index1.getCurrent()[0]*summations[0] + index1.getCurrent()[1];

						// TODO: use the symmetric structure of H
						result[column*sqr(summations[0]) + row] = this->componentScalarProduct(direction, index1,
																							   &this->v[direction][index2*dimension],
																							   dimension) ;// store the matrix colum wise
					}
				}
				return result;
			}

			void fill_A_mu(const TensorChainRepresentation<T> &representation, int mu, std::vector<T> &destination) const {
				int mu_plus_1 = (mu + 1) % this->d;

				int R_mu = representation.summations[mu],
					R_mu_plus_1 = representation.summations[mu_plus_1],
					r_mu = this->summations[mu],
					r_mu_plus_1 = this->summations[mu_plus_1];

				int dimension = this->componentDimensions[mu];

				if (destination.size() < static_cast<unsigned int>(R_mu*R_mu_plus_1*r_mu*r_mu_plus_1)) {
					destination.resize(R_mu*R_mu_plus_1*r_mu*r_mu_plus_1);
				}

				for (int n = 0; n < R_mu; n++) {
					for (int m = 0; m < R_mu_plus_1; m++) {
						for (int k = 0; k < r_mu; k++) {
							const int row = n+k*R_mu;

							for (int l = 0; l < r_mu_plus_1; l++) {
								const int column = m+l*R_mu_plus_1;

								destination[row+column*R_mu*r_mu] = Blas<T>::dot(dimension,
										&(representation.getV()[mu][dimension*(n+m*R_mu)]), 1, &(this->v[mu][dimension*(k+l*r_mu)]), 1);
							}
						}
					}
				}
			}

			void fill_B_mu(int mu, std::vector<T> &destination) const {
				fill_A_mu(*this, mu, destination);
			}

			/*
			 * Turns a R_mu2 r_mu2 x R_mu1 r_mu1 into a
			 *         r_mu1 r_mu2 x R_mu1 R_mu2 matrix.
			 */
			void reshape_into(const TensorChainRepresentation<T> &representation, const int mu1,
					          const int mu2,
					          const std::vector<T> &source, std::vector<T> &destination) const {

				int R_mu1 = representation.summations[mu1],
				    R_mu2 = representation.summations[mu2],
					r_mu1 = this->summations[mu1],
					r_mu2 = this->summations[mu2];

				//std::cout << "reshape von " << R_mu2* r_mu2 << "x" << R_mu1*r_mu1 << " nach "  << r_mu1*r_mu2 <<"x"<< R_mu1*R_mu2 << std::endl;

				if (destination.size() < static_cast<unsigned int>(R_mu1*R_mu2*r_mu1*r_mu2)) {
					destination.resize(R_mu1*R_mu2*r_mu1*r_mu2);
				}
				if (source.size() < static_cast<unsigned int>(R_mu1*R_mu2*r_mu1*r_mu2)) {
					throw std::invalid_argument("source has wrong size.");
				}

				for (int m = 0; m < R_mu2; m++) {
					for (int l = 0; l < r_mu2; l++) {
						const int row_s = l*R_mu2 + m;

						for (int n = 0; n < R_mu1; n++) {
							const int column_d = m*R_mu1+n;

							for (int k = 0; k < r_mu1; k++) {
								const int column_s = k*R_mu1+n,
										  row_d = l*r_mu1 + k;

								destination[column_d*r_mu1*r_mu2 + row_d] =
										source[column_s*R_mu2*r_mu2 + row_s];
							}

						}
					}
				}
			}


			void reshape_into(const int mu1, const int mu2, const std::vector<T> &source,
					           std::vector<T> &destination) const {
				reshape_into(*this, mu1, mu2, source, destination);
			}

			void reshape_vector_indices(std::vector<T>& vector, const std::vector<int> &indices, const int dimension) {
				int size = vector.size();

				int TSIZE = sizeof(T);

				std::vector<int> min(indices);
				std::sort(min.begin(), min.end());

				for (int n = 0, i = indices.size(); n < i; n++) {
					const int j = min[n];

					std::memmove(&vector[(j+1)*dimension], &vector[j*dimension], (size-j-1)*TSIZE*dimension);
					Blas<T>::scal(dimension, 0.0, &vector[j*dimension], 1);
				}
			}

			void performSingleALSStep(const TensorChainRepresentation<T> &representation,
									  const int mu, std::vector<T> dummy_B,
									  std::vector<T> dummy_A, als_memory_t<T> * als_memory,
									  std::vector<int>& dependentIndices) {
				const int mu_plus_1 = (mu + 1) % this->d;

				std::memcpy(&als_memory->dummy_B_copy[0], &dummy_B[0],
						    this->TSIZE*sqr(this->summations[mu]*this->summations[mu_plus_1]));

				if (invert_gramian_matrix(dummy_B, this->summations[mu]*this->summations[mu_plus_1], dependentIndices) == 0) {
					// dummy_B = dummy_B^{-1}

					delete_rows(dummy_A, this->summations[mu]*this->summations[mu_plus_1],
							    representation.summations[mu]*representation.summations[mu_plus_1], dependentIndices);

					const int rankReduction = dependentIndices.size();

					Blas<T>::gemm('N', 'T', this->componentDimensions[mu],
								  this->summations[mu]*this->summations[mu_plus_1]-rankReduction,
								  representation.summations[mu]*representation.summations[mu_plus_1],
								  1, &(representation.v[mu][0]), &dummy_A[0],
								  0, &als_memory->dummy[0]); //dummy = representation.v[mu]*dummy_A^T;


						Blas<T>::gemm('N', 'T', this->componentDimensions[mu],
									  this->summations[mu]*this->summations[mu_plus_1]-rankReduction,
									  this->summations[mu]*this->summations[mu_plus_1]-rankReduction,
									  1, &als_memory->dummy[0], &dummy_B[0], 0, &(this->v[mu][0]));	//this->v[mu] = dummy*dummy_B^T;

						als_memory->norm_sqr_aou = Blas<T>::dot(this->componentDimensions[mu]*(this->summations[mu]*this->summations[mu_plus_1]-rankReduction),
														        &(this->v[mu][0]), 1, &als_memory->dummy[0], 1);

						/* move the entries in v[mu] to the correct position (considering the zeros)*/
						//std::cout << this->v[mu] << std::endl;
						reshape_vector_indices(this->v[mu], dependentIndices, this->componentDimensions[mu]);

						Blas<T>::gemm('N', 'T', this->componentDimensions[mu],
									  this->summations[mu]*this->summations[mu_plus_1],
									  this->summations[mu]*this->summations[mu_plus_1],
									  1, &(this->v[mu][0]), &als_memory->dummy_B_copy[0],
									  0, &als_memory->dummy[0]); //als_memory->dummy = this->v[mu]*dummy_B_copy^T;



						als_memory->norm_sqr_u = Blas<T>::dot(this->componentDimensions[mu]*this->summations[mu]*this->summations[mu_plus_1],
															  &(this->v[mu][0]), 1, &als_memory->dummy[0], 1);


					//std::cout << "norm_after=" << als_memory->norm_sqr_u - 2* als_memory->norm_sqr_aou + als_memory->norm_sqr_a << std::endl;
				} else {
					std::cout << "SINGULAR" << std::endl;
				}
			}

			/*
			 * Complexity n*r^2, has to be performed only one time, independent
			 * of the number of iterations.
			 */
			void detemine0vectors(int mu, std::vector<int> & dependentIndices) {
				int mu_plus_1 = (mu+1) % this->d;

				int dimension = this->componentDimensions[mu];

				int r_mu1 = this->summations[mu],
					r_mu2 = this->summations[mu_plus_1];

				for (int j = 0; j < r_mu2; j++) {
					for (int i = 0; i < r_mu1; i++) {
						bool zero = true;

						const int index = j*r_mu1+i;

						for (int n = 0; n < dimension; n++) {
							if (this->v[mu][dimension*index+n] != 0.0) {
								zero = false;
								break;
							}
						}
						if (zero) {
							dependentIndices.push_back(index);
						}
					}
				}
			}

			int performALS(const TensorChainRepresentation<T> &representation, T eps = 1e-4) {
#ifdef ARGUMENT_CHECKS_ON
				checkCompatibility(representation);
#endif
				// we assume that the summations of this are smaller than the summations of representation

				int r_max_sqr = getSubsequentMaxRankProduct(*this);
				int R_max_sqr = getSubsequentMaxRankProduct(representation);

				als_memory_t<T> als_memory;

				als_memory.dummy_B_copy.resize(r_max_sqr);
				als_memory.dummy_v.resize(std::sqrt(r_max_sqr)*this->getMaxComponentDimension());

				std::vector< std::vector<int> > dependentIndices(this->d);

				/* determine 0-vectors*/
				for (int n = 0; n < this->d; n++) {
					//detemine0vectors(n, &linear_dependences[n]);
				}

				/* allocate memory that is needed for the method */
				/* dummy variables */
				std::vector<T> dummy_A(R_max_sqr),
							   dummy_B(r_max_sqr);

				std::vector< std::vector<T> > A_pre(this->d-1);
				std::vector< std::vector<T> > B_pre(this->d-1);

				A_pre[this->d-2].resize(representation.summations[this->d-1]*this->summations[this->d-1]*
						                representation.summations[0]*this->summations[0]);
				B_pre[this->d-2].resize(sqr(this->summations[this->d-1])*sqr(this->summations[0]));

				std::vector<T> A_tilde(this->getMaxRank(1)*representation.getMaxRank(1)*this->summations[0]*representation.summations[0]);
				std::vector<T> B_tilde(sqr(this->getMaxRank(1)*this->summations[0]));

				int dimProduct = this->componentDimensions[0]*representation.summations[0]*representation.summations[1];

				for (int mu = 0; mu < this->d; mu++) {
					const int mu_plus_1 = (mu + 1) % this->d;

					const int buffer = this->componentDimensions[mu]*representation.summations[mu]*representation.summations[mu_plus_1];
					if (buffer > dimProduct) {
						dimProduct = buffer;
					}
				}
				als_memory.dummy.resize(std::max(std::max(R_max_sqr, dimProduct), (int) std::max(A_tilde.size(), B_tilde.size())));

				for (int ell = this->d -3; ell >-1; ell--) {
					A_pre[ell].resize(representation.summations[ell+1]*this->summations[ell+1]*
									  representation.summations[0]*this->summations[0]);
					B_pre[ell].resize(sqr(this->summations[0])*sqr(this->summations[ell+1]));
				}
				/* memory allocation end */

				als_memory.norm_sqr_a = representation.l2norm_sqr();
				als_memory.norm_sqr_u = l2norm_pre_sqr(B_pre[0], B_tilde);
				als_memory.norm_sqr_aou = l2norm_pre_sqr(representation, A_pre[0], A_tilde);

				T norm_dist = std::sqrt(std::abs((als_memory.norm_sqr_u + als_memory.norm_sqr_a - 2*als_memory.norm_sqr_aou)/als_memory.norm_sqr_a)),
				  norm_dist0 = 0.0;

				int stepcount = 0;

				while (std::abs(norm_dist - norm_dist0) > eps) {
					stepcount++;
					norm_dist0 = norm_dist;
					/* preprocess begin */
					fill_A_mu(representation, this->d-1, A_pre[this->d-2]);
					fill_B_mu(this->d-1, B_pre[this->d-2]);
					for (int ell = this->d -3; ell >-1; ell--) {
						fill_A_mu(representation, ell+1, dummy_A);
						Blas<T>::gemm('N', 'N', representation.summations[ell+1]*this->summations[ell+1],
									  representation.summations[0]*this->summations[0],
									  representation.summations[ell+2]*this->summations[ell+2],
									  1.0, &dummy_A[0], &A_pre[ell+1][0], 0.0, &A_pre[ell][0]); //A_pre[ell] = dummy_A * A_pre[ell+1];
						fill_B_mu(ell+1, dummy_B);
						Blas<T>::gemm('N', 'N', sqr(this->summations[ell+1]),
									  sqr(this->summations[0]),
									  sqr(this->summations[ell+2]),
									  1, &dummy_B[0], &B_pre[ell+1][0], 0, &B_pre[ell][0]); //B_pre[ell] = dummy_B * B_pre[ell+1];
					}
					/* preprocess end */

					reshape_into(representation, 0, 1, A_pre[0], dummy_A); //dummy_A = reshape(A_pre[0]));
					reshape_into(0, 1, B_pre[0], dummy_B); //dummy_B = reshape(B_pre[0]));

					performSingleALSStep(representation, 0, dummy_B, dummy_A, &als_memory, dependentIndices[0]);

					fill_A_mu(representation, 0, A_tilde);
					fill_B_mu(0, B_tilde);

					for (int ell = 1; ell < this->d-1; ell++) {
						Blas<T>::gemm('N', 'N', representation.summations[ell+1]*this->summations[ell+1],
									  representation.summations[ell]*this->summations[ell],
									  representation.summations[0]*this->summations[0],
									  1, &(A_pre[ell][0]), &A_tilde[0], 0, &als_memory.dummy[0]); // dummy = A_pre[ell]*A_tilde
						reshape_into(representation, ell, ell+1, als_memory.dummy, dummy_A); //dummy_A = reshape(dummy);

						Blas<T>::gemm('N', 'N', sqr(this->summations[ell+1]),
									  sqr(this->summations[ell]),
									  sqr(this->summations[0]),
									  1, &(B_pre[ell][0]), &B_tilde[0], 0, &als_memory.dummy[0]); // dummy = B_pre[ell]*B_tilde;

						reshape_into(ell, ell+1, als_memory.dummy, dummy_B); //dummy_B = reshape(dummy)
						performSingleALSStep(representation, ell, dummy_B, dummy_A, &als_memory, dependentIndices[ell]);

						// now, this step is done and we calculate the stuff for the next step

						fill_A_mu(representation, ell, dummy_A);
						std::memcpy(&als_memory.dummy[0], &A_tilde[0], this->TSIZE*this->summations[0]*this->summations[ell]*
									representation.summations[0]*representation.summations[ell]);	// dummy = A_tilde
						Blas<T>::gemm('N', 'N', representation.summations[0]*this->summations[0],
									  representation.summations[ell+1]*this->summations[ell+1],
									  representation.summations[ell]*this->summations[ell],
									  1, &als_memory.dummy[0],
									  &dummy_A[0], 0, &A_tilde[0]);//A_tilde = dummy*dummy_A;
						fill_B_mu(ell, dummy_B);
						std::memcpy(&als_memory.dummy[0], &B_tilde[0], this->TSIZE*sqr(this->summations[0])*
									sqr(this->summations[ell])); // dummy = B_tilde
						Blas<T>::gemm('N', 'N', sqr(this->summations[0]),
									  sqr(this->summations[ell+1]),
									  sqr(this->summations[ell]), 1, &als_memory.dummy[0],
									  &dummy_B[0], 0, &B_tilde[0]); //B_tilde = dummy*dummy_B;
					}

					reshape_into(representation, this->d-1, 0, A_tilde, dummy_A); //dummy_A = reshape(A_tilde);
					reshape_into(this->d-1, 0, B_tilde, dummy_B); //dummy_B = reshape(B_tilde));
					performSingleALSStep(representation, this->d-1, dummy_B, dummy_A, &als_memory, dependentIndices[this->d-1]);
					//norm_sqr_u = l2norm_post_sqr(B_tilde, dummy_B);
					//norm_sqr_aou = l2norm_post_sqr(representation, A_tilde, dummy_A);
					norm_dist = std::sqrt((als_memory.norm_sqr_u +als_memory.norm_sqr_a - 2*als_memory.norm_sqr_aou) / als_memory.norm_sqr_a);
					std::cout << stepcount << "\t " << norm_dist << std::endl;
				}
				return stepcount;
			}

			void performSingleDMRGStep(const TensorChainRepresentation &representation,
					                   const int mu, const std::vector<T> &dummy_A,
					                   std::vector<T> &dummy_B, dmrg_memory_t<T> * dmrg_memory,
					                   const T eps) {

				const int mu_plus_1 = (mu + 1) % this->d,
						  mu_plus_2 = (mu + 2) % this->d;

				int dim_prod = this->componentDimensions[mu]* this->componentDimensions[mu_plus_1];

				T partialNormSumSQR = Blas<T>::asum(this->summations[mu]*this->summations[mu_plus_2],
												    &dummy_B[0],
													this->summations[mu]*this->summations[mu_plus_2]+1);
				if (dummy_B.size() < static_cast<unsigned int>(sqr(this->summations[mu]*this->summations[mu_plus_2]))) {
					throw std::invalid_argument("dummy_B does not have the required size.");
				}

				std::vector<T> dummy_B_copy(dummy_B);

				if (invert_matrix(dummy_B, this->summations[mu]*this->summations[mu_plus_2], &(dmrg_memory->minvert_result)) != 0) {// dummy_B = dummy_B^{-1}
					std::cout << "SINGULAR" << std::endl;
					return;
				}
				fill_a_mu_mu_plus_1(representation, mu, dmrg_memory->a_ell_ell_plus_1); // generate a_{0,1} \in R^{n_0*n_1 x R_0*R_2}
				dmrg_memory->dummy.resize(dim_prod*
										  this->summations[mu]*this->summations[mu_plus_2]);

				Blas<T>::gemm('N', 'T', dim_prod,
						      this->summations[mu]*this->summations[mu_plus_2],
						      representation.summations[mu]*representation.summations[mu_plus_2],
						      1, &(dmrg_memory->a_ell_ell_plus_1[0]), &dummy_A[0], 0, &(dmrg_memory->dummy[0])); // dummy = a_ell_ell_plus_1*dummy_A^T
				dmrg_memory->u_ell_ell_plus_1.resize(this->componentDimensions[mu]*this->componentDimensions[mu_plus_1]*
				                                     this->summations[mu]*this->summations[mu_plus_2]);

				Blas<T>::gemm('N', 'T', dim_prod,
							  this->summations[mu]*this->summations[mu_plus_2],
							  this->summations[mu]*this->summations[mu_plus_2],
							  1, &(dmrg_memory->dummy[0]), &dummy_B[0], 0, &(dmrg_memory->u_ell_ell_plus_1[0])); // u_ell_ell_plus_1 = dummy*dummy_B^T

				// now, we can test wheather the optimization was successful or not

				T norm_sqr_aou = Blas<T>::dot(dim_prod*this->summations[mu]*this->summations[mu_plus_2],
										      &dmrg_memory->u_ell_ell_plus_1[0], 1,
										      &dmrg_memory->dummy[0], 1);

				Blas<T>::gemm('N', 'T', dim_prod,
										this->summations[mu]*this->summations[mu_plus_2],
										this->summations[mu]*this->summations[mu_plus_2],
										1, &(dmrg_memory->u_ell_ell_plus_1[0]), &dummy_B_copy[0], 0, &(dmrg_memory->dummy[0])); // dummy = a_ell_ell_plus_1*dummy_A^T
				T norm_sqr_u = Blas<T>::dot(dim_prod*this->summations[mu]*this->summations[mu_plus_2],
									        &dmrg_memory->u_ell_ell_plus_1[0], 1,
									        &dmrg_memory->dummy[0], 1);

				if (dmrg_memory->minvert_result.approx) {
					/*
					T sum1 = 0, sum2 = 0;

					using namespace VectorOperators;
					for (int e = 0; e < this->summations[mu]; e++) {
						for (int r = 0; r < this->summations[mu_plus_2]; r++) {

							std::vector<T> test(dim_prod), buffer(dim_prod);
							for (int q = 0; q < representation.summations[mu]; q++) {
								for (int w = 0; w < representation.summations[mu_plus_2]; w++) {
									std::memcpy(&buffer[0], &(dmrg_memory->a_ell_ell_plus_1[(w*representation.summations[mu]+q)*dim_prod]), this->TSIZE * dim_prod);

									sum1 += Blas<T>::dot(dim_prod,
														 &(dmrg_memory->a_ell_ell_plus_1[(w*representation.summations[mu]+q)*dim_prod]),1,
														 &(dmrg_memory->u_ell_ell_plus_1[(r*this->summations[mu]+e)*dim_prod]), 1)*
														 dummy_A[(w*representation.summations[mu]+q)*this->summations[mu]*this->summations[mu_plus_2] + (r*this->summations[mu]+e)];
									buffer *= -dummy_A_copy[(w*representation.summations[mu]+q)*this->summations[mu]*this->summations[mu_plus_2] + (r*this->summations[mu]+e)];
									test += buffer;
								}
							}
							for (int q = 0; q < this->summations[mu]; q++) {
								for (int w = 0; w < this->summations[mu_plus_2]; w++) {
									std::memcpy(&buffer[0], &(dmrg_memory->u_ell_ell_plus_1[(w*this->summations[mu]+q)*dim_prod]), this->TSIZE * dim_prod);

									sum2 += Blas<T>::dot(dim_prod,
														 &(dmrg_memory->u_ell_ell_plus_1[(w*this->summations[mu]+q)*dim_prod]), 1,
														 &(dmrg_memory->u_ell_ell_plus_1[(r*this->summations[mu]+e)*dim_prod]), 1)*
														 dummy_B_copy[(w*this->summations[mu]+q)*this->summations[mu]*this->summations[mu_plus_2] + (r*this->summations[mu]+e)];

									buffer *= dummy_B_copy[(w*this->summations[mu]+q)*this->summations[mu]*this->summations[mu_plus_2] + (r*this->summations[mu]+e)];
									test += buffer;
								}
							}

							std::cout << " " << l2_norm(test);
						}


					}
					norm_sqr_u = sum2;
					norm_sqr_aou = sum1;
					*/

					if (std::abs(norm_sqr_u - 2*norm_sqr_aou + dmrg_memory->norm_sqr_a) >= std::abs(dmrg_memory->norm_sqr_u - 2* dmrg_memory->norm_sqr_aou + dmrg_memory->norm_sqr_a)) {
						std::cout << "new norm would be " << std::sqrt(std::abs(norm_sqr_u - 2*norm_sqr_aou + dmrg_memory->norm_sqr_a)/dmrg_memory->norm_sqr_a) << std::endl;
						return;
					}
				}
				std::cout << "new norm should be " << std::sqrt(std::abs(norm_sqr_u - 2*norm_sqr_aou + dmrg_memory->norm_sqr_a)/dmrg_memory->norm_sqr_a) << std::endl;
				dmrg_memory->norm_sqr_u = norm_sqr_u;
				dmrg_memory->norm_sqr_aou = norm_sqr_aou;
				reshape_u_ell_ell_plus_1(*this, mu, dmrg_memory->u_ell_ell_plus_1, dmrg_memory->dummy_u); // reshape u_ell_ell_plus_1 into dummy_u

				int orig_new_rank = std::min(this->summations[mu]*this->componentDimensions[mu],
											 this->summations[mu_plus_2]*this->componentDimensions[mu_plus_1]);

				dmrg_memory->S.resize(orig_new_rank);
				dmrg_memory->U.resize(orig_new_rank*this->summations[mu]*this->componentDimensions[mu]);
				dmrg_memory->VT.resize(orig_new_rank*this->summations[mu_plus_2]*this->componentDimensions[mu_plus_1]);

				T workspaceLength = 0.0;

				Lapack<T>::gesvd('S', 'S', this->summations[mu]*this->componentDimensions[mu],
						         this->summations[mu_plus_2]*this->componentDimensions[mu_plus_1],
						         0, this->summations[mu]*this->componentDimensions[mu],
						         0, 0, this->summations[mu]*this->componentDimensions[mu], 0,
						         orig_new_rank, &workspaceLength, -1);

				dmrg_memory->workspace.resize((int) workspaceLength);

				Lapack<T>::gesvd('S', 'S', this->summations[mu]*this->componentDimensions[mu],
						         this->summations[mu_plus_2]*this->componentDimensions[mu_plus_1],
						         &(dmrg_memory->dummy_u[0]), this->summations[mu]*this->componentDimensions[mu],
						         &(dmrg_memory->S[0]), &(dmrg_memory->U[0]), this->summations[mu]*this->componentDimensions[mu],
						         &(dmrg_memory->VT[0]), orig_new_rank, &(dmrg_memory->workspace[0]), (int) workspaceLength);
									// SVD of the dummy_u

				int rank = orig_new_rank;

				T partialSVDNormSQR = 0.0;

				for (int n = 1; n < rank; n++) {
					partialSVDNormSQR += dmrg_memory->S[rank-n];
					if (partialNormSumSQR * partialSVDNormSQR > sqr(eps)) {
						rank -= n -1;
						break;
					}
				}

				/* adjusts rank and sets new values */
				apply_svd_result(*this, dmrg_memory->S, dmrg_memory->U, dmrg_memory->VT, rank, mu);
			}

			void performDMRG(const TensorChainRepresentation<T> &representation, T eps = 10e-3) {
#ifdef ARGUMENT_CHECKS_ON
				checkCompatibility(representation);
#endif
				int d = this->d;

				/* dummy variables */
				std::vector<T> dummy_A;
				std::vector<T> dummy_B;

				dmrg_memory_t<T> dmrg_memory;

				dmrg_memory.norm_sqr_a = representation.l2norm_sqr();

				std::vector<T> A_tilde(this->summations[1]*this->getMaxRank(2)*representation.summations[1]*representation.getMaxRank(2));
				std::vector<T> B_tilde(sqr(this->summations[1]*this->getMaxRank(2)));


				dmrg_memory.dummy.resize(std::max(A_tilde.size(), B_tilde.size()));

				if (representation.summations[0] > 1) { // so in case of TT, this step is omitted
					dummy_A.resize(getSubsequentMaxRankProduct(representation));
					dummy_B.resize(getSubsequentMaxRankProduct(*this));

					fill_B_mu(1, B_tilde);
					fill_A_mu(representation, 1, A_tilde);
					for (int mu = 2; mu < d-1; mu++) {
						fill_B_mu(mu, dummy_B);
						Blas<T>::gemm('N', 'N', sqr(this->summations[1]),
									  sqr(this->summations[mu+1]),
									  sqr(this->summations[mu]),
									  1.0, &B_tilde[0], &dummy_B[0], 0.0, &(dmrg_memory.dummy[0]));
						std::memcpy(&B_tilde[0], &(dmrg_memory.dummy[0]),
								    this->TSIZE*sqr(this->summations[1])*sqr(this->summations[mu+1]));
						fill_A_mu(representation, mu, dummy_A);
						Blas<T>::gemm('N', 'N', representation.summations[1]*this->summations[1],
									  representation.summations[mu+1]*this->summations[mu+1],
									  representation.summations[mu]*this->summations[mu],
									  1.0, &A_tilde[0], &dummy_A[0], 0.0, &(dmrg_memory.dummy[0]));
						std::memcpy(&A_tilde[0], &(dmrg_memory.dummy[0]),
									this->TSIZE*representation.summations[1]*this->summations[1]*
									representation.summations[mu+1]*this->summations[mu+1]);
					}
					reshape_into(representation, d-1, 1, A_tilde, dummy_A); //dummy_A = reshape(A_tilde);
					reshape_into(d-1, 1, B_tilde, dummy_B); //dummy_B = reshape(B_tilde)
					// compute norm here
					performSingleDMRGStep(representation, d-1, dummy_A, dummy_B, &dmrg_memory, eps);

					std::cout << std::sqrt(std::abs((dmrg_memory.norm_sqr_u+dmrg_memory.norm_sqr_a-2*dmrg_memory.norm_sqr_aou)/dmrg_memory.norm_sqr_a)) << " vs " << rel_l2dist(representation) << std::endl;

					dummy_A.resize(std::max(representation.summations[0]*representation.getMaxRank(1), getSubsequentMaxRankProduct(representation)));
					dummy_B.resize(std::max(this->summations[0]*this->getMaxRank(1), getSubsequentMaxRankProduct(*this)));
				}



				/* preprocess begin */
				std::vector< std::vector<T> > A_pre(d-2);
				std::vector< std::vector<T> > B_pre(d-2);
				A_pre[d-3].resize(representation.summations[d-1]*this->summations[d-1]*
						          representation.summations[0]*this->summations[0]);
				fill_A_mu(representation, d-1, A_pre[d-3]);
				B_pre[d-3].resize(sqr(this->summations[d-1])*sqr(this->summations[0]));
				fill_B_mu(d-1, B_pre[d-3]);
				for (int ell = d -3; ell > 0; ell--) {
					fill_A_mu(representation, ell+1, dummy_A);
					A_pre[ell-1].resize(representation.summations[ell+1]*this->summations[ell+1]*
							            representation.summations[0]*this->summations[0]);
					Blas<T>::gemm('N', 'N', representation.summations[ell+1]*this->summations[ell+1],
							      representation.summations[0]*this->summations[0],
							      representation.summations[ell+2]*this->summations[ell+2],
							      1.0, &dummy_A[0], &A_pre[ell][0], 0.0, &A_pre[ell-1][0]); //A_pre[ell] = dummy_A * A_pre[ell+1];
					fill_B_mu(ell+1, dummy_B);
					B_pre[ell-1].resize(sqr(this->summations[ell+1])*sqr(this->summations[0]));
					Blas<T>::gemm('N', 'N', sqr(this->summations[ell+1]),
								  sqr(this->summations[0]),
								  sqr(this->summations[ell+2]),
								  1, &dummy_B[0], &B_pre[ell][0], 0, &B_pre[ell-1][0]); //B_pre[ell] = dummy_B * B_pre[ell+1];
				}
				/* preprocess end */

				reshape_into(representation, 0, 2, A_pre[0], dummy_A); //dummy_A = reshape(A_pre[0]));
				reshape_into(0, 2, B_pre[0], dummy_B); //dummy_B = reshape(B_pre[0]));

				performSingleDMRGStep(representation, 0, dummy_A, dummy_B, &dmrg_memory, eps);

				std::cout << std::sqrt(std::abs((dmrg_memory.norm_sqr_u+dmrg_memory.norm_sqr_a-2*dmrg_memory.norm_sqr_aou)/dmrg_memory.norm_sqr_a)) << " vs " << rel_l2dist(representation) << std::endl;

				A_tilde.resize(this->summations[0]*this->summations[1]*representation.summations[0]*representation.summations[1]);
				B_tilde.resize(sqr(this->summations[0]*this->summations[1]));

				fill_A_mu(representation, 0, A_tilde);
				fill_B_mu(0, B_tilde);
				for (int mu = 1; mu < d -2; mu++) {
					const int mu_plus_2 = (mu + 2) % d;

					dmrg_memory.dummy.resize(representation.summations[mu_plus_2]*this->summations[mu_plus_2]*
								              representation.summations[mu]*this->summations[mu]);
					Blas<T>::gemm('N', 'N', representation.summations[mu_plus_2]*this->summations[mu_plus_2],
							      representation.summations[mu]*this->summations[mu],
							      representation.summations[0]*this->summations[0],
								  1, &(A_pre[mu][0]), &A_tilde[0], 0, &(dmrg_memory.dummy[0])); // dummy = A_pre[ell]*A_tilde
					dummy_A.resize(dmrg_memory.dummy.size());
					reshape_into(representation, mu, mu_plus_2, dmrg_memory.dummy, dummy_A); //dummy_A = reshape(dummy);
					dmrg_memory.dummy.resize(sqr(this->summations[mu_plus_2])*sqr(this->summations[mu]));
					Blas<T>::gemm('N', 'N', sqr(this->summations[mu_plus_2]),
								  sqr(this->summations[mu]),
								  sqr(this->summations[0]),
								  1, &(B_pre[mu][0]), &B_tilde[0], 0, &(dmrg_memory.dummy[0])); // dummy = B_pre[ell]*B_tilde;
					dummy_B.resize(dmrg_memory.dummy.size());
					reshape_into(mu, mu_plus_2, dmrg_memory.dummy, dummy_B); //dummy_B = reshape(dummy)

					performSingleDMRGStep(representation, mu, dummy_A, dummy_B, &dmrg_memory, eps);
					std::cout << std::sqrt(std::abs((dmrg_memory.norm_sqr_u+dmrg_memory.norm_sqr_a-2*dmrg_memory.norm_sqr_aou)/dmrg_memory.norm_sqr_a)) << " vs " << rel_l2dist(representation) << std::endl;

					dummy_A.resize(this->summations[mu]*representation.summations[mu]*
							       this->summations[mu+1]*representation.summations[mu+1]);
					fill_A_mu(representation, mu, dummy_A);
					dmrg_memory.dummy.resize(A_tilde.size());
					std::memcpy(&(dmrg_memory.dummy[0]), &A_tilde[0], this->TSIZE*this->summations[0]*this->summations[mu]*
								representation.summations[0]*representation.summations[mu]);	// dummy = A_tilde
					A_tilde.resize(representation.summations[0]*this->summations[0]*
							       representation.summations[mu+1]*this->summations[mu+1]);
					Blas<T>::gemm('N', 'N', representation.summations[0]*this->summations[0],
								  representation.summations[mu+1]*this->summations[mu+1],
								  representation.summations[mu]*this->summations[mu],
								  1, &(dmrg_memory.dummy[0]),
								  &dummy_A[0], 0, &A_tilde[0]);//A_tilde = dummy*dummy_A;
					dummy_B.resize(sqr(this->summations[mu])*sqr(this->summations[mu+1]));
					fill_B_mu(mu, dummy_B);
					dmrg_memory.dummy.resize(B_tilde.size());
					std::memcpy(&(dmrg_memory.dummy[0]), &B_tilde[0], this->TSIZE*sqr(this->summations[0])*
								sqr(this->summations[mu])); // dummy = B_tilde
					B_tilde.resize(sqr(this->summations[0])*sqr(this->summations[mu+1]));
					Blas<T>::gemm('N', 'N', sqr(this->summations[0]),
								  sqr(this->summations[mu+1]),
								  sqr(this->summations[mu]), 1, &(dmrg_memory.dummy[0]),
								  &dummy_B[0], 0, &B_tilde[0]); //B_tilde = dummy*dummy_B;


				}
				dummy_A.resize(A_tilde.size());
				reshape_into(representation, d-2, 0, A_tilde, dummy_A); //dummy_A = reshape(A_tilde);
				dummy_B.resize(B_tilde.size());
				reshape_into(d-2, 0, B_tilde, dummy_B); //dummy_B = reshape(B_tilde)
				performSingleDMRGStep(representation, d-2, dummy_A, dummy_B, &dmrg_memory, eps);

				std::cout << std::sqrt(std::abs((dmrg_memory.norm_sqr_u+dmrg_memory.norm_sqr_a-2*dmrg_memory.norm_sqr_aou)/dmrg_memory.norm_sqr_a)) << " vs " << rel_l2dist(representation) << std::endl;

			}

			int getSubsequentMaxRankProduct(const TensorChainRepresentation<T> &representation) const {
				int result = this->summations[0]*this->summations[1]*representation.summations[0]*
										representation.summations[1];

				int buffer;

				for (int mu = 1; mu < this->d; mu++) {
					const int mu_plus_1 = (mu+1) % this->d;

					buffer = this->summations[mu]*this->summations[mu_plus_1]*representation.summations[mu]*
							   representation.summations[mu_plus_1];
					if (buffer > result) {
						result = buffer;
					}
				}
				return result;
			}

			T l2norm_post_sqr(const TensorChainRepresentation<T> &representation,
					      const std::vector<T> &post, std::vector<T> &buffer) const {
				int r = representation.summations[0]*this->summations[0];

				/* now, dummy_1 contains B_[d] */
				fill_A_mu(representation, this->d-1, buffer);

				T norm = 0.0;

				for (int n = 0; n < this->summations[this->d-1]; n++) {
					for (int k = 0; k < representation.summations[this->d-1]; k++) {
						norm += Blas<T>::dot(r, &buffer[n*representation.summations[this->d-1]+k],
											this->summations[this->d-1]*representation.summations[this->d-1],
											&post[(n*representation.summations[this->d-1]+k)*r], 1);
					}
				}
				return norm;
			}

			T l2norm_post(const TensorChainRepresentation<T> &representation,
								      const std::vector<T> &post, std::vector<T> &buffer) const {
				return std::sqrt(l2norm_post_sqr(representation, post, buffer));
			}

			T l2norm_post_sqr(const std::vector<T> &post, std::vector<T> &buffer) const {
				return l2norm_post_sqr(*this, post, buffer);
			}

			T l2norm_post(const std::vector<T> &post, std::vector<T> &buffer) const {
				return std::sqrt(l2norm_post_sqr(*this, post, buffer));
			}

			T l2norm_sqr(const TensorChainRepresentation<T> &representation) const {
				int r = this->summations[1]*representation.summations[1];

				int s = std::max(this->getMaxRank(2)*representation.getMaxRank(2), this->summations[0]*representation.summations[0]);

				int r_max_sqr = getSubsequentMaxRankProduct(representation);

				std::vector<T> dummy_1(r*s), dummy_2(r_max_sqr), dummy_3(s*r);

				fill_A_mu(representation, 1, dummy_1);

				for (int mu = 2; mu < this->d; mu++) {
					const int mu_plus_1 = (mu+1) % this->d;

					fill_A_mu(representation, mu, dummy_2);
					Blas<T>::gemm('N', 'N', r, this->summations[mu_plus_1]*representation.summations[mu_plus_1],
							      this->summations[mu]*representation.summations[mu], 1.0,  &dummy_1[0], &dummy_2[0],
							      0.0, &dummy_3[0]);
					std::memcpy(&dummy_1[0], &dummy_3[0],
							    this->TSIZE*r*this->summations[mu_plus_1]*representation.summations[mu_plus_1]);
				}
				return l2norm_pre_sqr(representation, dummy_1, dummy_2);
			}

			T l2norm_sqr() const {
				return l2norm_sqr(*this);
			}

			T l2norm() const {
				return std::sqrt(l2norm_sqr(*this));
			}

			T l2norm_pre_sqr(const TensorChainRepresentation<T> &representation,
					     const std::vector<T> &pre, std::vector<T> &buffer) const {
				int r = representation.summations[1]*this->summations[1];

				fill_A_mu(representation, 0, buffer);

				T norm = 0.0;

				for (int n = 0; n < this->summations[0]; n++) {
					for (int k = 0; k < representation.summations[0]; k++) {
						norm += Blas<T>::dot(r, &buffer[n*representation.summations[0]+k],
											 representation.summations[0]*this->summations[0],
								             &pre[(n*representation.summations[0]+k)*r], 1);
					}
				}
				return norm;
			}

			T l2norm_pre(const TensorChainRepresentation<T> &representation,
						 const std::vector<T> &pre, std::vector<T> &buffer) const {
				return std::sqrt(l2norm_pre_sqr(representation, pre, buffer));
			}

			T l2norm_pre_sqr(const std::vector<T> &pre, std::vector<T> &buffer) const {
				return l2norm_pre_sqr(*this, pre, buffer);
			}

			/*
			 * pre contains the matrix-product except the first one
			 */
			T l2norm_pre(const std::vector<T> &pre, std::vector<T> &buffer) const {
				return std::sqrt(l2norm_pre_sqr(*this, pre, buffer));
			}

			T l2dist(const TensorChainRepresentation<T> &representation) const {
				return std::sqrt(std::abs(l2dist_sqr(representation)));
			}

			T l2dist_sqr(const TensorChainRepresentation<T> &representation) const {
				return l2norm_sqr() -2*l2norm_sqr(representation)+ representation.l2norm_sqr();
			}

			T rel_l2dist_sqr(const TensorChainRepresentation<T> &representation) const {
				T norm = representation.l2norm_sqr();

				return std::abs((l2norm_sqr() -2*l2norm_sqr(representation)+ norm) / norm);
			}

			T rel_l2dist(const TensorChainRepresentation<T> &representation) const {
				return std::sqrt(rel_l2dist_sqr(representation));
			}

			/*
			 * Adjusts edge mu_plus_1 !
			 */
			void acaApproximation(const int mu, int newRank) {
				const int mu_plus_1 = (mu+1) % this->d,
						  mu_plus_2 = (mu+2) % this->d,
						  dim1 = this->componentDimensions[mu],
						  dim2 = this->componentDimensions[mu_plus_1],
						  dim_prod = dim1*dim2;

				std::vector<T> matrix(dim_prod*this->summations[mu]*this->summations[mu_plus_2]);

				std::vector<T> buffer(std::max(dim_prod, newRank*this->summations[mu_plus_2]*dim2));

				for (int n = 0; n < this->summations[mu_plus_1]; n++) {
					for (int k = 0; k < this->summations[mu]; k++) {
						for (int l = 0; l < this->summations[mu_plus_2]; l++) {
							Blas<T>::scal(dim_prod, 0.0, &buffer[0], 1); // reset the buffer
							Blas<T>::ger(dim1, dim2,
									     1.0, &(this->v[mu][dim1*(n*this->summations[mu]+k)]), 1, &(this->v[mu_plus_1][dim2*(l*this->summations[mu_plus_1]+n)]), 1,
									     &buffer[0], dim1);
							for (int h = 0; h < dim2; h++) {
								Blas<T>::axpy(dim1, 1.0,
										      &buffer[h*dim1], 1,
										      &matrix[(dim2*l+h)*this->summations[mu]*dim1 + dim1*k], 1);
							}
						}
					}
				}

				newRank = aca(matrix, dim1*this->summations[mu], dim2*this->summations[mu_plus_2],
					          newRank, this->v[mu], buffer);
				this->setSummation(mu_plus_1, newRank, false);
				for (int k = 0; k < newRank; k++) {
					for (int l = 0; l < this->summations[mu_plus_2]; l++) {
						Blas<T>::copy(dim2, &buffer[dim2*(k*this->summations[mu_plus_2]+l)], 1,
								      &(this->v[mu_plus_1][dim2*(l*newRank+k)]), 1);
					}
				}
			}
	};

	/*
	 * Returns n_mu n_mu_p1 x R_mu R_mu_p1 matrix
	 */
	template<typename T>
	void fill_a_mu_mu_plus_1(const TensorChainRepresentation<T> &representation, int mu,
							 std::vector<T> &destination) {
		int mu_plus_1 = (mu + 1) % representation.getD(),
			mu_plus_2 = (mu + 2) % representation.getD();

		int n_mu = representation.getComponentDimension(mu),
			n_mu_plus_1 = representation.getComponentDimension(mu_plus_1);

		int R_mu = representation.getSummation(mu),
		    R_mu_plus_1 = representation.getSummation(mu_plus_1),
		    R_mu_plus_2 = representation.getSummation(mu_plus_2);
		destination.resize(n_mu*n_mu_plus_1*R_mu*R_mu_plus_2);
		Blas<T>::scal(destination.size(), 0.0, &destination[0], 1); // reset data
		for (int i_ell_plus_1 = 0; i_ell_plus_1 < R_mu_plus_1; i_ell_plus_1++) {
			for (int i_ell = 0; i_ell < R_mu; i_ell++) {
				for (int i_ell_plus_2 = 0; i_ell_plus_2 < R_mu_plus_2; i_ell_plus_2++) {
					Blas<T>::ger(n_mu, n_mu_plus_1, 1, &(representation[mu][n_mu*(i_ell_plus_1*R_mu+i_ell)]),
							     1, &(representation[mu_plus_1][n_mu_plus_1*(i_ell_plus_2*R_mu_plus_1+i_ell_plus_1)]),
							     1, &destination[n_mu*n_mu_plus_1*(i_ell_plus_2*R_mu+i_ell)], n_mu);
				}

			}
		}

	}

	template<typename T>
	void reshape_u_ell_ell_plus_1(const TensorChainRepresentation<T> &representation, int mu,
			                      const std::vector<T> &source, std::vector<T> &destination) {
		int mu_plus_1 = (mu + 1) % representation.getD(),
			mu_plus_2 = (mu + 2) % representation.getD();

		int n_mu = representation.getComponentDimension(mu),
			n_mu_plus_1 = representation.getComponentDimension(mu_plus_1);

		int r_mu = representation.getSummation(mu),
		    r_mu_plus_2 = representation.getSummation(mu_plus_2);

		destination.resize(n_mu*r_mu*n_mu_plus_1*r_mu_plus_2);
		for (int z_ell = 0; z_ell < n_mu; z_ell++) {
			for (int z_ell_plus_1 = 0; z_ell_plus_1 < n_mu_plus_1; z_ell_plus_1++) {
				for (int j_ell = 0; j_ell < r_mu; j_ell++) {
					//use memcpy
					for (int j_ell_plus_2 = 0; j_ell_plus_2 < r_mu_plus_2; j_ell_plus_2++) {
						const int row_s = n_mu*z_ell_plus_1+z_ell,
							      column_s = j_ell_plus_2*r_mu+j_ell,
							      row_d = n_mu*j_ell+z_ell,
							      column_d = n_mu_plus_1*j_ell_plus_2+z_ell_plus_1;

						destination[column_d*n_mu*r_mu+row_d] = source[column_s*n_mu*n_mu_plus_1+row_s];
					}
				}
			}
		}
	}

	template<typename T>
	void apply_svd_result(TensorChainRepresentation<T> &representation, const std::vector<T> &S, std::vector<T> &U,
			              std::vector<T> &VT, const int rank, const int mu) {
		const int mu_plus_1 = (mu + 1) % representation.getD(),
				  mu_plus_2 = (mu + 2) % representation.getD();

		const int n_mu = representation.getComponentDimension(mu),
				  n_mu_plus_1 = representation.getComponentDimension(mu_plus_1);

		const int r_mu = representation.getSummation(mu),
				  r_mu_plus_2 = representation.getSummation(mu_plus_2);

		const int orig_new_rank = std::min(n_mu*r_mu, n_mu_plus_1*r_mu_plus_2);

		for (int l = 0; l < rank; l++) {
			Blas<T>::scal(n_mu_plus_1*r_mu_plus_2, S[l], &VT[l], orig_new_rank); // put the factors in the right term, not in the left one
		}

		representation.setSummation(mu_plus_1, rank, false);

		std::memcpy(&representation(mu,0), &U[0], n_mu*r_mu*rank*sizeof(T));

		for (int l = 0; l < rank; l++) { // here we simulate the storage order
			for (int j_mu_plus_2 = 0; j_mu_plus_2 < r_mu_plus_2; j_mu_plus_2++) {
				Blas<T>::copy(n_mu_plus_1, &VT[l+orig_new_rank*j_mu_plus_2*n_mu_plus_1], orig_new_rank,
						      &representation(mu_plus_1, n_mu_plus_1*(j_mu_plus_2*rank+l)), 1);
			}
		}
	}

	template<typename T>
	TensorChainRepresentation<T> createRandomTensorChain(std::vector<int> &summations, std::vector<int> &componentDimensions, T (*randomNumberGenerator)()) {
		int d = componentDimensions.size();

		std::vector< std::vector<T> > v_tc(d);

		for (int n = 0; n < d; n++) {
			v_tc[n].resize(summations[n]*summations[(n+1) % d]*componentDimensions[n]);
			for (unsigned int k = 0; k < v_tc[n].size(); k++) {
				v_tc[n][k] = randomNumberGenerator();
			}
		}

		TensorChainRepresentation<T> tensorChain(summations, componentDimensions, v_tc);

		return tensorChain;
	}

	template<typename T>
	TensorChainRepresentation<T> createRandomTensorChain(int tcRank, std::vector<int> &componentDimensions, T (*randomNumberGenerator)()) {
		std::vector<int> summations(componentDimensions.size(), tcRank);

		return createRandomTensorChain(summations, componentDimensions, randomNumberGenerator);
	}

	template<typename T>
	TensorChainRepresentation<T> createRandomTensorChain(int rank, int d, int componentDimension, T (*randomNumberGenerator)()) {
		std::vector<int> componentDimensions(d, componentDimension);

		return createRandomTensorChain(rank, componentDimensions, randomNumberGenerator);
	}



}


#endif /* __TENSORCHAINREPRESENTATION_HPP */
