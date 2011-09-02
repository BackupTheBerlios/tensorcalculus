/*
 * Copyright (C) 2010, 2011 Stefan Handschuh
 *
 * init(const char * filename) based on DKTS::readDataFrom(const IString& file)
 * which is Copyright (C) 2008 Mike Espig
 *
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

#ifndef __TENSORREPRESENTATION_HPP
#define __TENSORREPRESENTATION_HPP

#define ARGUMENT_CHECKS_ON

#ifndef SCAL_BORDER
#define SCAL_BORDER 10e-4
#endif

#include "BlasInterface.hpp"
#include "Vector/VectorOperators.hpp"
#include "Matrix/MatrixOperators.hpp"
#include "LapackInterface2.hpp"
#include "Utilities/Index.hpp"
#include "Utilities/Utilities.hpp"
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept> 
#include <cmath>
#include <time.h>
#include <omp.h>
#include <iostream>
#include <algorithm>



namespace TensorCalculus {
		
	int vector2index(const std::vector<int> &summations, const std::vector<int> &vector);
	
	int partialVector2index(const std::vector<int> &vector, const std::vector<int> &summations,
			                const std::vector< std::vector<int> > &incidenceMatrix, int mu);
	
	std::vector<int> index2vector(const std::vector<int> &summations, int index);
	
	std::vector<int> index2partialVector(const std::vector<int> &summations, int index,
			                             const std::vector<int> &incidenceRow);

	int partialIndex(const std::vector<int> &summations,
			         const std::vector< std::vector<int> > &incidenceMatrix, int index, int mu);
	
	int partialVector2index(const std::vector<int> &summations, const std::vector<int> &vector,
			                const std::vector<int> incidenceRow);

	std::vector<int> getPartialSummations(const std::vector<int> &summations,
			                              const std::vector<int> &incidenceRow);
	
	std::vector<int> getPartialSummations(const std::vector<int> &summations,
			                              const std::vector<int> &incidenceRow, int k);
	
	long getComponentProduct(const std::vector<int> &vector);
	
	std::vector<int> genereateGenericVector(int length, int offset = 0);

	void fillGenericVector(std::vector<int> &dest, int length, int offset = 0, bool resize = true);

	template<typename T> T componentScalarProduct(const std::vector< std::vector<T> > &v,
				                                      int direction, int i, const T * vector,
							                          const int dimension);

	template<typename T, typename U>
		std::vector<T>& tensorScalarProductComponentsDirection(const std::vector< std::vector<int> > &incidenceMatrix,
				                                               const std::vector<int> &summations,
				                                               const std::vector< std::vector<T> > &v,
				                                               const U &representation,
				                                               std::vector<T> &result,
				                                               /*std::vector<T> &factors,*/
				                                               const int direction,
				                                               const int componentDimension);
	template <typename T, typename U>
		std::vector< std::vector<T> >& tensorScalarProductComponents(
				const std::vector< std::vector<int> > &incidenceMatrix,
	            const std::vector<int> &summations,
	            const std::vector< std::vector<T> > &v,
	            const U &representation,
	            std::vector< std::vector<T> >& result,
	            const std::vector<int> componentDimensions);


	template <typename T> class TensorRepresentation {
		protected:
			std::vector< std::vector<T> > v;
			
			std::vector<int> summations;
			
			std::vector<int> maxSummations;

			std::vector<long> partialSummationProducts;

			std::vector<int> componentDimensions;

			int d;
				
			int TSIZE;

		private:
			std::vector< std::vector<T> > w;

			std::vector< std::vector<int> > incidenceMatrix; 
			
			std::vector< std::vector<int> > incidenceTable;
			
			std::vector< std::vector<T> > basis;

			int L;
			
			int summationsCount;

			int summationCountProduct;

		private:

		
		protected:

			/* NOT optimal and not finished
			void init(const TensorRepresentation<T> &representation,
					  std::vector< std::vector<int> > &incidenceMatrix) {
				TSIZE = sizeof(T);
				d = representation.getD();
#ifdef ARGUMENT_CHECKS_ON
				if (d != incidenceMatrix.size()) {
					throw std::invalid_argument("Incidence-matrix and ... .");
				}
				if (representation.getL() > 0) {
					throw std::invalid_argument("Not implemented for L > 0.");
				}

#endif

				for (int n = 0; n < d; n++) {
					for (int i = 0, j = incidenceMatrix[n].size(); i < j; i++) {
						if (indexOf(summations, incidenceMatrix[n][i]) == -1) {
							// use summations as a buffer
							summations.push_back(incidenceMatrix[n][i]);
						}
					}
				}
				summationsCount = summations.size();


				componentDimensions = std::vector<int>(representation.componentDimensions);

				Index index(representation.summations);

				long indexSize = index.size();

				for (int n = 0; n < summationsCount; n++) {
					summations[n] = indexSize;
				}
				partialSummationProducts.resize(d);
				v.resize(d);
				for (int n = 0; n < d; n++) {
					int nodeDegree = incidenceMatrix[n].size();

					partialSummationProducts[n] = std::pow(indexSize, nodeDegree);
					v[n].resize(partialSummationProducts[n]*componentDimensions[n]);
					std::vector<int> nodeSummations(incidenceMatrix[n].size(), indexSize);
					std::vector<int> nodePosition(incidenceMatrix[n].size());
					for (index.begin(); !index.end(); ++index) {

						const int pos = TensorCalculus::partialVector2index(representation.summations,
								                                            index,
								                                            incidenceMatrix[n]);
						const int pso = TensorCalculus::
						std::memcpy(&v[n][pso*componentDimensions[n]],
								    &representation[n][pos*componentDimensions[n]],
								    componentDimensions[n]*TSIZE);


					}
				}

			}

			*/

			void init(const std::vector<int> &summations, const std::vector< std::vector<T> > &v,
					  const std::vector<int> &componentDimensions,
					  const std::vector< std::vector<int> > &incidenceMatrix) {
				init(summations, v, componentDimensions, incidenceMatrix,
				     std::vector< std::vector<T> >(0));
			}
			
			void init(const std::vector<int> &summations, const std::vector< std::vector<T> > &v,
					  const std::vector<int> &componentDimensions,
					  const std::vector< std::vector<int> > &incidenceMatrix,
					  const std::vector< std::vector<T> > &w) {

				L = w.size();
 				d = v.size();
 				TSIZE = sizeof(T);
#ifdef ARGUMENT_CHECKS_ON
				if (static_cast<unsigned int>(d) != componentDimensions.size()) {
					throw std::invalid_argument("The number of components does not match the "
							                    "number of component-dimensions.");
				}
				if (incidenceMatrix.size() > 0 && incidenceMatrix.size() != static_cast<unsigned int>(d + L)) {
					throw std::invalid_argument("The number of components does not match the "
							                    "number of component-connections.");
				}
#endif
				summationsCount = summations.size();
#ifdef ARGUMENT_CHECKS_ON
				if (summationsCount == 0) {
					throw std::invalid_argument("There is no summation given so the tensor would "
							                    "be zero.");
				}
#endif				
				summationCountProduct = 1;
				
				for (int n = 0; n < summationsCount; n++) {
					summationCountProduct *= summations[n];
				}
#ifdef ARGUMENT_CHECKS_ON
				if (summationCountProduct == 0) {
					throw std::invalid_argument("There is no summation given so the tensor would "
							                    "be zero.");
				}
#endif					
				partialSummationProducts.resize(d+L);
				for (int n = 0; n < d + L; n++)	{
					int countProduct = 1;
					
					if (incidenceMatrix.size() == 0) {
						countProduct = summationCountProduct;
					} else { 
						const std::vector<int> &incidenceRow = incidenceMatrix[n];
						
						for (int k = 0, l = incidenceRow.size(); k < l; k++) {
#ifdef ARGUMENT_CHECKS_ON
							if (summationsCount <= incidenceRow[k]) {
								throw std::invalid_argument("The incidence matrix does not match "
										                    "the summation indices.");
							}
#endif
							countProduct *= summations[incidenceRow[k]];
						}
					}
#ifdef ARGUMENT_CHECKS_ON
					if (n < d) {
						if (v[n].size() != static_cast<unsigned int>(componentDimensions[n]*countProduct)) {
							throw std::invalid_argument("The component vector v at position n "
									                    "does not have the correct size.");
						}
					} else {
						if (w[n-d].size() != static_cast<unsigned int>(countProduct)) {
							throw std::invalid_argument("The component vector w at position n "
									                    "does not have the correct size.");
						}
					}
#endif	
					partialSummationProducts[n] = countProduct;
				}
			
				(*this).summations = summations;
				(*this).v = v;
				(*this).w = w;
				(*this).componentDimensions = componentDimensions;
				if (incidenceMatrix.size() == 0) {
					(*this).incidenceMatrix.resize(d+L);
					incidenceTable.resize(1);
					incidenceTable[0].resize(d+L);
					for (int n = 0; n < d+L; n++) {
						(*this).incidenceMatrix[n].resize(1);	
						incidenceTable[0][n] = n;
					}		
					
				} else {
					(*this).incidenceMatrix = incidenceMatrix;
					for (int n = 0; n < d+L; n++) {
						for (int k = 0, l = incidenceMatrix[n].size(); k < l; k++) {
							if (incidenceTable.size() < static_cast<unsigned int>(incidenceMatrix[n][k]+1)) {
								incidenceTable.resize(incidenceMatrix[n][k]+1);
							}
							
							if (!contains(incidenceTable[incidenceMatrix[n][k]], n)) {
								incidenceTable[incidenceMatrix[n][k]].push_back(n);
							}	
						}
					}
				}
				
			}
			
			void init(const int r, const std::vector< std::vector<T> > &v,
					  const std::vector<int> &componentDimensions,
					  const std::vector< std::vector<int> > &incidenceMatrix,
					  const std::vector< std::vector<T> > &w) {
				int size = v.size() + w.size();

				std::vector<int> summations(size, r);
				init(summations, v, componentDimensions, incidenceMatrix, w);
			}

			void init(int r, const std::vector<std::vector<T> > &v,
					  const std::vector<int> &componentDimensions,
					  const std::vector<int> &incidenceVectorV) {
				std::vector<int> summations(1);
				
				summations[0] = r;
				
				int d = v.size();
				
				if (incidenceVectorV.size() == 0) {
					init(summations, v, componentDimensions, std::vector< std::vector<int> >(0));
				} else {
					std::vector< std::vector<int> > incidenceMatrix(d);
				
					for (int n = 0; n < d; n++) {
						incidenceMatrix[n].resize(1);
						incidenceMatrix[n][0] = incidenceVectorV[n];
					}
					init(summations, v, componentDimensions, incidenceMatrix);	
				}
			}

			void init(const std::vector<int> &summations,
					  const std::vector<int> &componentDimensions,
					  const std::vector< std::vector<int> > &incidenceMatrix) {
				int d = componentDimensions.size();

				std::vector< std::vector<T> > v(d);

				// hier wird nicht die Konsitenz der Anzahl der componentDimensions
				// mit der adjazenzmatrix überprüft
				for (int n = 0; n < d; n++) {
					int size = 1;

					for (int k = 0, l = incidenceMatrix[n].size(); k < l; k++) {
						size *= summations[incidenceMatrix[n][k]];
					}
					v[n].resize(size*componentDimensions[n]);
				}
				init(summations, v, componentDimensions, incidenceMatrix);
			}

			
			
			void init(const char * filename) {
				/* based on DKTS::readDataFrom(const IString& file) */
				std::ifstream file(filename);

				long d = 0;

				char t[64];

				file.setf(std::ios::scientific, std::ios::floatfield);

				file >> t >> t >> t >> t;
				file >> d;

				std::vector<int> componentDimensions(d);

				for (int mu = 0; mu < d; mu++) {
					file >> t >> t;
					file >> componentDimensions[mu];
				}

				std::vector< std::vector<T> > v(d);

				file >> t;

				std::string r = t;

				if ("Rank" == r) { // old format --> CP
					int k;

					file >> t >> t;
					file >> k;

					for (int mu = 0; mu < d; mu++) {
					  	v[mu].resize(k*componentDimensions[mu]);
					}
					for (int j = 0; j < k; j++) {
					    for (int mu = 0; mu < d; mu++) {
					        for (int i = 0, n = componentDimensions[mu]; i < n; i++) {
					        	file >> v[mu][j*n+i];
							}
					    }
					}
					file.close();
					init(k, v, componentDimensions, std::vector<int>(0));
				} else { // new format
					file >> t;

					int z;

					std::vector< std::vector<int> > incidenceMatrix(d);

					std::vector<int> edges;

					for (int n = 0; n < d; n++) {
						file.getline(t, 64);
						std::stringstream g(t);
						while(!(g >> z).fail()) {
							incidenceMatrix[n].push_back(z);
							if (!contains(edges, z)) {
								edges.push_back(z);
							}
						}
						file >> t >> t;
					}

					int summationCount = edges.size();

					std::vector< std::vector<T> > v(d);

					std::vector<int> summations(summationCount);
					for (int n = 0; n < summationCount-1; n++) {
						file >> summations[n];
						file >> t >> t;
					}
					file >> summations[summationCount-1];

					int summationProduct;

					for (int mu = 0; mu < d; mu++) {
						summationProduct = 1;

						for (int n = 0, k = incidenceMatrix[mu].size(); n < k; n++) {
							summationProduct *= summations[incidenceMatrix[mu][n]];
						}
						v[mu].resize(componentDimensions[mu]*summationProduct);
					}
					for (int mu = 0; mu < d; mu++) {
						std::vector<int> directionSummations =
								TensorCalculus::getPartialSummations(summations,
										                             incidenceMatrix[mu]);

						z = componentDimensions[mu];
						for (Index index(directionSummations); !index.end(); ++index) {
							for (int i = 0; i < z; i++) {
								file >> v[mu][index*z+i];
							}
						}
					}
					file.close();
					init(summations, v, componentDimensions, incidenceMatrix);
				}
			}

		public:
			TensorRepresentation(const std::vector<int> &summations,
					             const std::vector<std::vector<T> > &v,
					             const std::vector<int> &componentDimensions,
					             const std::vector< std::vector<int> > &incidenceMatrix) {
				init(summations, v, componentDimensions, incidenceMatrix);
			}

			TensorRepresentation(const std::vector<int> &summations,
					             const std::vector<std::vector<T> > &v,
					             const std::vector<int> &componentDimensions,
					             const std::vector< std::vector<int> > &incidenceMatrix,
					             const std::vector<std::vector<T> > &w) {
				init(summations, v, componentDimensions, incidenceMatrix, w);
			}

			TensorRepresentation(const std::vector<int> &summations,
					             const std::vector<int> &componentDimensions,
					             const std::vector< std::vector<int> > &incidenceMatrix) {
				init(summations, componentDimensions, incidenceMatrix);
			}

			TensorRepresentation(int r, const std::vector<std::vector<T> > &v,
					             const std::vector<int> &componentDimensions,
					             const std::vector<int> &incidenceVector) {
				init(r, v, componentDimensions, incidenceVector);
			}

			TensorRepresentation(int r, const std::vector<std::vector<T> > &v,
					             const std::vector<int> &componentDimensions) {
				init(r, v, componentDimensions, std::vector<int>(0));
			}

			TensorRepresentation(int r, const std::vector< std::vector<T> > &v,
								 const std::vector<int> &componentDimensions,
								 const std::vector< std::vector<T> > &w) {
				init(r, v, componentDimensions, std::vector< std::vector<int> >(0), w);
			}

			TensorRepresentation(int r, const std::vector< std::vector<T> > &v,
					             int componentDimension) {
				int componentCount = v.size();

				std::vector<int> componentDimensions(componentCount);

				for (int n = 0; n < componentCount; n++) {
					componentDimensions[n] = componentDimension;
				}
				init(r, v, componentDimensions, std::vector<int>(0));
			}

			TensorRepresentation(const char * filename) {
				init(filename);
			}

			TensorRepresentation() {
				
			}
			
			virtual void writeToFile(const char * filename) {
				std::ofstream file(filename);

				file.setf(std::ios::scientific, std::ios::floatfield);
				file << "TensorNetwork in d = " << d << std::endl;
				for (int n = 0; n < d; n++) {
					file << "n[" << n << "] = " << componentDimensions[n] << std::endl;
				}
				for (int n = 0; n < d; n++) {
					file << "i[" << n << "] =";
					for (int k = 0, j = incidenceMatrix[n].size(); k < j; k++) {
						file << ' ' << incidenceMatrix[n][k];
					}
					file << std::endl;
				}
				for (int n = 0, j = summations.size(); n < j; n++) {
					file << "s[" << n << "] = " << summations[n] << std::endl;
				}
				file.precision(20);
				for (int n = 0; n < d; n++) {
					std::vector<int> nodeSummations =
							TensorCalculus::getPartialSummations(summations, incidenceMatrix[n]);

					for (Index index(nodeSummations); !index.end(); ++index) {
						for (int k = 0, j = componentDimensions[n]; k < j; k++) {
							file << v[n][index*j + k] << std::endl;
						}
					}
				}
				file.close();
			}

			/*
			 *
			 */
			const T evaluateAt(const std::vector<int> &index) const {
#ifdef ARGUMENT_CHECKS_ON
				//TODO implement
				/*
				 * we need to check:
				 *
				 *
				 */
#endif
				T result = 0;
				
				for (Index index2(summations); !index2.end(); ++index2) {
					T product = 1.0;
					
					for (int n = 0; n < d; n++) {
						if (product == 0) {
							break;
						}
						product *= v[n][partialVector2index(summations, index2,
								        incidenceMatrix[n])*componentDimensions[n]+index[n]];
						
					}	
					// TODO use w
					result += product;
				}
				return result;
			}

			const std::vector<T>& operator [] (int i) const {
				return v[i];
			}
			
			T& operator () (int i, int j) {
				return v[i][j];
			}

			const int getSummation(const int index) const {
#ifdef ARGUMENT_CHECKS_ON
				if (summationsCount <= index) {
					throw std::out_of_range("The requested summation index does not exist.");
				}
#endif		
				return summations[index]; 
			}
			
			const std::vector<int>& getIncidenceRow(const int index) const {
#ifdef ARGUMENT_CHECKS_ON
				if (d + L <= index) {
					throw std::out_of_range("The requested incidence row does not exist.");
				}
#endif	
				return incidenceMatrix[index];
			}
			
			/*
			 * Fills the coefficients matrix with the given vector x at position v.
			 * 
			 * v can be either in full representation (0,0,1,0,1) or in sparse representation (1,1)
			 */ 
			void fill(int mu, const std::vector<int> &v, const std::vector<T> &x) {
#ifdef ARGUMENT_CHECKS_ON
				if (d <= mu) {
					throw std::out_of_range("The requested myu index does not exist.");
				} else if (x.size() != componentDimensions[mu]) {
					throw std::invalid_argument("The the dimension of x does not match the "
							                    "component dimension.");
				}
				if (v.size() == summationsCount) {  // in case of a complete v is given
					for (int n = 0; n < summationsCount; n++) {
						if (summations[n] <= v[n]) {
							throw std::out_of_range("The position exceeds the summation.");
						}
					}
				} else {  // in case of v is given in sparse representation
					if (v.size() != incidenceMatrix[mu].size()) {
						throw std::invalid_argument("Invalid v index position descriptor.");
					} else {
						for (int n = 0, i = v.size(); n < i; n++) {
							if (summations[incidenceMatrix[mu][n]] <= v[n]) {
								throw std::out_of_range("The position exceeds the summation.");
							}
						}
					}
				}
#endif			
				std::memcpy(this.v[mu][partialVector2index(summations, v, incidenceMatrix, mu)
				                       *componentDimensions[mu]], &x[0], x.size()*TSIZE);
			}
			
			int getNodeDegree(int mu) const {
#ifdef ARGUMENT_CHECKS_ON
				if (mu < 0 || mu >= d) {
					throw std::invalid_argument("Mu is out of range");
				}
#endif
				return incidenceMatrix[mu].size();
			}

			/*
			 * Returns the nodes that are connected to edge k.
			 */
			std::vector<int> getAffectedNodes(int k) const {
				std::vector<int> result;

				for (int n = 0; n < d + L; n++) {
					if (indexOf(incidenceMatrix[n], k) > -1) {
						result.push_back(n);
					}
				}
				return result;
			}

			/*
			 * Returns the nodes that are not connected to edge k.
			 */
			std::vector<int> getNotAffectedNodes(int k) const {
				std::vector<int> result;

				for (int n = 0; n < d + L; n++) {
					if (indexOf(incidenceMatrix[n], k) == -1) {
						result.push_back(n);
					}
				}
				return result;
			}

			std::vector<int> getNotAffectedNodes(int edge, const std::vector<int> &fixedNodes) const {
				std::vector<int> result(fixedNodes);

				for (int mu = 0; mu < d + L; mu++) {
					if (indexOf(incidenceMatrix[mu], edge) == -1 && indexOf(result, mu) == -1) {
						result.push_back(mu);
					}
				}
				return result;
			}

			/*
			 * Returns all nodes that are connected to the edge except the given nodes.
			 */
			std::vector<int> getNodesExcept(int edge, const std::vector<int> &fixedNodes) const {
				std::vector<int> result;

				for (int mu = 0; mu < d + L; mu++) {
					if (indexOf(fixedNodes, mu) == -1 && indexOf(incidenceMatrix[mu], edge) > -1
							&& indexOf(result, mu) == -1) {
						result.push_back(mu);
					}
				}
				return result;
			}


			/*
			 * Returns the uncommon summations of node k1 and k2.
			 */
			std::vector<int> getUncommonSummations(int k1, int k2) const {
				std::vector<int> result;

				result.reserve(incidenceMatrix[k1].size() + incidenceMatrix[k2].size());
				for (int k = 0; k < summationsCount; k++) {
					if (indexOf(incidenceMatrix[k1], k) > -1 xor
						indexOf(incidenceMatrix[k2], k) > -1) {
						result.push_back(summations[k]);
					}
				}
				return result;
			}

			/*
			 * Returns true, if node k1 and node k2 have an edge in common, false otherwise
			 */
			bool haveCommonEdge(int k1, int k2) const {
				return getEdgeBetween(k1, k2) > -1;
			}

			int getMaxRank(int offset = 0) const {
				int result = summations[offset];

				for (int n = offset+1; n < summationsCount; n++) {
					if (summations[n] > result) {
						result = summations[n];
					}
				}
				return result;
			}

			int getMaxComponentDimension() const {
				int result = componentDimensions[0];

				for (int n = 1; n < d; n++) {
					if (componentDimensions[n] > result) {
						result = componentDimensions[n];
					}
				}
				return result;
			}

			long getMaxRankProduct() const {
				long result = partialSummationProducts[0];

				for (int n = 1; n < d; n++) {
					if (partialSummationProducts[n] > result) {
						result = partialSummationProducts[n];
					}
				}
				return result;
			}

			/*
			 * Returns the incidence row of node k except the e-numbered index.
			 *
			 * It is assumed, that e is in k-th incidencerow
			 */
			std::vector<int> getIncidenceRowExceptEdge(int k, int e) const {
				int l = incidenceMatrix[k].size();

				std::vector<int> result(l -1);

				for (int n = 0, i = 0; n < l; n++) {
					if (incidenceMatrix[k][n]!= e) {
						result[i++] = incidenceMatrix[k][n];
					}
				}
				return result;
			}

			/*
			 * Returns the summations of node k except the e-th summation.
			 *
			 * It is assumed, that e is in k-th incidencerow
			 */
			std::vector<int> getSummationsExceptEdge(const int k, const int e) const {
				int l = incidenceMatrix[k].size();

				std::vector<int> result(l -1);

				for (int n = 0, i = 0; n < l; n++) {
					if (incidenceMatrix[k][n]!= e) {
						result[i++] = summations[incidenceMatrix[k][n]];
					}
				}
				return result;
			}

			std::vector<int> getEdgesExceptEdge(const std::vector<int> &nodes,
					                                 const int e) const {
				std::vector<int> result;

				for (int n = 0, i = nodes.size(); n < i; n++) {
					for (int k = 0, l = incidenceMatrix[nodes[n]].size(); k < l; k++) {
						if (incidenceMatrix[nodes[n]][k] != e && !contains(result, incidenceMatrix[nodes[n]][k])) {
							result.push_back(incidenceMatrix[nodes[n]][k]);
						}
					}
				}
				return result;
			}

			std::vector<int> getEdges(const std::vector<int> &nodes) const {
				std::vector<int> result;

				for (int n = 0, i = nodes.size(); n < i; n++) {
					for (int k = 0, l = incidenceMatrix[nodes[n]].size(); k < l; k++) {
						if (!contains(result, incidenceMatrix[nodes[n]][k])) {
							result.push_back(incidenceMatrix[nodes[n]][k]);
						}
					}
				}
				return result;
			}

			std::vector<int> getSummationsOfNode(int mu) const {
				int size = incidenceMatrix[mu].size();

				std::vector<int> result(size);
				for (int n = 0; n < size; n++) {
					result[n] = summations[incidenceMatrix[mu][n]];
				}
				return result;
			}

			std::vector<int> getSummationsExceptNodes(const int k1, const int k2) const {
				std::vector<int> result;

				result.reserve(summations.size() - incidenceMatrix[k1].size() -
						       incidenceMatrix[k2].size() + 1);
				for (int k = 0; k < summationsCount; k++) {
					if (indexOf(incidenceMatrix[k1], k) == -1 &&
						indexOf(incidenceMatrix[k2], k) == -1) {
						result.push_back(summations[k]);
					}
				}
				return result;
			}

			/*
			 * Returns the unique items of getIncidenceRow(k1) and getIncidenceRow(k2)
			 */
			std::vector<int> getIncidenceRow(const int k1, const int k2) const {
				std::vector<int> result(incidenceMatrix[k1]);

				int size = incidenceMatrix[k2].size();

				result.reserve(size + incidenceMatrix[k1].size());
				for (int k = 0; k < size; k++) {
					if (indexOf(incidenceMatrix[k1], incidenceMatrix[k2][k]) == -1) {
						result.push_back(incidenceMatrix[k2][k]);
					}
				}
				return result;
			}

			/*
			 * Returns the edge between node k1 and node k2, -1 if there is no edge
			 */
			int getEdgeBetween(const int k1, const int k2) const {
				for (int n = 0, i = incidenceMatrix[k1].size(); n < i; n++) {
					if (indexOf(incidenceMatrix[k2], incidenceMatrix[k1][n]) > -1) {
						return incidenceMatrix[k1][n];
					}
				}
				return -1;
			}
			
			TensorRepresentation& hadamardAdd(const TensorRepresentation& representation) {
#ifdef ARGUMENT_CHECKS_ON
				if (summations != representation.summations) {
					throw std::invalid_argument("Summations do not match.");
				} 
				if (componentDimensions != representation.componentDimensions) {
					throw std::invalid_argument("Component dimensions do not match.");
				}
				if (incidenceMatrix != representation.incidenceMatrix) {
					throw std::invalid_argument("Incidence matrix does not match.");
				}
#endif
				for (int n = 0; n < d; n++) {
					Blas<T>::axpy(v[n].size(), 1.0, &representation.v[n][0], 1, &v[n][0], 1);
				}
				for (int n = 0; n < L; n++) {
					Blas<T>::axpy(w[n].size(), 1.0, &representation.w[n][0], 1, &w[n][0], 1);
				}
				return *this;
			}
			
			TensorRepresentation& setAsHadamardSumOf(const TensorRepresentation& representation1,
					                                 const TensorRepresentation& representation2) {
				// maybe add some checks
				(*this).summations = representation1.summations;
				(*this).v = representation1.v;
				(*this).w = representation1.w;
				(*this).componentDimensions = representation1.componentDimensions;
				(*this).incidenceMatrix = representation1.incidenceMatrix;
				return hadamardAdd(representation2);
			}
			
			virtual TensorRepresentation& add(const TensorRepresentation& representation) {
				throw std::invalid_argument("No implementation of the summation found.");
			};
			
			TensorRepresentation& scale(const T& factor) {
				T ballancedFactor = pow(factor, (T) 1 / (T) (d+L));
		
				if (ballancedFactor < SCAL_BORDER) {
					Blas<T>::scal(v[0].size(), factor, &v[0][0], 1); // since d >= 1
				} else {
					for (int i = 0, n = d; i < n; i++) {
						Blas<T>::scal(v[i].size(), ballancedFactor, &v[i][0], 1); 
					}
					for (int i = 0, n = L; i < n; i++) {
						Blas<T>::scal(w[i].size(), ballancedFactor, &w[i][0], 1); 
					}	
				}
				return *this;
			}
			
			const T normalize() {
				T norm = l2norm(*this);

				scale(1.0/norm);
				return norm;
			}

			const int getL() const { 
				return L;
			}
			
			const int getD() const {
				return d;
			}
			
			template <typename U>
			void checkCompatibility(const U &representation) const {
				if (d != representation.getD()) {
					throw std::invalid_argument("The tensor component count is not equal.");
				}
				if (componentDimensions != representation.getComponentDimensions()) {
					throw std::invalid_argument("The component dimensions are not equal.");
				}
			}

			const std::vector<int>& getSummations() const {
				return summations;	
			}
			
			const int getSummationsCount() const {
				return summationsCount;
			}				
			
			const std::vector<int>& getComponentDimensions() const {
				return componentDimensions;
			}
			
			int getComponentDimension(int k) const {
#ifdef ARGUMENT_CHECKS_ON
				if (k < 0 || static_cast<unsigned int>(k) >= componentDimensions.size()) {
					throw std::invalid_argument("The given k is out of range.");
				}
#endif
				return componentDimensions[k];
			}

			const std::vector< std::vector<T> >& getV() const {
				return v;	
			}
			
			void setV(const std::vector< std::vector<T> > &v) {
				// TODO add checks since this is a critical operation

				this->v = std::vector< std::vector<T> >(v);
			}

			const std::vector< std::vector<T> >& getW() const {
				return w;	
			}
			
			const std::vector< std::vector<int> >& getIncidenceMatrix() const {
				return incidenceMatrix;
			}
			
			T scalarProduct(const TensorRepresentation<T> &representation) const {
#ifdef ARGUMENT_CHECKS_ON
				(*this).checkCompatibility(representation);
#endif
				T result = 0;
				
				for (int n = 0, i = d; n < i; n++) {
					result += Blas<T>::dot(v[n].size(), &v[n][0], 1, &representation.v[n][0], 1);
				}
				for (int n = 0, i = L; n < i; n++) {
					result += Blas<T>::dot(w[n].size(), &w[n][0], 1, &representation.w[n][0], 1);
				}
				return result;
			}
			
			const T componentScalarProduct(int direction, int i, const T * vector,
					                       const int dimension) const {
				return TensorCalculus::componentScalarProduct(this->v, direction, i, vector, dimension);
			}

			/*
			 * Computes the inner products of this representation
			 * and the given on in the specified direction and stores
			 * the results in the result variable.
			 */
			template <typename U>
			std::vector<T>& tensorScalarProductComponentsDirection(const U &representation,
					                                               std::vector<T> &result,
					                                               /*std::vector<T> &factors,*/
					                                               const int direction) const {
				// we assume that everything is checked
				return TensorCalculus::tensorScalarProductComponentsDirection(incidenceMatrix,
						summations, v, representation, result, direction, componentDimensions[direction]);
			} 
			
			/*
			 * Computes all single inner products of this representation
			 * and the given one and stores them in the result variable.
			 */
			template <typename U>
			std::vector< std::vector<T> >& tensorScalarProductComponents(const U &representation,
					                               std::vector< std::vector<T> >& result) const {
				return TensorCalculus::tensorScalarProductComponents(incidenceMatrix,
						summations, v, representation, result, componentDimensions);
			}
			
			template <typename U>
			T tensorScalarProduct(const U &representation) const {
#ifdef ARGUMENT_CHECKS_ON
				(*this).checkCompatibility(representation);
#endif
				std::vector< std::vector<T> > scalarComponents;
				
				tensorScalarProductComponents(representation, scalarComponents);
				
				return tensorScalarProduct(scalarComponents, representation);
			}
			
			T tensorScalarProduct(const std::vector< std::vector<T> > &scalarComponents,
					              const TensorRepresentation<T> &representation) const {
				T product, result = 0;
				
				std::vector< std::vector<int> > position(d);
				
				std::vector< std::vector<int> > partialSummations(d);
				
				std::vector< std::vector<int> > position_index2(d);
				
				for (int mu = 0; mu < d; mu++) {
					partialSummations[mu].resize(2);
					partialSummations[mu][0] = 1;
					partialSummations[mu][1] = 1;
							
					for (int m = 0, l = incidenceMatrix[mu].size(); m < l; m++) {
						partialSummations[mu][0] *= summations[incidenceMatrix[mu][m]];
					}
					for (int m = 0, l = representation.incidenceMatrix[mu].size(); m < l; m++) {
						partialSummations[mu][1] *=
								representation.summations[representation.incidenceMatrix[mu][m]];
					}
					position[mu].resize(2);
					position_index2[mu].resize(representation.summationCountProduct);
					for (Index index2(representation.summations); !index2.end(); ++index2) {
						position_index2[mu][index2] =
								partialVector2index(representation.summations,
										            index2, representation.incidenceMatrix, mu);
					}
				}
				
				std::vector< std::vector<T> > summands(summationCountProduct*
						                 representation.summationCountProduct);

				for (Index index1(summations); !index1.end(); ++index1) {
					for (int mu = 0; mu < d; mu++) {	
						position[mu][0] = partialVector2index(summations, index1, incidenceMatrix, mu);
					}	
					for (int index2 = 0; index2 < representation.summationCountProduct; ++index2) {
						product = 1.0;
						for (int mu = 0; mu < d; mu++) {								
							position[mu][1] = position_index2[mu][index2];

							const T factor = scalarComponents[mu][vector2index(
									      partialSummations[mu], position[mu])];

							summands[summationCountProduct*index2 + index1].push_back(factor);
							product *= factor;
						}
						// TODO use w
						result += product;
					}	
				}


				/*
				T result2 = sum_prod(summands);// *norm;

				for (Index index1(summations); !index1.end(); ++index1) {
					for (Index index2(representation.getSummations()); !index2.end(); ++index2) {
						product = 1.0;
						for (int mu = 0; mu < d; mu++) {
							const int dimension = this->getComponentDimension(mu);

							product *= Blas<T>::dot(dimension,
								&((*this).v[mu][dimension*partialVector2index(summations, index1,
																		  incidenceMatrix, mu)]),
								1, &(representation.getV()[mu][dimension*partialVector2index(
									representation.getSummations(), index2,
									representation.getIncidenceMatrix(), mu)]), 1);
						}
						result2 += product;
					}
				}

				std::cout << "diff=" << result2 << " - " << result << " = " << result2 - result << std::endl;
				*/

				return result;
			}


			T partiall2norm_sqr(const std::vector< std::vector<T> > &scalarComponents, int edge,
					        int k1, int k2) const {
				int _k1, _k2;

				/* rearange such that _k1 <= _k2 */
				if (k1 > k2) {
					_k1 = k2;
					_k2 = k1;
				} else {
					_k1 = k1;
					_k2 = k2;
				}

				T result = 0;

				std::vector<int> fixedIndices = getIncidenceRow(k1, k2);

				std::vector<int> affectedSummations =
						TensorCalculus::getPartialSummations(summations, fixedIndices);

				std::vector< std::vector<int> > position(d);

				std::vector< std::vector<int> > partialSummations(d);

				for (int mu = 0; mu < d; mu++) {  // don't care about the two entries that are not needed
					partialSummations[mu].resize(2);
					partialSummations[mu][0] = componentProduct(
							TensorCalculus::getPartialSummations(summations, incidenceMatrix[mu]));
					partialSummations[mu][1] = partialSummations[mu][0];
					position[mu].resize(2);
				}

				for (Index index(affectedSummations); !index.end(); ++index) {
					Index index2(summations, fixedIndices, index);

					for (Index index1(summations, fixedIndices, index); !index1.end(); ++index1) {
						for (int mu = 0; mu < d; mu++) { // don't care about the two entries that are not needed
							position[mu][0] = partialVector2index(summations, index1,
									                              incidenceMatrix, mu);
						}
						for (index2.begin(); !index2.end(); ++index2) {
							T product = 1;

							for (int mu = 0; mu < _k1; mu++) {
								position[mu][1] = partialVector2index(summations, index2,
										                              incidenceMatrix, mu);
								product *= scalarComponents[mu][vector2index(partialSummations[mu],
										                                     position[mu])];
							}
							for (int mu = _k1+1; mu < _k2; mu++) {
								position[mu][1] = partialVector2index(summations, index2,
										                              incidenceMatrix, mu);
								product *= scalarComponents[mu][vector2index(partialSummations[mu],
										                                     position[mu])];
							}
							for (int mu = _k2+1; mu < d; mu++) {
								position[mu][1] = partialVector2index(summations, index2,
										                              incidenceMatrix, mu);
								product *= scalarComponents[mu][vector2index(partialSummations[mu],
										                                     position[mu])];
							}
							//std::cout << product << std::endl;
							result += product;
						}
					}
				}
				return result;
			}

			void toDot() const {
				std::ofstream fout;
				
				fout.open("representation.dot");
    			fout << "graph G {" << std::endl << "    ratio=1;" << std::endl;
    			// Write the number of numbers
    			for (int n = 0; n < summationsCount; n++) {
    				bool found = false;
    				
    				if (incidenceMatrix.size() > 0) {
    					std::vector<int> nodes;
    					
    					for (int i = 0; i < d+L; i++) {
    						
	    					for (int k = 0, l = incidenceMatrix[i].size(); k < l; k++) {
	    						if (incidenceMatrix[i][k] == n) {
	    							// v or w?
	    							nodes.reserve(nodes.size()+1);
	    							nodes.push_back(i < d ? i : -i);
	    						}
	    					}
	    				}	
	    				
	    				int node_size = nodes.size();

	    				switch (node_size) {
	    					case 0: {
	    						break;
	    					}
	    					case 1: {
	    						fout << "    " << (nodes[0] < 0 ? "w" : "v") << "_"
	    							 << std::abs(nodes[0]) <<";" << std::endl;
	    						break;
	    					}	
	    					case 2: {
	    						fout << "    " << (nodes[0] < 0 ? "w" : "v") << "_"
                                     << std::abs(nodes[0]) <<" -- " << (nodes[1] < 0 ? "w" : "v")
                                     << "_" << std::abs(nodes[1]) << ";" << std::endl;
	    						break;	
	    					}
	    					default: {
	    						fout << "    x_" << n << "[shape=point]";
	    						for (int i = 0; i < node_size; i++) {
	    							fout <<" x_" << n << "-- " << (nodes[i] < 0 ? "w" : "v")
	    								 << "_" << std::abs(nodes[i]) << std::endl;
	    						}	
	    						break;
	    					}
	    				}
    				} else {
    					fout << "    v_" << n; // irgendwie x_0 [shape=point] noch einbringen
	    				for (int i = 0; i < summationsCount; i++) {
	    					fout <<" -- x_" << i;
	    				}	
	    				fout << ";" << std::endl;
    				}
    				
    			}
    			fout << "}" << std::endl;
    			fout.close(); 
    			system("dot -Tps representation.dot -o representation.ps");
			}
			

			/**
			 * Prints matlab code to use the tensor representation there.
			 */
			void toMat(const char char_v = 'v', const char char_w = 'w') {
				for (int n = 0; n < L; n++) {
					std::vector<int> summations = getPartialSummations(n+d);

					for (Index index(summations); !index.end(); ++index) {
						std::cout << char_w << '_' << n << '_' << index.getPosition()
								  << " = " << w[n][index] << ';' << std::endl;
					}
				}
				for (int n = 0; n < d; n++) {
					std::vector<int> summations = getPartialSummations(n);

					const int componentDimension = componentDimensions[n];

					for (Index index(summations); !index.end(); ++index) {
						std::cout << char_v << '_' << n << '_' <<  index.getPosition()
								  << " = transpose([";
						for (int i = 0; i < componentDimension-1; i++) {
							std::cout << v[n][componentDimension*index+i] << ", ";
						}

						std::cout << v[n][componentDimension*index+componentDimension-1]
						          << "]);" << std::endl;
					}
				}

				std::cout << char_v << char_w << "_ten = ";

				for (Index index(summations); !index.end(); ++index) {
					for (int n = 0; n < L; n++) {
						std::cout << char_w << '_' << n << '_'
								  << partialVector2index(summations, index, incidenceMatrix[n+d])
								  << '*';

					}
					for (int n = d-1; n > 0; n--) {
						std::cout << "kron(" << char_v << '_' << n << '_'
								  << partialVector2index(summations, index, incidenceMatrix[n])
								  << ", ";
					}
					std::cout << char_v << '_' << 0 << '_'
							  << partialVector2index(summations, index, incidenceMatrix[0]);
					for (int n = 0; n < d-1; n++) {
						std::cout << ')';
					}
					std::cout << " + ";
				}
				std::cout << "0;" << std::endl;
				// "0" is needed because otherwise, we would end up with a "+" as the last sign
			}

			/*
			 * The following assumptions are made:
			 *  - index1 \in \times_{l \in I(k1) \cup I(k2)} N_{\leq r(l)}
			 *
			 *
			 */
			T tensorScalarProduct(const std::vector< std::vector<T> > &tensorScalarComponents,
					              const int k1, const int k2, const std::vector<int> &index1,
					              const std::vector<int> &index2,
					              const TensorRepresentation<T> &representation) const {
				int _k1, _k2;
				
				if (k1 > k2) {
					_k1 = k2;
					_k2 = k1;
				} else {
					_k1 = k1;
					_k2 = k2;
				}
				
				T result = 0;
				
				std::vector<int> incidenceRow1 = k1 == k2 ?
						representation.incidenceMatrix[k1] :
						representation.getIncidenceRow(k1, k2);
				
				Index countIndex2(representation.summations, incidenceRow1, index2);
				
				std::vector< std::vector<long> > offsetCountIndex2(d);
				
				for (int mu = 0; mu < d; mu++) {
					// don't care about the two entries that are not needed
					offsetCountIndex2[mu].resize(countIndex2.size());
				}

				for (countIndex2.begin(); !countIndex2.end(); ++countIndex2) {
					for (int mu = 0; mu < d; mu++) {
						// don't care about the two entries that are not needed
						offsetCountIndex2[mu][countIndex2] = partialSummationProducts[mu]*
								partialVector2index(representation.summations, countIndex2,
										            representation.incidenceMatrix[mu]);
					}
				}
				std::vector<int> incidenceRow2 = &representation != this && (k1 >= d || k2 >= d) ?
						genereateGenericVector(summations.size()) :
						(k1 == k2 ? incidenceMatrix[k1] : getIncidenceRow(k1, k2));
				//std::vector<T> summands;
				for (Index countIndex1(summations, incidenceRow2 , index1);
					 !countIndex1.end(); ++countIndex1) {
					std::vector<long> offset(d);
					
					for (int mu = 0; mu < d; mu++) {
						// don't care about the two entries that are not needed
						offset[mu] = partialVector2index(summations, countIndex1,
								                         incidenceMatrix[mu]);
					}

					T factor = 1.0;

					for (int mu = 0; mu < L; mu++) {
						if (&representation != this || (mu + d != k1 && mu + d != k2)) {
							factor *= w[mu][partialVector2index(summations, countIndex1,
									                            incidenceMatrix[mu+d])];
						}
					}
					for (countIndex2.begin(); !countIndex2.end(); ++countIndex2) {
						T product = factor;
						
						for (int mu = 0; mu < d; mu++) {
							if (mu != k1 && mu != k2) {
								// no tripple "for" possible because k1 >= d in some cases
								product *= tensorScalarComponents[mu][offset[mu] +
								                        offsetCountIndex2[mu][countIndex2]];
							}
						}

						for (int mu = 0; mu < representation.L; mu++) {
							if (mu + d != k1 && mu +d != k2) {
								product *= representation.w[mu][partialVector2index(
										representation.summations, countIndex2,
										representation.incidenceMatrix[mu+d])];
							}
						}
						//summands.push_back(product);
						result += product;
					}
				}
				//T ONE = 1.0;
				//T test = Blas<T>::dot(summands.size(), &summands[0], 1, &ONE, 0);
				return result;
			}
			
			T tensorScalarProduct(const std::vector< std::vector<T> > &tensorScalarComponents,
					              const int k, const std::vector<int> &index1,
					              const std::vector<int> &index2,
					              const TensorRepresentation<T> &representation) const {
				return tensorScalarProduct(tensorScalarComponents, k, k, index1, index2,
						                   representation);
			}
			
			void ballance(std::vector<int> & affectedNodes) {
				using namespace VectorOperators; // for *=
				for (int n = 0; n < summationsCount; n++) {
					if (summations[n] > 1) {
						throw std::invalid_argument("Only implemented for rank one tensors.");
					}
				}

				T product = 1.0;

				int i = affectedNodes.size();

				std::vector<T> norms(i);

				for (int n = 0; n < i; n++) {
					norms[n] = affectedNodes[n] < d ? l2_norm(v[affectedNodes[n]]) : std::abs(w[affectedNodes[n]-d][0]);
					product *= norms[n];
				}

				product = std::pow(product, 1.0/i);
				for (int n = 0; n < i; n++) {
					if (affectedNodes[n] < d) {
						v[affectedNodes[n]] *= (product / norms[n]);
					} else {
						w[affectedNodes[n]-d] *= (product / norms[n]);
					}
				}
			}

			const std::vector< std::vector<T> > & getBasis() const {
				return basis;
			}

			/**
			 * Computes the basis vectors and stores the coefficients in v.
			 */
			void computeOrthonormalBasis(T eps = 0.0) {
				if (basis.size() > 0) {
					return; // we have already computed the basis and set up the coefficients
				}


				int maxRankProduct = getMaxRankProduct();

				int maxDimension = getMaxComponentDimension();

				int s_count = std::min(maxDimension, maxRankProduct);

				double lwork;

				Lapack<T>::gesvd('S', 'S', maxDimension, maxRankProduct, 0, maxDimension, 0, 0,
						         maxDimension, 0, s_count, &lwork, -1);

				std::vector<T> workspace((int) lwork);

				std::vector<T> S(s_count); // storage can be optimized
				std::vector<T> VT(s_count*maxRankProduct); // storage can be optimized

				basis.resize(d);

				for (int n = 0; n < d; n++) {
					s_count = std::min(componentDimensions[n], (int) partialSummationProducts[n]);

					basis[n].resize(componentDimensions[n]*s_count);

					Lapack<T>::gesvd('S', 'S', componentDimensions[n], partialSummationProducts[n],
							         &v[n][0], componentDimensions[n], &S[0], &basis[n][0],
							         componentDimensions[n], &VT[0], s_count, &workspace[0],
							         (int) lwork);
					using namespace VectorOperators; std::cout << S << std::endl;
					/* determine the rank */
					T sum = S[s_count-1];

					int eff_s_count = s_count;

					for (int k = s_count - 2; k > -1 && sum < eps; k--) {
						sum += S[k];
						eff_s_count--;
					}
					v[n].resize(s_count*partialSummationProducts[n]);
					// put singular values into vt in order to take U as orthonormal basis
					for (int k = 0; k < eff_s_count; k++) {
						Blas<T>::scal(partialSummationProducts[n], S[k], &VT[k], s_count);
					}
					std::memcpy(&v[n][0], &VT[0], TSIZE*v[n].size());
					componentDimensions[n] = s_count;
				}
			}

			void includeOrthonormalBasis() {
				if (basis.size() == 0) {
					return; // no basis present
				}

				int maxSize = partialSummationProducts[0] * basis[0].size() / componentDimensions[0];

				for (int n = 1; n < d; n++) {
					const int b = partialSummationProducts[n] * basis[n].size() / componentDimensions[n];

					if (b > maxSize) {
						maxSize = b;
					}
				}
				std::vector<T> buffer(maxSize);
				for (int n = 0; n < d; n++) {
					const int realComponentDimension = basis[n].size() / componentDimensions[n];

					Blas<T>::gemm('N', 'N', realComponentDimension,
							      (int) partialSummationProducts[n], componentDimensions[n],
							      1.0, &basis[n][0], &v[n][0], 0.0, &buffer[0]);
					componentDimensions[n] = realComponentDimension;
					v[n].resize(realComponentDimension*partialSummationProducts[n]);
					std::memcpy(&v[n][0], &buffer[0],
							    TSIZE*partialSummationProducts[n]*realComponentDimension);
				}
				basis.resize(0);
			}

			/**
			 * Yes, we want to copy the basis parameter
			 */
			void setOrthonormalBasis(const std::vector< std::vector<T> > basis) {
				// to be sure not to have already one, basistransormation easier
				includeOrthonormalBasis();
				this->basis = basis;
			}

			int performALS(const TensorRepresentation &representation, T eps = 10e-6, TensorRepresentation *subtrahend = 0, std::vector<int> * changeNodes = 0) {
				std::vector< std::vector<T> > scalarComponentsUU, scalarComponentsBU,
				                              scalarComponentsAU;
				#pragma omp parallel sections
				{
					#pragma omp section
					{
						tensorScalarProductComponents((*this), scalarComponentsUU);
					}
					#pragma omp section
					{
						representation.tensorScalarProductComponents((*this), scalarComponentsBU);
					}
					#pragma omp section
					{
						if (subtrahend != 0) {
							(*subtrahend).tensorScalarProductComponents((*this), scalarComponentsAU);
						}
					}
				}

				T norm = l2norm(representation);

				T norm_this = 0.0, norm_this_0 = l2norm(*this);

				int stepCount = 0;

				int k; // node

				std::vector<int> changeV, changeW;

				if (changeNodes == 0) {
					fillGenericVector(changeV, d);
					fillGenericVector(changeW, L);
				} else {
					for (int n = 0, l = changeNodes->size(); n < l; n++) {
						if ((*changeNodes)[n] < d) {
							changeV.push_back((*changeNodes)[n]);
						} else {
							changeW.push_back((*changeNodes)[n]-d);
						}
					}
				}

				while (std::abs(norm_this - norm_this_0) > eps)
				{
					stepCount++;
					norm_this_0 = norm_this;
					norm_this = l2norm(*this);
				for (int j = 0, i = changeW.size(); j < i; j++) {
					k = changeW[j];

					std::vector<int> partialSummation = getPartialSummations(k+d);

					Index index2(partialSummation);

					int count = index2.size();

					std::vector<T> matrix(count*count);

					std::vector<T> vector(count);

					Index index3(representation.summations);

					for (Index index1(partialSummation); !index1.end(); ++index1) {
						#pragma omp parallel sections
						{
							#pragma omp section
							{
								for (index2.begin(); !index2.end(); ++index2) {
									matrix[index1*count+index2] =
											tensorScalarProduct(scalarComponentsUU, k+d,
													            index2, index1, (*this));
									// here we generate that matrix entry that is at (index1, index2)
								}
							}
							#pragma omp section
							{
								for (index3.begin(); !index3.end(); ++index3) {
									vector[index1] +=
											representation.tensorScalarProduct(scalarComponentsBU,
													                           k+d, index3, index1,
													                           (*this));
								}
							}
						}
					}
					if (invert_matrix(matrix, count) == 0) {
						// we can copy the matrix as it is into the representation
						Blas<T>::gemv('N', count, count, 1, &matrix[0], count, &vector[0],
								      1, 0, &w[k][0], 1);
					}
				}

				for (int j = 0, i = changeV.size(); j < i; j++) { // k is the node we are optimizing
					k = changeV[j];
					int dimension = componentDimensions[k];

					std::vector<int> partialSummation = getPartialSummations(k);

					Index index2(partialSummation);

					int count = index2.size();

					std::vector<T> matrix(count*count);

					std::vector<T> vector(dimension*count);

					std::vector<int> partialSummationB = representation.getPartialSummations(k);

					Index index3(partialSummationB);

					std::vector<int> partialSummationsA;

					if (subtrahend != 0) {
						partialSummationsA = (*subtrahend).getPartialSummations(k);
					}

					Index index4(partialSummationsA);

					for (Index index1(partialSummation); !index1.end(); ++index1) {
						#pragma omp parallel sections
						{
							#pragma omp section
							{
								// Construction of the H-Matrix
								for (index2.begin(); !index2.end(); ++index2) {
									matrix[index1*count+index2] =
											tensorScalarProduct(scalarComponentsUU, k, index2,
													            index1, (*this));
									// here we generate that matrix entry that is at (index1, index2)
								}
							}

							#pragma omp section
							{
								// Construction of b-vector
								/*
								 * this is a trick - we actually do not need to separate the
								 * summations of the fixed representation
								 */
								for (index3.begin(); !index3.end(); ++index3) {
									const T factor =
											representation.tensorScalarProduct(scalarComponentsBU,
																	k, index3, index1, (*this));

									Blas<T>::axpy(dimension, factor,
											      &representation.v[k][index3*dimension],
											      1, &vector[index1*dimension], 1);
								}
								if (subtrahend != 0) {
									for (index4.begin(); !index4.end(); ++index4) {
										const T factor = -(*subtrahend).tensorScalarProduct(scalarComponentsAU,
																k, index4, index1, (*this));

										Blas<T>::axpy(dimension, factor, &(*subtrahend)[k][index4*dimension],
													  1, &vector[index1*dimension], 1);
									}
								}
							}
						}
					}

					if (invert_matrix(matrix, count) == 0) {
						Blas<T>::gemm('N', 'T', dimension, count, count, 1, &vector[0],
									  &matrix[0], 0, &v[k][0]);
						/*
						 * We can also use
						 * #pragma omp parallel for
						for (int n = 0; n < dimension; n++) {
							Blas<T>::gemv('T', count, count, 1, &matrix[0], count, &vector[n],
									      dimension, 0, &v[k][n], dimension);
						}
						 *
						 * but this is only slightly faster than the above version even w/o
						 * parallelization.
						 *
						 * Is this also true for release code?
						 */

						#pragma omp parallel sections
						{
							#pragma omp section
							{
								tensorScalarProductComponentsDirection((*this),
										scalarComponentsUU[k], k);
							}
							#pragma omp section
							{
								representation.tensorScalarProductComponentsDirection((*this),
										scalarComponentsBU[k], k);
							}
							#pragma omp section
							{
								if (subtrahend != 0) {
									(*subtrahend).tensorScalarProductComponentsDirection((*this),
											scalarComponentsAU[k], k);
								}
							}
						}
						std::cout << stepCount << ". step finished" << std::endl;
					} else {
						std::cout << "SINGULAR" << std::endl;
					}
				}

				}
				return stepCount;
			}

			/*
			 * Updates the k-th rank, but does NOT rearrange the data such that all nodes with
			 * edge k have to be overwritten.
			 */
			void setSummation(int k, int rank, bool reorder = true) {
#ifdef ARGUMENT_CHECKS_ON
				if (k >= summationsCount || k < 0) {
					throw std::invalid_argument("K is out of range.");
				}
				if (rank < 1) {
					throw std::invalid_argument("The rank has to be positive.");
				}
#endif
				for (int n = 0; n < d; n++) {
					if (indexOf(incidenceMatrix[n], k) > -1) {
						v[n].resize((v[n].size()/summations[k])*rank);
						if (reorder && rank > summations[k]) { // rearrange on rank increment
							std::vector<T> v_n(v[n].size());
							std::vector<int> partialSummations = TensorCalculus::getPartialSummations(summations, incidenceMatrix[n]);
							std::vector<int> newSummations(partialSummations);
							newSummations[indexOf(incidenceMatrix[n], k)] = rank;
							for (Index index(partialSummations); !index.end(); ++index) {

								std::memcpy(&v_n[TensorCalculus::vector2index(newSummations, index)*componentDimensions[n]], &v[n][index*componentDimensions[n]], TSIZE*componentDimensions[n]);
							}
							std::memcpy(&v[n][0], &v_n[0], TSIZE*v_n.size());
						}
						partialSummationProducts[n] /= summations[k];
						partialSummationProducts[n] *= rank;
					}
				}
				for (int n = 0; n < L; n++) {
					if (indexOf(incidenceMatrix[n+d], k) > -1) {
						w[n].resize((w[n].size()/summations[k])*rank);
						if (reorder && rank > summations[k]) { // rearrange on rank increment
							std::vector<T> w_n(w[n].size());
							std::vector<int> partialSummations = TensorCalculus::getPartialSummations(summations, incidenceMatrix[n+d]);
							std::vector<int> newSummations(partialSummations);
							newSummations[indexOf(incidenceMatrix[n+d], k)] = rank;
							for (Index index(partialSummations); !index.end(); ++index) {
								w_n[TensorCalculus::vector2index(newSummations, index)] = w[n][index];
							}
							std::memcpy(&w[n][0], &w_n[0], TSIZE*w_n.size());
						}

						partialSummationProducts[n+d] /= summations[k];
						partialSummationProducts[n+d] *= rank;
					}
				}
				summationCountProduct /= summations[k];
				summationCountProduct *= rank;
				summations[k] = rank;
			}

			void setAllSummations(int rank) {
				for (int n = 0; n < summationsCount; n++) {
					setSummation(n, rank);
				}
			}

			/*
			 * Returns the summations of node k.
			 */
			std::vector<int> getPartialSummations(int k) const {
#ifdef ARGUMENT_CHECKS_ON
				if (k >= d+L) {
					throw std::invalid_argument("K is out of range.");
				}
#endif
				return TensorCalculus::getPartialSummations(summations, incidenceMatrix[k]);
			}

			/*
			 * Evaluates the full tensor.
			 */
			std::vector<T> evaluate() const {
				using namespace VectorOperators;

				long size = getComponentProduct(componentDimensions);

				std::vector<T> result(size);

				for (Index index(summations); !index.end(); ++index) {
					std::vector<T> buffer(size);
					std::vector<T> buffer2(size);

					long dimensionProduct = componentDimensions[0];

					Blas<T>::copy(dimensionProduct,
							      &v[0][dimensionProduct*partialVector2index(summations, index,
							    		                                     incidenceMatrix[0])],
							      1, &buffer[0], 1);

					for (int n = 1; n < d; n++) {
						// just for debugging
						const int dimension2 = componentDimensions[n];

						Blas<T>::ger(dimensionProduct, dimension2, 1, &buffer[0], 1,
								     &v[n][dimension2*partialVector2index(summations, index,
								    		                              incidenceMatrix[n])],
								     1, &buffer2[0], dimensionProduct);
						dimensionProduct *= dimension2;
						Blas<T>::copy(dimensionProduct, &buffer2[0], 1, &buffer[0], 1);
						Blas<T>::scal(dimensionProduct, 0, &buffer2[0], 1);
					}

					if (L > 0) {
						T product = 1.0;

						for (int n = 0; n < L; n++) {
							product *= w[n][partialVector2index(summations, index,
									                            incidenceMatrix[n+d])];
						}
						buffer *= product;
					}

					result += buffer;
				}
				return result;
			}

			/*
			 * Returns the storage in bytes that is needed for this tensor representation.
			 */
			long getStorage() const {
				long buffer = 0;

				for (int n = 0; n < d; n++) {
					buffer += (long)v[n].size();
				}
				for (int n = 0; n < L; n++) {
					buffer += (long)w[n].size();
				}
				return buffer*TSIZE;
			}

			void setPointApproximation(const TensorRepresentation<T> &representation) {
#ifdef ARGUMENT_CHECKS_ON
				(*this).checkCompatibility(representation);
#endif
				for (Index index(summations); !index.end(); ++index) {
					std::cout << index.getCurrent() << std::endl;
				}
			}

	};
	
	template<typename T, typename Generator>
	TensorRepresentation<T> createRandomTensorRepresentation
			(const std::vector<int> &summations, const std::vector<int> componentDimensions,
			 const std::vector< std::vector<int> > incidenceMatrix, Generator generator) {
		int d = componentDimensions.size();

		std::vector< std::vector<T> > v(d);

		for (int n = 0; n < d; n++) {
			const long valueCount =
					componentProduct(getPartialSummations(summations, incidenceMatrix[n]))
					* componentDimensions[n];
			v[n].resize(valueCount);
			for (int k = 0; k < valueCount; k++) {
				v[n][k] = generator();
			}
		}
		return TensorRepresentation<T>(summations, v, componentDimensions, incidenceMatrix);
	}

	template<typename T, typename Generator>
	TensorRepresentation<T> createRandomTensorRepresentation
			(int r, const std::vector<int> &componentDimensions,
			 const std::vector< std::vector<int> > &incidenceMatrix, Generator generator) {
		std::vector<int> summationNumbers;

		int d = componentDimensions.size();

		for (int n = 0; n < d; n++) {
			for (int k = 0, i = incidenceMatrix[n].size(); k < i; k++) {
				if (!contains(summationNumbers, incidenceMatrix[n][k])) {
					summationNumbers.push_back(incidenceMatrix[n][k]);
				}
			}
		}

		int i = summationNumbers.size();

		std::vector<int> summations(i);
		for (int n = 0; n < i; n++) {
			summations[n] = r;
		}
		return createRandomTensorRepresentation<T>(summations, componentDimensions,
				                                       incidenceMatrix, generator);
	}

	template<typename T> TensorRepresentation<T> operator +
				(const TensorRepresentation<T>& representation1,
				 const TensorRepresentation<T>& representation2) {
		TensorRepresentation<T> result(representation1);
		
		result += representation2;
		return result;
	}

	template<typename T> TensorRepresentation<T>& operator +=
				(TensorRepresentation<T>& representation1,
				 const TensorRepresentation<T>& representation2) {
		return representation1.add(representation2);
	}

	template<typename T> TensorRepresentation<T> operator *
				(const T& factor, const TensorRepresentation<T>& representation) {
		TensorRepresentation<T> result(representation);
		
		result *= factor;
		return result;
	}

	template<typename T> TensorRepresentation<T>& operator *=
				(TensorRepresentation<T>& representation, const T& factor) {
		return representation.scale(factor);
	}
	
	template<typename T> T tensorScalarProduct(const TensorRepresentation<T> &representation1,
			                                   const TensorRepresentation<T> &representation2) {
		return representation1.tensorScalarProduct(representation2);	
	}
	
	template<typename T> T l2norm(const TensorRepresentation<T> &representation1,
			                      const TensorRepresentation<T> &representation2) {
		using namespace VectorOperators;
		return l2_norm(representation1.evaluate()- representation2.evaluate());
	}
	
	template<typename T> T l2norm(const TensorRepresentation<T> &representation) {
		using namespace VectorOperators;
		return l2_norm(representation.evaluate());
	}
	
	template<typename T> T sqr(T t) {
		return t * t;
	}

	template<typename T> T const_one() {
		return 1.0;
	}

	template<typename T> T componentScalarProduct(const std::vector< std::vector<T> > &v,
			                                      int direction, int i, const T * vector,
						                          const int dimension) {
		return Blas<T>::dot(dimension, vector, 1, &v[direction][i*dimension],1);
	}

	/*
	 * Computes the inner products of this representation
	 * and the given on in the specified direction and stores
	 * the results in the result variable.
	 */
	template <typename T, typename U>
	std::vector<T>& tensorScalarProductComponentsDirection(const std::vector< std::vector<int> > &incidenceMatrix,
			                                               const std::vector<int> &summations,
			                                               const std::vector< std::vector<T> > &v,
			                                               const U &representation,
			                                               std::vector<T> &result,
			                                               /*std::vector<T> &factors,*/
			                                               const int direction,
			                                               const int componentDimension) {
		int summation0 = 1;
		int summation1 = 1;

		for (int k = 0, l = incidenceMatrix[direction].size(); k < l; k++) {
			summation0 *= summations[incidenceMatrix[direction][k]];
		}

		const std::vector<int> incidenceRow = representation.getIncidenceRow(direction);

		for (int k = 0, l = incidenceRow.size(); k < l; k++) {
			summation1 *= representation.getSummation(incidenceRow[k]);
		}

		result.resize(summation0*summation1);

		for (long n = 0; n < summation0; n++) {
			for (long i = 0; i < summation1; i++) {
				result[i*summation0+n] = representation.componentScalarProduct(direction,
						                     i, &(v[direction][n*componentDimension]), componentDimension);
			}
		}
		return result;
	}

	/*
	 * Computes all single inner products of this representation
	 * and the given one and stores them in the result variable.
	 */
	template <typename T, typename U>
	std::vector< std::vector<T> >& tensorScalarProductComponents(
			const std::vector< std::vector<int> > &incidenceMatrix,
            const std::vector<int> &summations,
            const std::vector< std::vector<T> > &v,
            const U &representation,
            std::vector< std::vector<T> >& result,
            const std::vector<int> componentDimensions) {
		int d = v.size();
#ifdef ARGUMENT_CHECKS_ON
		// TODO implement
#endif

		result.resize(d);
		#pragma omp parallel for
		for (int n = 0; n < d; n++) {
			tensorScalarProductComponentsDirection(
					incidenceMatrix,
					summations,
					v,
					representation,
					result[n],
					n,
					componentDimensions[n]);
		}
		return result;
	}

};


#endif /* __TENSORREPRESENTATION_HPP */
