/*
 * Copyright (C) 2010 Stefan Handschuh
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
#include <time.h>
#include "Representation/TensorRepresentation.hpp"
#include "Representation/TensorChainRepresentation.hpp"
#include "Utilities/Index.hpp"

void testIndexFunctions();

void testScalarProduct();

void testTC();

int main() {
	testTC();
	return 0;
}

void testTC() {

	// CP-Tensor with rank 2, component dimension 3 and a 4 fold tensor product
	
	std::vector< std::vector<double> > x1(4);
	std::vector< std::vector<double> > x2(4);
	srand(time(NULL));
	for (int n = 0; n < 4; n++) {
		x1[n].resize(6);  // each summand-component has dimension 3 and there are 2 summand-components
		x2[n].resize(6);
		for (int i = 0; i < 6; i++) {
			x1[n][i] = (double) rand() / (double) RAND_MAX;
			x2[n][i] = (double) rand() / (double) RAND_MAX;	
		}	
	}
	
	
	TensorCalculus::TensorRepresentation<double> cpRepresentation1(2, x1, 3);	
	TensorCalculus::TensorRepresentation<double> cpRepresentation2(2, x2, 3);
	//cpRepresentation1.toDot();

	/*
	TensorCalculus::TensorRepresentation<double> cpRepresentation3 = cpRepresentation1 + cpRepresentation2;
	
	for (int n = 0; n < 4; n++) {
		for (int i = 0; i < 6; i++) {
			if (cpRepresentation1[n][i] + cpRepresentation2[n][i] != cpRepresentation3[n][i]) {
				throw std::invalid_argument("Incorrect summation.");	
			}
		}	
	}
	
	*/
	// Tensor Chain format
	int chainLength = 5;
	
	int componentDimension = 3;
	
	int summation = 2;
	
	std::vector< std::vector<double> > x_tc(chainLength);
	
	std::vector<int> summations(chainLength);
	
	std::vector< std::vector<int> > incidenceMatrix(chainLength);
	
	std::vector<int> componentDimensions(chainLength);
	
	for (int n = 0; n < chainLength; n++) {
		x_tc[n].resize(componentDimension*pow(summation,2));
		
		incidenceMatrix[n].resize(2);
		incidenceMatrix[n][0] = n;
		incidenceMatrix[n][1] = (1+n) % chainLength;
		
		summations[n] = summation;
		
		componentDimensions[n] = componentDimension;
	}
	

	//TensorCalculus::TensorRepresentation<double> tcRepresentation(summations, x_tc, componentDimensions, incidenceMatrix);
	

	// PEPS format
	std::vector< std::vector<double> > v_peps(5);
	componentDimensions = std::vector<int>(5);
	summations = std::vector<int>(12);
	
	for (int n = 0; n < 12; n++) {
		summations[n] = summation;
	}
	
	v_peps[0].resize(pow(summation,2)*componentDimension);
	v_peps[1].resize(pow(summation,3)*componentDimension);
	v_peps[2].resize(pow(summation,2)*componentDimension);
	v_peps[3].resize(pow(summation,3)*componentDimension);
	v_peps[4].resize(pow(summation,2)*componentDimension);
	
	for (int n = 0; n < 5; n++) {
		componentDimensions[n] = componentDimension;
	}
	
	std::vector< std::vector<double> > w_peps(4);
	
	w_peps[0].resize(pow(summation,3));
	w_peps[1].resize(pow(summation,4));
	w_peps[2].resize(pow(summation,2));
	w_peps[3].resize(pow(summation,3));
	
	incidenceMatrix.resize(9);
	incidenceMatrix[0].resize(2);
	incidenceMatrix[0][0] = 0;
	incidenceMatrix[0][1] = 7;
	incidenceMatrix[1].resize(3);
	incidenceMatrix[1][0] = 0;
	incidenceMatrix[1][1] = 8;
	incidenceMatrix[1][2] = 1;
	incidenceMatrix[2].resize(2);
	incidenceMatrix[2][0] = 1;
	incidenceMatrix[2][1] = 2;
	incidenceMatrix[3].resize(3);
	incidenceMatrix[3][0] = 2;
	incidenceMatrix[3][1] = 9;
	incidenceMatrix[3][2] = 3;
	incidenceMatrix[4].resize(2);
	incidenceMatrix[4][0] = 3;
	incidenceMatrix[4][1] = 4;
	incidenceMatrix[5].resize(3);
	incidenceMatrix[5][0] = 7;
	incidenceMatrix[5][1] = 11;
	incidenceMatrix[5][2] = 6;
	incidenceMatrix[6].resize(4);
	incidenceMatrix[6][0] = 8;
	incidenceMatrix[6][1] = 9;
	incidenceMatrix[6][2] = 10;
	incidenceMatrix[6][3] = 11;
	incidenceMatrix[7].resize(2);
	incidenceMatrix[7][0] = 6;
	incidenceMatrix[7][1] = 5;
	incidenceMatrix[8].resize(3);
	incidenceMatrix[8][0] = 5;
	incidenceMatrix[8][1] = 10;
	incidenceMatrix[8][2] = 4;
	
	TensorCalculus::TensorRepresentation<double> pepsRepresentation(summations, v_peps, componentDimensions, incidenceMatrix, w_peps);
	
	//pepsRepresentation.toDot();
	
	//TensorCalculus::TensorChainRepresentation<double> tensorChain1(4,2);
	
	//TensorCalculus::TensorChainRepresentation<double> tensorChain2(4,2);
	
	//tensorChain1 += tensorChain2;
	//testIndexFunctions();


}

void testIndexFunctions() {
	std::vector<int> summations(3);
	
	summations[0] = 5;
	summations[1] = 3;
	summations[2] = 7;
	
	std::vector<int> position(3);
	
	position[0] = 4;
	position[1] = 1;
	position[2] = 1;
	
	std::vector<int> partialSummations(2);
	
	partialSummations[0] = 0;
	partialSummations[1] = 1;
	
	std::vector<int> partialPosition(2);
	partialPosition[0] = 4;
	partialPosition[1] = 1;
	
	std::cout << TensorCalculus::partialVector2index(partialPosition, summations, partialSummations) << std::endl;
	//std::cout << partialPosition <<std::endl;
	std::cout << TensorCalculus::vector2index(summations, position) << std::endl;
	
	TensorCalculus::index2vector(summations, 9);
	std::cout << TensorCalculus::vector2index(summations, TensorCalculus::index2vector(summations, TensorCalculus::vector2index(summations, position))) <<std::endl;
	
}

void testScalarProduct() {
	int length = 3;
	
	std::vector< std::vector<double> > x1(length);
	std::vector< std::vector<double> > x2(length);
	
	int dimension = 10;
	
	int summation = 10;
	
	srand(time(NULL));
	
	std::vector<int> summations(length);
	
	for (int n = 0; n < length; n++) {
		x1[n].resize(dimension*summation*summation);  
		x2[n].resize(dimension*summation*summation);
		summations[n] = summation;
		
		for (int i = 0; i < dimension*summation*summation; i++) {
			x1[n][i] = (double) rand() / (double) RAND_MAX;	
			x2[n][i] = x1[n][i] ;//(double) rand() / (double) RAND_MAX;	
		}	
	}
	TensorCalculus::TensorChainRepresentation<double> tcRepresentation1(summation, dimension, x1);	
	TensorCalculus::TensorChainRepresentation<double> tcRepresentation2(summation, dimension, x2);
	
	double checkResult = 0;
	
	std::vector<int> partialSummations(2);
	partialSummations[0] = summation*summation;
	partialSummations[1] = summation*summation;
	std::vector<int> dummySummations(2);
	dummySummations[0] = summation;
	dummySummations[1] = summation;
	
	
	
	std::vector<int> p1(2);
	std::vector<int> p2(2);
	
	double temp1, temp2, temp3, product = 0;
	
	float c1 = clock();
	/*
	int count = std::pow(summation, length);
	
	std::vector< std::vector<double> > scalarComponents = tcRepresentation1.tensorScalarProductComponents(tcRepresentation2);
	
	std::vector<int> index1, index2;
	
	std::vector<int> position(2);
	
	for (Index index1(tcRepresentation1.getSummations()); !index1.end(); ++index1) {
		for (Index index2(tcRepresentation2.getSummations()); !index2.end(); ++index2) {
			product = 1;
			for (int mu = 0; mu < length; mu++) {
				position[0] = TensorCalculus::partialVector2index(tcRepresentation1.getSummations(), index1, tcRepresentation1.getIncidenceMatrix(), mu);
				position[1] = TensorCalculus::partialVector2index(tcRepresentation2.getSummations(), index2, tcRepresentation2.getIncidenceMatrix(), mu);
							//std::cout << vector2string(position) << " " << scalarComponents[mu][vector2index(partialSummations, position)] << "*";
				product *= scalarComponents[mu][TensorCalculus::vector2index(partialSummations, position)];
			}
			checkResult += product;
		}	
	}
	std::cout <<  std::endl << checkResult <<  " in " << clock() - c1 <<  std::endl;
	
	c1 = clock();
	*/
	checkResult = 0;
	
	double temp_res;
	
	for (int n = 0; n < summation; n++) {
		
		for (int k = 0; k < summation; k++) {
			p1[0] = n;
			p1[1] = k;
			int index1 = TensorCalculus::vector2index(dummySummations, p1);
			
			for (int l = 0; l < summation; l++) {
				p1[0] = l;
				p1[1] = n;
				int index2 = TensorCalculus::vector2index(dummySummations, p1);
				p1[0] = k;
				p1[1] = l;
				int index3 = TensorCalculus::vector2index(dummySummations, p1);
				for (int n2 = 0; n2 < summation; n2++) {
					for (int k2 = 0; k2 < summation; k2++) {
						
						p2[0] = n2;
						p2[1] = k2;
						temp2 = TensorCalculus::Blas<double>::dot(dimension, &(x1[1][index1*dimension]), 1, &(x2[1][TensorCalculus::vector2index(dummySummations, p2)*dimension]), 1);
						temp_res = 0;
						for (int l2 = 0; l2 < summation; l2++) {
							
							p2[0] = l2;
							p2[1] = n2;
							temp1 = TensorCalculus::Blas<double>::dot(dimension, &(x1[0][index2*dimension]), 1, &(x2[0][TensorCalculus::vector2index(dummySummations, p2)*dimension]), 1);
							
							p2[0] = k2;
							p2[1] = l2;
							temp3 = TensorCalculus::Blas<double>::dot(dimension, &(x1[2][index3*dimension]), 1, &(x2[2][TensorCalculus::vector2index(dummySummations, p2)*dimension]), 1);
							
							temp_res += temp1*temp3;
						}
						checkResult += temp_res*temp2;			
					}
				}
			}
		}
	} 
	std::cout << "Norm = " << TensorCalculus::l2norm(tcRepresentation1) << "   " << TensorCalculus::tensorScalarProduct(tcRepresentation1, tcRepresentation2)*TensorCalculus::tensorScalarProduct(tcRepresentation1, tcRepresentation2) << std::endl;
	std::cout <<  std::endl << checkResult <<  " in " << clock() - c1 <<  std::endl;
	for (int n = 0; n < 20; n++) {
		c1 = clock();
		checkResult = TensorCalculus::tensorScalarProduct(tcRepresentation1, tcRepresentation2);
		std::cout << checkResult << " in " << clock() - c1 << std::endl;
	}
	

}
