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

/*
 * This is a test for checking the tensorrepresentation algorithms
 * for tensors that have a component tensor.
 */


#include "Representation/TensorRepresentation.hpp"
#include "Representation/CPTensorRepresentation.hpp"
#include "Representation/TensorChainRepresentation.hpp"
#include "Utilities/Random.hpp"
#include "Vector/VectorOperators.hpp"
#include "DGKTSDCG.hpp"
#include "Representation/MPSRepresentation.hpp"

using namespace TensorCalculus;

using namespace VectorOperators;

int main() {

	double real_values[3][2][2]
	= { { { 1.0, 8.0 / (25.0 * M_PI * M_PI) },
	{ -8.0 / (9.0 * M_PI * M_PI), -2.0 / (M_PI * M_PI) } },
	{ { 8.0 / (M_PI * M_PI), 2.0 / (9.0 * M_PI * M_PI) },
	{ -2.0 / (M_PI * M_PI), -9080.0 / (3969.0 * M_PI * M_PI) } },
	{ { (M_PI * M_PI + 4.0) / (2.0 * M_PI * M_PI), 4664.0 / (11025.0 * M_PI * M_PI) },
	{ -568.0 / (225.0 * M_PI * M_PI), -(225.0 * M_PI * M_PI + 1936.0) / (1800.0 * M_PI * M_PI) } } };

	std::vector<int> dims(3);

	dims[0] = 3;
	dims[1] = 2;
	dims[2] = 2;
	std::vector<double> values(3*2*2);

	for (int n = 0; n < 3; n++) {
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				values[j*2*3 + i*3 +n] = real_values[n][i][j];
			}
		}
	}

	CPTensorRepresentation<double> cp0(dims, values);

	CPTensorRepresentation<double> cp4("data/tensors/quadappr_wow.ten");

	for (int n1 = 0; n1 < cp4.getD(); n1++) {
		for (int n2 = n1+1; n2 < cp4.getD(); n2++) {
			std::cout << "(" << n1 << ", " <<n2 << ") = " << std::endl;
		}
	}



	DKTS dkts4, dkts5;

	dkts4.readDataFrom("data/tensors/quadappr_wow.ten");

	DGKTSDCG cg4;
	//std::cout << l2norm(cp4) << std::endl;

	cg4.setAccuracy(1e-10);
	//cg4.truncate2Eps(1e-7, dkts4, dkts5);


	TensorRepresentation<double> cp_p = createRandomCPTensor<double>(50, cp4.getD(), std::vector<int>(cp4.getComponentDimensions()), TensorCalculus::Random<double>(0.0, 0.0));


	int f = 0;
	//cp0.performALS(cp0);
	//std::cout << ++f << '\t' << l2_norm(cp0.evaluate() - cp_p.evaluate()) / l2_norm(cp0.evaluate()) << std::endl;

	TensorRepresentation<double> tc0("data/tensors/TensorChain/berlin.tc");

	std::vector<int> nSummations(tc0.getSummations());

	for (unsigned int n = 0; n < nSummations.size(); n++) {
		nSummations[n]++;
	}
	TensorRepresentation<double> tc1 = createRandomTensorRepresentation<double>(nSummations,
			   tc0.getComponentDimensions(), tc0.getIncidenceMatrix(), TensorCalculus::Random<double>(0.0, 0.0));


	CPTensorRepresentation<double> cp1("data/tensors/H2O STO-3G/aoint.eps.ten");

	for (int n1 = 0; n1 < cp1.getD(); n1++) {
		for (int n2 = n1+1; n2 < cp1.getD(); n2++) {
			std::cout << "(" << n1 << ", " <<n2 << ") = " << std::endl;
		}
	}
	cp1.computeOrthonormalBasis();

	DKTS dkts, dkts2;

	dkts.readDataFrom("data/tensors/H2O 6-31G/v_ijkl.ten");

	DGKTSDCG cg;

	cg.setAccuracy(1e-5);
	cg.truncate2Eps(1e-8, dkts, dkts2);
	//TensorRepresentation<double> cp1 = createRandomCPTensor(rank, d, componentDimensions /* 'dimension' is also possible */, Utilities<double>::rand);

	int d = cp1.getD();

	CPTensorRepresentation<double> cp2 = createRandomCPTensor<double>(2, d, cp1.getComponentDimensions(), TensorCalculus::Random<double>());

	std::cout << cp2[0] << std::endl;

	std::vector<int> componentDimensions(d);

	for (int n = 0; n < d; n++) {
		componentDimensions[n] = cp1.getComponentDimension(n);
	}

	int rank;
	std::cout << "rank=";
	std::cin >> rank;

	std::vector< std::vector<double> > v(d);

	TensorCalculus::Random<double> random;
  
  for (int n = 0; n < d; n++) {
		v[n].resize(rank*componentDimensions[n]);
		for (int r = 0, j = v[n].size(); r < j; r++) {
			v[n][r] = random();
		}
	}

	int L = 1; // for the beginning;

	std::vector< std::vector<double> > w(L);

	for (int n = 0; n < L; n++) {
		w[n].resize(rank*rank);
		for (int r = 0; r < rank*rank; r++) {
			w[n][r] = random();
		}
	}

	std::vector<int> summations(2);
	summations[0] = rank;
	summations[1] = rank;

	std::vector< std::vector<int> > incidenceMatrix(d+L);
	for (int n = 0; n < d; n++) {
		incidenceMatrix[n].resize(1);
	}
	incidenceMatrix[0][0] = 1;
	incidenceMatrix[2][0] = 1;
	if (L > 0) {
		incidenceMatrix[4].resize(2);
		incidenceMatrix[4][1] = 1;
	}







	TensorRepresentation<double> componentTensorRepresentation(summations, v, componentDimensions, incidenceMatrix, w);

	std::cout << "Storage: cp=" << cp1.getStorage() << "\t\t ct=" <<  componentTensorRepresentation.getStorage()<< std::endl;

	std::cout << "Norm   : cp=" << l2_norm(cp1.evaluate()) << "\t ct=" << l2_norm(componentTensorRepresentation.evaluate()) << std::endl;

	std::cout << "RDist. : " << l2_norm(cp1.evaluate() - componentTensorRepresentation.evaluate()) / l2_norm(cp1.evaluate()) << " \twith norm=" << l2_norm(componentTensorRepresentation.evaluate()) << std::endl;


	//componentTensorRepresentation.performALS(cp1, 1e-7);
	//std::cout << "RDist. : " << l2_norm(cp1.evaluate() - componentTensorRepresentation.evaluate()) / l2_norm(cp1.evaluate()) << " \twith norm=" << l2_norm(componentTensorRepresentation.evaluate()) << std::endl;

	componentTensorRepresentation.performALS(cp1, 1e-8);
	std::cout << "RDist. : " << l2_norm(cp1.evaluate() - componentTensorRepresentation.evaluate()) / l2_norm(cp1.evaluate()) << " \twith norm=" << l2_norm(componentTensorRepresentation.evaluate()) << std::endl;


	MPSRepresentation<double> mps;

	mps.read("data/tensors/MPStensor/h20_6-31G_mpsChem.MPSten");
	//cp1.toMat('a');
	//componentTensorRepresentation.toMat();
	//for (rank = 20; rank < 25; rank++)
	componentDimensions = mps.getComponentDimensions();
	rank = 12;
	{
		std::vector< std::vector<double> > v(d);

		for (int n = 0; n < d; n++) {
			v[n].resize(rank*componentDimensions[n]);
			for (int r = 0, j = v[n].size(); r < j; r++) {
				v[n][r] = random();
			}
		}

		int L = 1; // for the beginning;

		std::vector< std::vector<double> > w(L);

		for (int n = 0; n < L; n++) {
			w[n].resize(rank*rank);
			for (int r = 0; r < rank*rank; r++) {
				w[n][r] = random();
			}
		}

		std::vector<int> summations(2);
		summations[0] = rank;
		summations[1] = rank;

		std::vector< std::vector<int> > incidenceMatrix(d+L);
		for (int n = 0; n < d; n++) {
			incidenceMatrix[n].resize(1);
		}
		incidenceMatrix[2][0] = 1;
		incidenceMatrix[3][0] = 1;
		if (L > 0) {
			incidenceMatrix[4].resize(2);
			incidenceMatrix[4][1] = 1;
		}


		TensorRepresentation<double> componentTensorRepresentation(summations, v, componentDimensions, incidenceMatrix, w);
		std::cout << componentTensorRepresentation.getComponentDimensions() << std::endl;
		componentTensorRepresentation.computeOrthonormalBasis();
				std::cout << componentTensorRepresentation.getComponentDimensions() << std::endl;
		componentTensorRepresentation.performALS(mps);
		//componentTensorRepresentation.performALS(mps);

		//cp2.performALS(cp1);

		std::vector<double> f = mps.evaluate();

		std::cout << rank << "\t" << l2_norm(f - componentTensorRepresentation.evaluate()) / l2_norm(f) << std::endl;
	}
	//componentTensorRepresentation.toMat();
}
