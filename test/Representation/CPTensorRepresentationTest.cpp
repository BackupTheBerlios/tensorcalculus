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

#include "Representation/CPTensorRepresentation.hpp"
#include "Representation/TensorChainRepresentation.hpp"
#include "Representation/TensorTrainRepresentation.hpp"
#include "LapackInterface2.hpp"
#include "DKTS.hpp"
#include "DGKTSDCG.hpp"
#include "Utilities/Utilities.hpp"

using namespace TensorCalculus;




void testALS() {
	using namespace VectorOperators;

	TensorCalculus::CPTensorRepresentation<double> cpTensor("data/tensors/H2O 6-31G/v_ijkl.ten");

	TensorRepresentation<double> tensorChain = createRandomCPTensor(72, cpTensor.getD(), cpTensor.getComponentDimensions(), Utilities<double>::rand);
	
	double norm = l2norm(cpTensor);

	std::cout << "Vorher  " << l2norm(tensorChain, cpTensor) /norm << std::endl;

	while (true) {
	tensorChain.performALS(cpTensor);
	std::cout << "TensorChain=" << l2_norm(tensorChain.evaluate() -  cpTensor.evaluate()) / norm << std::endl;
	}std::cout << "Storage: TC=" << tensorChain.getStorage() << ", CP=" << cpTensor.getStorage() << std::endl;
	std::cout << "Nachher " << l2norm(tensorChain, cpTensor) /norm << std::endl;

}

/*
void testDMRG() {
	using namespace VectorOperators;

	TensorCalculus::CPTensorRepresentation<double> cpTensor("data/tensors/H2O 6-31G/v_ijkl.ten");

	std::vector<int> componentDimensions(cpTensor.getComponentDimensions());

	TensorRepresentation<double> tensorChain = createRandomTensorChain(1, componentDimensions, Utilities<double>::rand);

	TensorTrainRepresentation<double> tensorTrain = createRandomTensorTrain(1, componentDimensions, Utilities<double>::rand);

	double norm = cpTensor.normalize();

	//cpTensor.scale(100);
	std::cout << " ----- DMRG -----" << std::endl;
	std::cout << "Vorher  " << l2norm(tensorTrain, cpTensor) << std::endl;
	tensorTrain.performDMRG(cpTensor);
	tensorTrain.performDMRG(cpTensor);
	std::cout << " ------ ALS  -----" << std::endl;
	while (true) {
		tensorChain.performALS(cpTensor);
		tensorChain.normalize();
	//std::cout << "TensorChain=" << l2_norm(tensorChain.evaluate() -  cpTensor.evaluate()) << std::endl;
	std::cout << "Storage: TC=" << tensorTrain.getStorage() << ", CP=" << cpTensor.getStorage() << std::endl;
	std::cout << "Nachher " << l2norm(tensorTrain, cpTensor) << std::endl;
	}
	std::cout << tensorChain.evaluate()- cpTensor.evaluate() << std::endl;
}
*/
int main() {

	DKTS dkts;

		dkts.readDataFrom("data/tensors/H2O 6-31G/v_ijkl.ten");

		DKTS dkts2(dkts.d(), 8, dkts.n());

		DGKTSDCG ka;

		std::cout << "DKTS-Norm0=" << dkts2.distanceTo(dkts) << std::endl;

		//dkts2.successiveApproximation(dkts, 10e-7);



		//ka.truncate2Eps(10e-7, dkts, dkts2, true, false);


		std::cout << "DKTS-Norm1=" << dkts2.distanceTo(dkts) << ", " << dkts2.k() << std::endl;


//		testDMRG();
	testALS();
	//
	/*
	TensorCalculus::CPTensorRepresentation<double> cpTensor("test.ten");
	
	
	
	
	TensorCalculus::TensorChainRepresentation<double> tensorChain(cpTensor);
	
	*/
	/*
	double sum = 0;
	
	for (Index index(cpTensor.getComponentDimensions()); !index.end(); ++index) {
		sum += tensorChain[index] - cpTensor[index];
		std::cout << tensorChain[index] << "   " << cpTensor[index] << std::endl;
		
	}
	*/
	
	/*
	int d = cpTensor.getD();
	
	std::vector< std::vector<double> > v(d);
	
	std::vector<int> summations(d);
	
	int summation = 4;
	
	for (int n = 0; n < d; n++) {
		summations[n] = 4;
		v[n].resize(summation*summation*cpTensor.getComponentDimensions()[n]);
		for (int k = 0; k < summation*summation*cpTensor.getComponentDimensions()[n]; k++) {
			v[n][k] = 1;
		}
	}
	TensorCalculus::TensorChainRepresentation<double> alsTc(summations, cpTensor.getComponentDimensions(), v);
	std::cout << TensorCalculus::l2norm(cpTensor) << " norm = " << TensorCalculus::l2norm(alsTc, cpTensor) << std::endl;
	for (int n = 0; n < 30; n++) {
		alsTc.performALS(cpTensor, 0);
		std::cout << TensorCalculus::l2norm(alsTc) << " norm = " << TensorCalculus::l2norm(alsTc, cpTensor) << std::endl;
	}
	*/
	
	return 0;
	
}
