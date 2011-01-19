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

#include "Representation/TensorChainRepresentation.hpp"
#include "Utilities/Utilities.hpp"

using namespace TensorCalculus;

using namespace VectorOperators;

int main() {

	int d = 4;

	int n = 3;

	int r = 3;

	TensorChainRepresentation<double> tc1 = createRandomTensorChain(r-1, d, n, Utilities<double>::rand);

	TensorChainRepresentation<double> tc2 = createRandomTensorChain(r, d, n, Utilities<double>::rand);

	std::cout << "Original norm=" << l2norm(tc1) << std::endl;
	std::cout << "Original norm=" << l2_norm(tc1.evaluate() -  tc2.evaluate()) << std::endl;
	tc1.performDMRG(tc2);
	//while(true) {
	//tc1.performALS(tc2);
	std::cout << "Original norm2=" << l2_norm(tc1.evaluate() -  tc2.evaluate()) << std::endl;
	//}
	return 0;
}
