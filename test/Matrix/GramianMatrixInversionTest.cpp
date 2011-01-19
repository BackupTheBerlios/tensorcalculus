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

#include "Matrix/MatrixOperators.hpp"

using namespace TensorCalculus;

using namespace VectorOperators;



int main() {
	int n = 5;

	double alpha = 1.245234314;

	std::vector<double> matrix(n*n);

	for (int i = 0; i < n*n; i++) {
		matrix[i] = ((double) std::rand())/RAND_MAX;
	}

	int dep = n-1;

	for (int i = 0; i < n; i++) {
		matrix[n*dep+i] = (1.0/alpha)*matrix[i];
	}

	std::vector<int> linear_dependence;

	std::cout << matrix << std::endl;
	std::cout << invert_gramian_matrix(matrix, n, linear_dependence) << std::endl;


	std::cout << "Matrix has size " << n-linear_dependence.size() << std::endl;

	std::cout << matrix << std::endl;
	return 0;
}
