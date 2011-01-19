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

#include "Representation/CPTensorRepresentation.hpp"

#include <time.h>

using namespace TensorCalculus;

int main() {

	CPTensorRepresentation<double> start("data/tensors/rect0b_start.ten");

	CPTensorRepresentation<double> quad("data/tensors/rect0b_quad.ten");

	CPTensorRepresentation<double> exact("data/tensors/rect0b_exact.ten");

	using namespace VectorOperators;
	std::cout << l2_norm(start.evaluate() - exact.evaluate()) << std::endl;
	clock_t start_clock, end;

	start_clock = clock();

	for (int n = 0; n < 4000; n++) {
		start.performALS(quad);

	}
	 end = clock();
	 std::cout << l2_norm(start.evaluate() - exact.evaluate())/l2_norm(start.evaluate()) << std::endl;
	 std::cout << ((double) end - start_clock) / CLOCKS_PER_SEC << std::endl;



	/*
	for (int n = 0; n < start.getSummation(0); n++) {
		for (int mu = 0; mu < start.getD(); mu++) {
			for (int j = 0; j < start.getComponentDimension(mu); j++) {
				std::cout << start[mu][n*start.getComponentDimension(mu)+j] << std::endl;
			}
		}
	}
	*/

	return 0;
}
