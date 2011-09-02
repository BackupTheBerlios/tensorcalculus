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

#include "Utilities/Random.hpp"
#include "Matrix/MatrixOperators.hpp"

using namespace TensorCalculus;
int main() {

	int n = 20, m = 13; // n = colcount, m = rowrount

	int rank = std::min(n,m);

	std::vector<double> cols(rank*m);
	std::vector<double> rows(rank*n);
	std::vector<double> matrix(n*m);

	Random<double> random;
  for (int k = 0; k < n*m; k++) {
		matrix[k] = random();
	}

	int eff_rank = aca(matrix, m, n, 1e-3, cols, rows);
	return 0;
}
