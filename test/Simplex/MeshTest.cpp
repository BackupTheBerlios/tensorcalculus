/*
 * Copyright (C) 2010 Philipp Wähnert
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
#include <fstream>

#include "Mesh/Mesh.hpp"

using namespace TensorCalculus;

int main(int argc, const char* argv[]) {
  if (argc < 2) {
    std::cout << "Usage: MeshTest <meshfile>" << std::endl;
    return 0;
  }
  std::ifstream file(argv[1]);
  if (!file) throw StreamFormatError("File doesn't exist");
  Mesh<double> m(file);
  std::cout << m << std::endl;
  
  for (int i = 0; i < m.countSimplices(); ++i) {
    const Mesh<double>::Simplex& simplex = m.getSimplex(i);
    double max_norm = 0.0;
    for (int j = 0; j <= m.getDimension(); ++j) {
      const double norm = Blas<double>::nrm2(m.getDimension(), &simplex.getBaseVector(simplex.getPoint(j))[0], 1);
      if (norm > max_norm) max_norm = norm;
    }
    // std::cout << max_norm << std::endl;
  }
  
  std::vector<double> v = m.getSimplex(0).getRank1ApproxBaseVectors(2);
  // Mesh<double> m2(cross_product(m, m));
  // std::cout << m2 << std::endl;
  return 0;
}
