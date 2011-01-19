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
#include <string>

#include "Mesh/Mesh.hpp"
#include "Mesh/MeshFunction.hpp"
#include "Vector/VectorOperators.hpp"
#include "Utilities/Factorials.hpp"

using namespace TensorCalculus;

int main(int argc, const char* argv[]) {
  if (argc != 2) {
    std::cout << "Usage: MeshFunctionTest <mesh>" << std::endl;
    return 0;
  }
  std::string basename(argv[1]);
  
  std::string mesh_filename(basename + ".txt");
  std::ifstream mesh_file(mesh_filename.c_str());
  if (!mesh_file) {
    std::cout << "Can't open mesh file \"" << mesh_filename << '"' << std::endl;
    return -1;
  }
  Mesh<double> mesh(mesh_file);
  
  std::string meshfunction_filename(basename + ".func");
  std::ifstream meshfunction_file(meshfunction_filename.c_str());
  if (!meshfunction_file) {
    std::cout << "Can't open function file \"" << meshfunction_filename << '"' << std::endl;
    return -2;
  }
  MeshFunction<double> mesh_function(meshfunction_file);
  
  double max_error = 0.0;
  const double factor = factorial<double>(mesh.getDimension());
  for (int i = 0; i < mesh.countSimplices(); ++i) {
    const Mesh<double>::Simplex& simplex = mesh.getSimplex(i);
    double max_norm = 0.0;
    double g = 0.0;
    for (int j = 0; j <= mesh.getDimension(); ++j) {
      const double norm = Blas<double>::nrm2(mesh.getDimension(), &simplex.getBaseVector(simplex.getPoint(j))[0], 1);
      if (norm > max_norm) max_norm = norm;
      
      g += Blas<double>::asum(mesh_function.countFunctions(), &mesh_function.getValues(simplex.getPoint(j))[0], 1);
    }
    const double error = simplex.getDeterminant() * simplex.getDiameter() * 10 * g * max_norm / factor;
    if (error > max_error) max_error = error;
    // std::cout << error << std::endl;
  }

  std::cout << max_error << std::endl;
  
  return 0;
}
