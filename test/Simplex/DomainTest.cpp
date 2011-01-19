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

#include "Mesh/Domain.hpp"
#include "Mesh/Mesh.hpp"
#include "Mesh/MeshFunction.hpp"

using namespace TensorCalculus;

int main(int argc, const char* argv[]) {
  if (argc < 2 || argc > 3) {
    std::cout << "Usage: DomainTest <meshfile> [<meshfunction>]" << std::endl;
    return 0;
  }
  std::ifstream mesh_file(argv[1]);
  if (!mesh_file) {
    std::cout << "Can't open meshfile" << std::endl;
    return -1;
  }
  Mesh<double> mesh(mesh_file);

  if (argc == 3) {
    std::ifstream mesh_function_file(argv[2]);
    if (!mesh_function_file) {
      std::cout << "Can't open meshfile" << std::endl;
      return -1;
    }
    MeshFunction<double> mesh_function(mesh_function_file);
    Domain<double> domain(mesh, mesh_function);
    Domain<double> cd = cross_product(domain, domain);
    return 0;
  } else {
    Domain<double> domain(mesh);
    return 0;
  }
}
