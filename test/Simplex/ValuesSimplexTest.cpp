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

#include "Mesh/ValuesSimplex.hpp"
#include "Utilities/Utilities.hpp"

using namespace TensorCalculus;

const int domain_dim = 3;
const int value_dim = 2;
const double values[domain_dim + 1][value_dim] =
  { { 0.0, 0.0 },
    { 0.5, 0.0 },
    { 0.5, 0.5 },
    { 0.0, 0.5 } };

int main(int argc, const char* argv[])
{
  {
    const int domain_dim = 3;
    const int value_dim = 2;
    const double values[domain_dim + 1][value_dim] =
      { { 0.0, 0.0 },
        { 0.5, 0.0 },
        { 0.5, 0.5 },
        { 0.0, 0.5 } };
    std::vector< std::vector<double> > values_vect(domain_dim + 1);
    for (int i = 0; i <= domain_dim; ++i) {
      values_vect[i].assign(&values[i][0], &values[i][0] + value_dim);
    }
    ValuesSimplex<double> simplex(values_vect, 1);
    std::cout << "Decomposition of ";
    simplex.print(std::cout);
    std::cout << '\n';
    
    std::vector< ValuesSimplex<double> > components = simplex.decompose();
    
    for (unsigned int i = 0; i < components.size(); ++i) {
      components[i].print(std::cout);
      std::cout << '\n';
    }
    std::cout << std::endl;
  }
  {
    const int domain_dim = 1;
    const int value_dim = 2;
    const double values[domain_dim + 1][value_dim] =
      { { 0.0,  1.0 },
        { 1.0, -1.0 } };
    std::vector< std::vector<double> > values_vect(domain_dim + 1);
    for (int i = 0; i <= domain_dim; ++i) {
      values_vect[i].assign(&values[i][0], &values[i][0] + value_dim);
    }
    ValuesSimplex<double> simplex(values_vect, 1);
    std::cout << "Decomposition of ";
    simplex.print(std::cout);
    std::cout << '\n';
    
    std::vector< ValuesSimplex<double> > components = simplex.decompose();
    
    for (unsigned int i = 0; i < components.size(); ++i) {
      components[i].print(std::cout);
      std::cout << '\n';
    }
    std::cout << std::endl;
    
    std::cout << "self-cross product\n";
    std::vector< ValuesSimplex<double> > parts = simplex.cross_product(simplex);
    for (unsigned int i = 0; i < parts.size(); ++i) {
      parts[i].print(std::cout);
      std::cout << '\n';
    }
    std::cout << std::endl;
    
    const int domain_dim2 = 2;
    const int value_dim2 = 1;
    const double values2[domain_dim2 + 1][value_dim2] =
      { { 0.0 },
        { 1.0 },
        { 2.0 } };
    std::vector< std::vector<double> > values_vect2(domain_dim2 + 1);
    for (int i = 0; i <= domain_dim2; ++i) {
      values_vect2[i].assign(&values2[i][0], &values2[i][0] + value_dim2);
    }
    ValuesSimplex<double> simplex2(values_vect2, 1);
    
    std::cout << "cross product with ";
    simplex2.print(std::cout);
    std::cout << '\n';
    std::vector< ValuesSimplex<double> > blub = simplex.cross_product(simplex2); // simplex2.cross_product(simplex);
    for (unsigned int i = 0; i < blub.size(); ++i) {
      blub[i].print(std::cout);
      std::cout << '\n';
    }
    std::cout << std::endl;
  }
  
  return 0;
}
