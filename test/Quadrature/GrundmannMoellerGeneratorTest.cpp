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
#include <vector>
#include <cmath>
#include <algorithm>

#include "Quadrature/GrundmannMoellerGenerator.hpp"
#include "Vector/VectorOperators.hpp"

using namespace TensorCalculus;

const int d = 2;

const int accuracy = 10;
const int null_rules = 2;

const int exponent = 100;
// const int exponents[d] = { 10, 5 };


double integral() {
  double result = 1;
  for (int i = 1; i <= d; ++i) {
    result *= (exponent + i);
  }
  return 1.0/result;
}

double func(const std::vector<double>& x) {
  return std::pow(x[0], exponent);
}

int main() {
  GrundmannMoellerGenerator<double> g(d);

  std::vector< Node<std::vector<double>, double > > nodes = g.nodes(accuracy);

  double value = 0.0;
  double sum = 0.0;

  std::cout << "Nodes: " << nodes.size() << '\n';
  for (int i = 0; i < nodes.size(); ++i) {
    const Node<std::vector<double>,double >& node = nodes[i];
    
    // using VectorOperators::operator <<;
    // std::cout << node.getAbscissa() << " ";

    value += node.getWeight() * func(node.getAbscissa());
    sum += node.getWeight();
  }

  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout.precision(20);
  std::cout << "Sum: " << sum << '\n';
  std::cout << "Error: " << (sum - 0.5) / 0.5 << '\n';
  std::cout << "Value: " << value << '\n';
  std::cout << "Integral: " << integral() << std::endl;
  return 0;
}
