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

#include "CG/Bisection.hpp"
#include "TestUtilities.hpp"
#include "Constants.hpp"
#include <cmath>
#include <iostream>

namespace TensorCalculus {

  template<>
  struct Constants<double> {
    static const double epsilon;
  };

}

const double TensorCalculus::Constants<double>::epsilon = 1e-8;

using namespace TensorCalculus;

double f(double x) {
  const double n = 10000;
  return x > 1.0/n ? -1 : -2*n*x + 1;
  // return std::cos(x);
}

int main() {

  Bisection<double, double(*)(double)> bisection(0.0, 1.0, &f);

  while(bisection.residuum() > Constants<double>::epsilon) {
    std::cout << "n = " << bisection.iterations()
              << " a = " << bisection.get_a() << " b = " << bisection.get_b()
              << std::endl;
    bisection.iterate();
  }

  // std::cout << M_PI_2 << std::endl;
  std::cout << bisection.value() << std::endl;
  // TEST_EQUAL_FLOAT(bisection.value(), M_PI_2, Constants<double>::epsilon);

  return 0;
}
