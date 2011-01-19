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

#include "CG/ThreePG.hpp"
#include "TestUtilities.hpp"
#include "Constants.hpp"
#include "Utilities/Utilities.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>

namespace TensorCalculus {

  template<>
  struct Constants<double> {
    static const double epsilon;
  };

}

const double TensorCalculus::Constants<double>::epsilon = 1e-8;

using namespace TensorCalculus;

double f(double x) {
  // const double n = 10;
  // return x > 1.0/n ? -1 : -2*n*x + 1;
  return std::cos(x);
}

int main() {

  TEST_FOR_EXCEPTION(ThreePG<double> test(0.25, 0.25), intervall_error);

  ThreePG<double> t(0.75, 0.75);

  // Function is not specified:
  TEST_FOR_EXCEPTION(t.iterate(), std::runtime_error);
  // Invalid intervall:
  TEST_FOR_EXCEPTION(t.initialize(1.0, 0.0, Constants<double>::epsilon, f),
                     intervall_error);
  // Invalid epsilon:
  TEST_FOR_EXCEPTION(t.initialize(0.0, 3.0, -1.0, f),
                     std::invalid_argument);
  // Signs of function at boundaries don't match:
  TEST_FOR_EXCEPTION(t.initialize(0.0, 1.0, 0.1, f),
                     intervall_error);

  t.initialize(0.0, 3.0, Constants<double>::epsilon, f);
  TEST_EQUAL(t.iterations(), 0);

  while (!t.iterate()) { }

  TEST_EQUAL_FLOAT(t.value(), M_PI_2, Constants<double>::epsilon);

  std::cout << "Tests passed!" << std::endl;
  return 0;
}
