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

#include <vector>
#include "Tensor/TensorBasisOperators.hpp"
#include "Tensor/TensorBasis.hpp"
#include "Tensor/CPTensor.hpp"
#include "TestUtilities.hpp"

using namespace TensorCalculus;

int main() {
  CPTensor<double> x(2, 2, 2);
  std::vector<int> n(2);
  n[0] = 3;
  n[1] = 2;
  std::vector<int> m(2);
  m[0] = 2;
  m[1] = 2;
  TensorBasis<double> b(2, n, m);
  std::vector<double>& u = b.getBasis(0);
  u[0] = 1.0;
  u[1] = 0.0;
  u[2] = 0.0;
  u[3] = 0.0;
  u[4] = 1.0;
  u[5] = 1.0;
  std::vector<double>& v = b.getBasis(1);
  v[0] = 1.0;
  v[1] = 0.0;
  v[2] = 1.0;
  v[3] = 1.0;

  return 0;
}
