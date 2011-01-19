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

#include "Tensor/CPTensor.hpp"
#include "Tensor/CPTensorOperators.hpp"
#include "Vector/VectorOperators.hpp"
#include "Constants.hpp"
#include "TestUtilities.hpp"
#include <cmath>

namespace TensorCalculus {

  template<>
  struct Constants<double> {
    static const double epsilon;
  };

}

const double TensorCalculus::Constants<double>::epsilon = 1e-8;

using namespace TensorCalculus;

class PrintTensor {
public:
  PrintTensor(const CPTensor<double>& x) : x(x) { }
  
  void operator () (const std::vector<int>& i) const {
    using VectorOperators::operator <<;
    std::cout << i << ": " << x.getTensorEntry(i) << std::endl;
  }
  
private:
  const CPTensor<double>& x;
};

template<typename F>
void test_tensor(const CPTensor<double>& x,
                 F& function) {
  const int d = x.get_d();
  const std::vector<int>& n = x.get_n();
  
  std::vector<int> i(d, 0);
  while (i[d - 1] < n[d - 1]) {
    function(i);
    ++i[0];
    for (int mu = 0; mu < d - 1; ++mu) {
      if (i[mu] >= n[mu]) {
        i[mu] = 0;
        ++i[mu + 1];
      }
    }
  }
}

int main() {
  std::vector<int> n(2);
  n[0] = 3;
  n[1] = 2;
  CPTensor<double> t(3, 2, n);

  std::vector<double>& v = t.getVectorOfDimension(0);
  v[0] = 1.0;
  v[4] = 1.0;
  v[8] = 1.0;

  std::vector<double>& w = t.getVectorOfDimension(1);
  w[0] = 1.0;
  w[2] = 1.0;
  w[5] = 1.0;

  CPTensor<double> s(t);
  scale(s, 4.0);
  TEST_EQUAL(s(0, 0, 0), 2.0);
  TEST_EQUAL(s(1, 0, 1), 2.0);
  TEST_EQUAL(s(2, 0, 2), 2.0);
  TEST_EQUAL(s(0, 1, 0), 2.0);
  TEST_EQUAL(s(1, 1, 0), 2.0);
  TEST_EQUAL(s(2, 1, 1), 2.0);

  CPTensor<double> u;
  hadamardProduct(t, s, u);
  TEST_EQUAL(u(0, 0, 0), 2.0);
  TEST_EQUAL(u(0, 0, 1), 0.0);
  TEST_EQUAL(u(0, 0, 2), 0.0);
  TEST_EQUAL(u(0, 1, 0), 2.0);
  TEST_EQUAL(u(0, 1, 1), 0.0);
  TEST_EQUAL(u(1, 0, 0), 0.0);
  TEST_EQUAL(u(1, 0, 1), 0.0);
  TEST_EQUAL(u(1, 0, 2), 0.0);  
  TEST_EQUAL(u(1, 1, 0), 2.0);
  TEST_EQUAL(u(1, 1, 1), 0.0);
  TEST_EQUAL(u(2, 0, 0), 0.0);
  TEST_EQUAL(u(2, 0, 1), 0.0);
  TEST_EQUAL(u(2, 0, 2), 0.0);
  TEST_EQUAL(u(2, 1, 0), 0.0);
  TEST_EQUAL(u(2, 1, 1), 0.0);
  TEST_EQUAL(u(8, 0, 0), 0.0);
  TEST_EQUAL(u(8, 0, 1), 0.0);
  TEST_EQUAL(u(8, 0, 2), 2.0);
  TEST_EQUAL(u(8, 1, 0), 0.0);
  TEST_EQUAL(u(8, 1, 1), 2.0);

  double norm = frobeniusNormOfSummand(t, 0);
  TEST_EQUAL(norm, 1.0);
  norm = frobeniusNormOfSummand(u, 1);
  TEST_EQUAL(norm, 0.0);

  norm = frobeniusNorm(t);
  TEST_EQUAL(norm, std::sqrt(3.0));

  double alpha = innerProduct(t, u);
  TEST_EQUAL(alpha, 12.0);

  CPTensor<double> m(t);
  TEST(equals(t, m));
  TEST(equals(m, t));
  m.setNull();
  TEST(!equals(m, t));

  // TODO: Add more tests...

  std::cout << "Tests passed!" << std::endl;

  PrintTensor pt(t);
  test_tensor(t, pt);

  return 0;
}
