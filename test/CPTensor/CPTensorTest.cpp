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

#include <fstream>

#include "Tensor/CPTensor.hpp"
#include "TestUtilities.hpp"
#include "Constants.hpp"

namespace TensorCalculus {

  template<>
  struct Constants<double> {
    static const double epsilon;
  };

}

const double TensorCalculus::Constants<double>::epsilon = 1e-8;

using namespace TensorCalculus;

int main(int argc, const char* argv[]) {
  std::vector<int> n(2);
  n[0] = 3;
  n[1] = 2;
  CPTensor<double> t(3, 2, n);
  TEST_EQUAL(t.get_d(), 2);
  TEST_EQUAL(t.get_r(), 3);
  TEST_EQUAL(t.get_n()[0], 3);
  TEST_EQUAL(t.get_n()[1], 2);

  std::vector<double>& v = t.getVectorOfDimension(0);
  TEST_EQUAL(v.size(), 9);
  v[0] = 1.0;
  v[4] = 1.0;
  v[8] = 1.0;
  TEST_EQUAL(t.getVector(0, 0, 0), 1.0);
  TEST_EQUAL(t.getVector(1, 0, 1), 1.0);
  TEST_EQUAL(t.getVector(2, 0, 2), 1.0);

  std::vector<double>& w = t.getVectorOfDimension(1);
  TEST_EQUAL(w.size(), 6);
  w[0] = 1.0;
  w[2] = 1.0;
  w[5] = 1.0;
  TEST_EQUAL(t.getVector(0, 1, 0), 1.0);
  TEST_EQUAL(t.getVector(1, 1, 0), 1.0);
  TEST_EQUAL(t.getVector(2, 1, 1), 1.0);

  WeakVector<double> a = t.getVectorOfRepresentants(0, 0);
  TEST_EQUAL(a.size(), 3);
  TEST_EQUAL(a[0], 1);
  TEST_EQUAL(a[1], 0);
  TEST_EQUAL(a[2], 0);

  WeakVector<double> b = t.getVectorOfRepresentants(1, 1);
  TEST_EQUAL(b.size(), 2);
  TEST_EQUAL(b[0], 1);
  TEST_EQUAL(b[1], 0);

  std::vector<int> i(2);
  i[0] = 0; i[1] = 0;
  TEST_EQUAL(t.getTensorEntry(i), 1.0);
  i[0] = 0; i[1] = 1;
  TEST_EQUAL(t.getTensorEntry(i), 0.0);
  i[0] = 1; i[1] = 0;
  TEST_EQUAL(t.getTensorEntry(i), 1.0);
  i[0] = 1; i[1] = 1;
  TEST_EQUAL(t.getTensorEntry(i), 0.0);
  i[0] = 2; i[1] = 0;
  TEST_EQUAL(t.getTensorEntry(i), 0.0);
  i[0] = 2; i[1] = 1;
  TEST_EQUAL(t.getTensorEntry(i), 1.0);

  if (argc > 1) {
    {
      std::ofstream file(argv[1]);
      if (!file) {
        std::cerr << "Can't open file " << argv[1] << '\n';
        return -1;
      }
      t.write(file);
    }
    {
      std::ifstream file(argv[1]);
      if (!file) {
        std::cerr << "Can't open file " << argv[1] << '\n';
        return -1;
      }
      const CPTensor<double> s(file);
      TEST_EQUAL(t.get_d(), s.get_d());
      TEST_EQUAL(t.get_r(), s.get_r());
      for (int mu = 0; mu < s.get_d(); ++mu) {
        TEST_EQUAL(t.get_n(mu), s.get_n(mu));
        
        const std::vector<double>& t_mu = t.getVectorOfDimension(mu);
        const std::vector<double>& s_mu = s.getVectorOfDimension(mu);
        
        for (unsigned int i = 0; i < s.getVectorOfDimension(mu).size(); ++i) {
          TEST_EQUAL(t_mu[i], s_mu[i]);
        }
      }
    }
  }

  std::cout << "Tests passed!" << std::endl;
  return 0;
}
