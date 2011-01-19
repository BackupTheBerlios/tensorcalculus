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

#ifndef __FACTORIALS_HPP
#define __FACTORIALS_HPP

#include <vector>
#include <algorithm>

namespace TensorCalculus {

  template<typename T>
  T factorial(unsigned int n) {
    T result = 1;
    for (unsigned int i = 2; i <= n; ++i) {
      result *= i;
    }
    return result;
  }

  template<typename T>
  T choose(unsigned int n, unsigned int k) {
    if (k > n) return 0;
    if (k > n/2) k = n-k;

    T result = 1;
    for (unsigned i = 0; i < k; ++i) {
      result *= static_cast<T>(n - i) / (k - i);
    }
    return result;
  }

  std::vector< std::vector<int> > partition(int sum, int num_parts) {
    std::vector< std::vector<int> > result;
    result.reserve(choose<int>(num_parts + sum - 1, sum));

    int t = sum;
    int h = 0;
    std::vector<int> next(num_parts, 0);
    next[0] = sum;
    result.push_back(next);
    while (next[num_parts - 1] != sum) {
      if (t > 1) h = 0;
      t = next[h];
      next[h] = 0;
      ++h;
      next[0] = t - 1;
      ++next[h];
      result.push_back(next);
    }

    return result;
  }

  bool unordered(const std::vector<int>& v) {
    for (unsigned int i = 0; i < v.size() - 1; ++i) {
      if (v[i] < v[i + 1]) return true;
    }
    return false;
  }

  std::vector< std::vector<int> > orderedPartition(const int sum, const int parts) {
    std::vector< std::vector<int> > result = partition(sum, parts);
    result.erase(std::remove_if(result.begin(), result.end(), unordered), result.end());
    return result;
  }


} // namespace TensorCalculus

#endif // __FACTORIALS_HPP
