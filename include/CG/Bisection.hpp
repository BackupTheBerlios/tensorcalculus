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

#ifndef __BISECT_HPP
#define __BISECT_HPP

#include <Utilities/Utilities.hpp>
#include <Constants.hpp>
#include <algorithm> // std::max, std::swap
#include <cmath> // std::abs

namespace TensorCalculus {

  template<typename T, typename F>
  class Bisection {

    private:
      T a;
      T b;
      T x;
      T residuum_;

      F func;

      unsigned int n;

    public:
      Bisection(T a, T b, F func)
        : a(a), b(b), func(func), n(0)
      {
        T a_value = func(a);
        T b_value = func(b);
        const bool a_pos = a_value >= 0;
        const bool b_pos = b_value >= 0;

        if (a_pos == b_pos)
          throw intervall_error("func(a) and func(b) "
                                "must have different signs");

        if (!b_pos)
        {
          std::swap(this->a, this->b);
          std::swap(a_value, b_value);
        }

        if (std::abs(a_value) < std::abs(b_value))
        {
          x = a;
          residuum_ = std::abs(a_value);
        } else {
          x = b;
          residuum_ = std::abs(b_value);
        }
      }

    void iterate(T epsilon = TensorCalculus::Constants<T>::epsilon) {
      if (residuum_ > epsilon)
      {
        x = 0.5 * (a + b);
        const T mid_value = func(x);

        if (mid_value <= 0) {
          a = x;
          residuum_ = std::abs(mid_value);
        } else if (mid_value > 0) {
          b = x;
          residuum_ = std::abs(mid_value);
        }
        ++n;
      }
    }

    const T residuum() const { return residuum_; }
    const T value() const { return x; }
    const unsigned int iterations() const { return n; }

    const T get_a() const { return a; }
    const T get_b() const { return b; }

  };

} // namespace TensorCalculus

#endif // _BISECT_HPP
