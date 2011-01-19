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

#ifndef __THREEPG_HPP
#define __THREEPG_HPP

#include <Utilities/Utilities.hpp>
// #include <algorithm> // std::max, std::swap
#include <cmath> // std::abs
#include "Function.hpp"

namespace TensorCalculus {

  template<typename T>
  class ThreePG {

  private:
    FunctionInterface<T, T> *func;
    T a;
    T b;
    T c;
    T R;

    bool c_sign;

    T fa;
    T fb;
    T fc;

    T C;
    T D;
    T residuum_;
    T epsilon;

    unsigned int n;

    void bisection_step(T& y, T& fy, bool& y_sign)
    {
      y = 0.5*(c + R);
      fy = (*func)(y);
      y_sign = fy >= 0;
      if (y_sign != c_sign) R = c;
    }

  public:

    ThreePG(T C, T D)
      : func(new FunctionInterface<T, T>), C(C), D(D), n(0)
    {
      if (C <= 0.5 || C >= 1.0 || D <= 0.5 || D >= 1.0)
        throw intervall_error("C, D must be out of (1/2, 1)");
    }

    template<typename F>
    ThreePG(T a, T b, T C, T D, T epsilon, F func)
      : C(C), D(D), n(0)
    {
      initialize(a, b, epsilon, func);
    }

    ~ThreePG()
    {
      delete func;
    }

    template<typename F>
    void initialize(T a, T b, T epsilon, F func)
    {
      if (a >= b)
        throw intervall_error("a < b must hold");
      if (epsilon <= 0)
        throw std::invalid_argument("Epsilon must be greater than zero");

      this->a = a;
      this->b = b;
      this->epsilon = epsilon;
      this->func = new Function<F, T, T>(func);
      n = 0;

      fa = (*func)(a);
      fb = (*func)(b);
      if (std::abs(fa) < epsilon)
      {
        c = a;
        residuum_ = std::abs(fa);
      } else if (std::abs(fb) < epsilon)
      {
        c = b;
        residuum_ = std::abs(fb);
      } else {
        const bool a_sign = fa >= 0;
        const bool b_sign = fb >= 0;

        if (a_sign == b_sign)
          throw intervall_error("func(a) and func(b) "
                                "must have different signs");

        c = b - fb*(b - a) / (fb - fa);
        fc = (*func)(c);
        residuum_ = std::abs(fc);
        c_sign = fc >= 0;

        if (c_sign != b_sign)
        {
          R = b;
        } else {
          R = a;
        }
      }
    }

    bool iterate()
    {
      if (residuum_ > epsilon)
      {
        const T Q = ( (c - a)*(c - a)*(fc - fb) - (c - b)*(c - b)*(fc - fa) ) /
          ( (b - a)*(c - a)*(c - b) );
        T y, fy;
        bool y_sign;
        if (std::abs(Q) < epsilon) {
          bisection_step(y, fy, y_sign);
        } else {
          y = c - fc / Q;
          fy = (*func)(y);
          if (!((y < R && y > c) || (y < c && y > R)))
          {
            // y not in (c,R) or (R,c)
            bisection_step(y, fy, y_sign);
          } else {
            y_sign = fy >= 0;
            T new_R;
            if (y_sign != c_sign)
            {
              new_R = c;
            } else {
              new_R = R;
            }
            if (std::abs(y - new_R) >= C*std::abs(R - c) ||
                std::abs(fy) >= D*std::abs(fc))
            {
              // Intervall (y, new_R) isn't significant shorter
              // than (c, R) or the f(y) significant better than
              // f(c)
              bisection_step(y, fy, y_sign);
            } else R = new_R;
          }
        }
        a = b;
        fa = fb;
        b = c;
        fb = fc;
        c = y;
        fc = fy;
        c_sign = y_sign;
        residuum_ = std::abs(fc);
        ++n;
        return false;
      } else return true;
    }

    const T value() const { return c; }
    const T residuum() const { return residuum_; }
    const unsigned int iterations() const { return n; }

    const T get_a() const { return a; }
    const T get_b() const { return b; }
    const T get_c() const { return c; }
    const T get_R() const { return R; }

  };

  template<typename T>
  ThreePG<T> threePG(T C, T D)
  {
    return ThreePG<T>(C, D);
  }

  template<typename T, typename F>
  ThreePG<T> threePG(T a, T b, T C, T D, T epsilon, F func)
  {
    return ThreePG<T>(a, b, C, D, epsilon, func);
  }

}

#endif // __THREEPG_HPP
