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

/*!
 * \file Functional.hpp
 * \brief Some usefull functors
 * \author Philipp Wähnert
 */

#ifndef __FUNCTIONAL_HPP
#define __FUNCTIONAL_HPP

#include <functional>
#include <cmath>

namespace TensorCalculus {

  template <typename F, typename G, typename H>
  class binary_compose
    : public std::unary_function<typename G::argument_type,
                                 typename F::result_type>
  {
  protected:
    F f;
    G g;
    H h;
  public:
    binary_compose(const F& __f, const G& __g, const H& __h)
      : f(__f), g(__g), h(__h)
    { }

    typename F::result_type
    operator()(const typename G::argument_type& x) const {
      return f(g(x), h(x));
    }
  };

  /// \brief Functor providing the absolute less comparison
  template<typename T>
  struct abs_less : public std::binary_function<T, T, bool>
  {
    bool operator() (T x, T y) {
      return std::abs(x) < std::abs(y);
    }
  };

} // namespace TensorCalculus

#endif // __FUNCTIONAL_HPP
