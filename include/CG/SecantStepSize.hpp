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

#ifndef __SECANTSTEPSIZE_HPP
#define __SECANTSTEPSIZE_HPP

#include <cmath>

namespace TensorCalculus {

  template<typename InnerProductSpace, template<typename> class InnerProductSpaceTraits = StandardInnerProductSpaceTraits>
  class SecantStepSize {
  public:
    SecantStepSize(typename InnerProductSpaceTraits<InnerProductSpace>::Scalars step,
                   typename InnerProductSpaceTraits<InnerProductSpace>::Scalars epsilon,
                   InnerProductSpaceTraits<InnerProductSpace> traits = InnerProductSpaceTraits<InnerProductSpace>())
      : step(step), epsilon(epsilon), traits(traits) { }
      
    template<typename Function>
    typename InnerProductSpaceTraits<InnerProductSpace>::Scalars
    operator () (Function df,
                 const typename InnerProductSpaceTraits<InnerProductSpace>::Vectors& df_x,
                 const typename InnerProductSpaceTraits<InnerProductSpace>::Vectors& x,
                 const typename InnerProductSpaceTraits<InnerProductSpace>::Vectors& d)
    {
      using std::abs;
      typename InnerProductSpaceTraits<InnerProductSpace>::Scalars a = 0;
      typename InnerProductSpaceTraits<InnerProductSpace>::Scalars phi_a = traits.innerProduct(df_x, d);
      
      typename InnerProductSpaceTraits<InnerProductSpace>::Scalars b = step;
      typename InnerProductSpaceTraits<InnerProductSpace>::Vectors y = x;
      traits.update(b, d, y);
      typename InnerProductSpaceTraits<InnerProductSpace>::Scalars phi_b = traits.innerProduct(df(y), d);
      if (abs(phi_b) < epsilon) {
        return b;
      }

      typename InnerProductSpaceTraits<InnerProductSpace>::Vectors z = x;
      
      const int max_n = 100;
      int n = 0;
      while (phi_b < 0 && n++ < max_n) {
        b += step;
        traits.update(step, d, y);
        phi_b = traits.innerProduct(df(y), d);
        if (abs(phi_b) < epsilon) {
          return b;
        }
      }
      if (n >= max_n) {
        std::cout << '*' << std::endl;
        return b;
      }
      
      typename InnerProductSpaceTraits<InnerProductSpace>::Scalars c
        = a - phi_a * (b - a) / (phi_b - phi_a);
      z = x;
      traits.update(c, d, z);
      typename InnerProductSpaceTraits<InnerProductSpace>::Scalars phi_c
        = traits.innerProduct(df(z), d);

      n = 0;
      while (abs(phi_c) >= epsilon && n < max_n) {
        if (phi_c > 0) {
          b = c;
          phi_b = phi_c;
        } else {
          a = c;
          phi_a = phi_c;
        }
        if (abs(phi_c) < epsilon) {
          return c;
        }
        c = a - phi_a * (b - a) / (phi_b - phi_a);
        z = x;
        traits.update(c, d, z);
        phi_c = traits.innerProduct(df(z), d);
        ++n;
      }
      // if (n >= max_n) {
      //   return 0;
      // }

      return c;
    }
  
  private:
    typename InnerProductSpaceTraits<InnerProductSpace>::Scalars step;
    typename InnerProductSpaceTraits<InnerProductSpace>::Scalars epsilon;
    InnerProductSpaceTraits<InnerProductSpace> traits;
  };

} // namespace TensorCalculus

#endif // __SECANTSTEPSIZE_HPP
