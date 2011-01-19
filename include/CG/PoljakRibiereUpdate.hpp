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

#ifndef __POLJAKRIBIEREUPDATE_HPP
#define __POLJAKRIBIEREUPDATE_HPP

#include <algorithm>

#include "StandardTraits.hpp"

namespace TensorCalculus {

  template<typename InnerProductSpace, template<typename> class InnerProductSpaceTraits = StandardInnerProductSpaceTraits>
  class PoljakRibiereUpdate {
  public:
    explicit PoljakRibiereUpdate(bool restart, InnerProductSpaceTraits<InnerProductSpace> traits = InnerProductSpaceTraits<InnerProductSpace>())
      : restart(restart), traits(traits) { }
    
    typename InnerProductSpaceTraits<InnerProductSpace>::Scalars
    operator () (const typename InnerProductSpaceTraits<InnerProductSpace>::Vectors& g_new,
                 const typename InnerProductSpaceTraits<InnerProductSpace>::Vectors& g) const
    {
      typename InnerProductSpaceTraits<InnerProductSpace>::Vectors diff_g = g_new;
      traits.update(-1, g, diff_g);
      // diff_g -= g;
      const typename InnerProductSpaceTraits<InnerProductSpace>::Scalars beta
        = traits.innerProduct(diff_g, g_new) / traits.innerProduct(g, g);
      using std::max;
      return max(beta, static_cast<typename InnerProductSpaceTraits<InnerProductSpace>::Scalars>(0));
    }
  
  private:
    bool restart;
    InnerProductSpaceTraits<InnerProductSpace> traits;
  };

  //template<typename Scalars>
  //struct PoljakRibiereUpdate {
    
    //template<typename InnerProductSpace>
    //Scalars operator () (const InnerProductSpace& g_new, const InnerProductSpace& g) const
    //{
      //InnerProductSpace diff_g = g_new;
      //diff_g -= g;
      //return innerProduct(diff_g, g_new) / innerProduct(g, g);
    //}
    
  //};
  
  //template<typename Scalars>
  //PoljakRibiereUpdate<Scalars> poljak_ribiere_update()
  //{
    //return PoljakRibiereUpdate<Scalars>();
  //}

  // template<typename VectorSpace>
  // struct PoljakRibiereUpdate {
  // private:
  //   typedef typename VectorSpace::vectors V;
  //   typedef typename VectorSpace::scalars K;
  // 
  // public:
  //   static K updateFactor(V g_new, V g) {
  //     V diff_g = g_new;
  //     diff_g -= g;
  //     return scalarProduct(diff_g, g_new) / scalarProduct(g, g);
  //   }
  // 
  // };

} // namespace TensorCalculus

#endif // __POLJAKRIBIEREUPDATE_HPP
