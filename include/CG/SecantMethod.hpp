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

#ifndef __SECANTMETHOD_HPP
#define __SECANTMETHOD_HPP

#include "StandardTraits.hpp"

namespace TensorCalculus {

  template<typename InnerProductSpace, template<typename> class InnerProductSpaceTraits = StandardInnerProductSpaceTraits>
  class SecantRule {
  public:
    explicit SecantRule(typename InnerProductSpaceTraits<InnerProductSpace>::Scalars sigma,
                        InnerProductSpaceTraits<InnerProductSpace> traits = InnerProductSpaceTraits<InnerProductSpace>())
      : sigma(sigma), traits(traits) { }
      
    template<typename Function>
    typename InnerProductSpaceTraits<InnerProductSpace>::Scalars
    operator () (Function df,
                 const typename InnerProductSpaceTraits<InnerProductSpace>::Vectors& df_x,
                 const typename InnerProductSpaceTraits<InnerProductSpace>::Vectors& x,
                 const typename InnerProductSpaceTraits<InnerProductSpace>::Vectors& d)
    {
      // const typename InnerProductSpaceTraits<InnerProductSpace>::Vectors df_x = df(x);
      const typename InnerProductSpaceTraits<InnerProductSpace>::Scalars df_xd = traits.innerProduct(df_x, d);
      typename InnerProductSpaceTraits<InnerProductSpace>::Vectors y = x;
      traits.update(sigma, d, y);
      typename InnerProductSpaceTraits<InnerProductSpace>::Vectors df_y = df(y);
      // df(y, df_y);
      const typename InnerProductSpaceTraits<InnerProductSpace>::Scalars df_yd = traits.innerProduct(df_y, d);
      // InnerProductSpace step = d;
      // step *= sigma;
      // y += step;
      
      return -sigma * df_xd / (df_yd - df_xd);
    }
    
    template<typename Function>
    typename InnerProductSpaceTraits<InnerProductSpace>::Scalars
    operator () (Function function,
                 const typename InnerProductSpaceTraits<InnerProductSpace>::Vectors& x,
                 const typename InnerProductSpaceTraits<InnerProductSpace>::Vectors& d)
    {
      return operator () (function, df(x), x, d);
    }
    
  private:
    typename InnerProductSpaceTraits<InnerProductSpace>::Scalars sigma;
    InnerProductSpaceTraits<InnerProductSpace> traits;
  };

  //template<typename Function, typename Scalars>
  //class SecantRule {
  //public:
    //SecantRule(Function df, Scalars sigma)
      //: df(df), sigma(sigma) { }

    //template<typename InnerProductSpace>
    //Scalars operator () (const InnerProductSpace& x, const InnerProductSpace& d) const
    //{
      //const InnerProductSpace df_x = df(x);
      //const Scalars scal_dfx_d = innerProduct(df_x, d);
      //InnerProductSpace y = x;
      //InnerProductSpace step = d;
      //step *= sigma;
      //y += step;
      //const InnerProductSpace df_y = df(y);

      //return -sigma*scal_dfx_d / (innerProduct(df_y, d) - scal_dfx_d);
    //}
    
  //private:
    //Function df;
    //Scalars sigma;
    
  //};
  
  //template<typename Function, typename Scalars>
  //SecantRule<Function, Scalars> secant_rule(Function df, Scalars sigma)
  //{
    //return SecantRule<Function, Scalars>(df, sigma);
  //}

  // template<typename VectorSpace>
  // class SecantMethod {
  // private:
  //   typedef typename VectorSpace::vectors V;
  //   typedef typename VectorSpace::scalars K;
  // 
  //   std::tr1::function<V (V)> df;
  //   K sigma;
  //   V x;
  //   V d;
  // 
  // 
  // public:
  //   SecantMethod(K sigma) : sigma(sigma) { }
  // 
  //   template<typename D>
  //   void initialize(D df) {
  //     this->df = df;
  //   }
  // 
  //   K stepSize(const V& x, const V& d) {
  //     const V df_x = df(x);
  //     const K scal_dfx_d = scalarProduct(df_x, d);
  //     V y = x;
  //     V step = d;
  //     step *= sigma;
  //     y += step;
  //     const V df_y = df(y);
  // 
  //     return -sigma*scal_dfx_d / (scalarProduct(df_y, d) - scal_dfx_d);
  //   }
  // };

} // namespace TensorCalculus

#endif // __SECANTMETHOD_HPP
