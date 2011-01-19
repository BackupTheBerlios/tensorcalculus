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

#ifndef __STANDARDTAITS_HPP
#define __STANDARDTAITS_HPP

// These standard templates can be specialized appropriately

namespace TensorCalculus {

  template<typename T>
  struct StandardVectorSpaceTraits {
    typedef T Vectors;
    typedef T Scalars;
    
    void update(const Scalars alpha, const Vectors& x, Vectors& y) const
    {
      y += alpha * x;
    }
    
    void scale(const Scalars alpha, Vectors& x) const
    {
      x *= alpha;
    }
    
    Vectors zero() const {
      return Vectors();
    }
    
  };
  
  template<typename T>
  struct StandardInnerProductSpaceTraits : public StandardVectorSpaceTraits<T> {
    using StandardVectorSpaceTraits<T>::update;
    using StandardVectorSpaceTraits<T>::scale;
    using StandardVectorSpaceTraits<T>::zero;
    
    typedef T Vectors;
    typedef T Scalars;
    
    Scalars innerProduct(const Vectors& x, const Vectors& y) const
    {
      return x * y;
    }
  };

  template<typename Monoid>
  struct StandardMonoidTraits {
    typedef Monoid MonoidType;
    
    MonoidType identityElement() const { return MonoidType(1); }
    void update(const MonoidType& x, MonoidType& y) const { y *= x; }
  };

} // namespace TensorCalculus

#endif // __STANDARDTAITS_HPP
