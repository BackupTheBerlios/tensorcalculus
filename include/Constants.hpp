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
 * \file Constants.hpp
 * \brief Declares the struct containing the essential constants
 *        used by the TensorCalculus library
 * \author Philipp Wähnert
 */

#ifndef __CONSTANTS_HPP
#define __CONSTANTS_HPP

namespace TensorCalculus {

  /// \struct Constants
  /// \brief Template containing several constants used by the 
  ///        TensorCalculus library.
  /// \details These constants must be defined after including this
  ///   header but before using the library. For example like this:
  /// \dontinclude DoubleConstants.cpp
  /// \skip #include
  /// \until epsilon
  /// \remark By linking to DoubleConstants.cpp you already get these 
  ///         constants defined to the usual values
  template<typename T>
  struct Constants
  {
    /// \brief Values below epsilon are accepted as small enough.
    ///        For example used by <tt>CPTensor::equals</tt>.
    static const T epsilon;
  };

} // namespace TensorCalculus

#endif // __CONSTANTS_HPP
