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
 * \file WeakVector.hpp
 * \brief Declares the WeakVector class
 * \author Philipp Wähnert
 */

#ifndef __WEAK_VECTOR_HPP
#define __WEAK_VECTOR_HPP

#include <cstddef>            // std::size_t, std::ptrdiff_t
#include <algorithm>          // std::copy
#include <iterator>           // std::reverse_iterator<T>, std::distance, ...
#include <vector>             // std::vector<T>
#include <limits>             // std::numeric_limits<T>
#include "Utilities/Utilities.hpp"      // remove_const<T>
#include <stdexcept>          // std::out_of_range
#include "BlasInterface.hpp"  // Blas<T>

// // #define RANGE_CHECKS_ON

/*!
 * \namespace TensorCalculus
 * \brief Namespace containing the TensorCalculus library
 */
namespace TensorCalculus {

  /*!
   * \class WeakVector
   * \brief Class template providing a similar interface like
   *        the \c std::vector<T>
   *
   * The WeakVector behaves like a reference. At construction time it
   * must be set to a special location given by a pointer and a length
   * of the range. It is copy but not standard constructable. After
   * construction its internal pointers can't be changed anymore and
   * every operation is exclusivly carried out on the data it points to.
   *
   * Constness of \c T means the constness of the data it points to. Then
   * the non-\c const and \c const versions of the public types - like
   * \c iterator and \c const_iterator - are the same. Assignment won't
   * work because this operation needs a non-\c const pointer \c _first.
   *
   * In constrast to that constness of a WeakVector itself means \c size
   * can't be changed
   *
   * It provides a similar interface like the \c std::vector<T> container
   * and thus it is ready-for-use for the standard algorithms
   *
   * \tparam T data type of its elements
   */
  template<typename T>
  class WeakVector
  {
  public:

    typedef T value_type;
    typedef T* pointer;
    typedef T& reference;
    typedef const T& const_reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef T* iterator;
    typedef const T* const_iterator;
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    /// \param ptr Pointer to the memory containing the data
    /// \param size Size of the data memory
    WeakVector(pointer ptr, int size = 0)
      : _first(ptr), _end(ptr+size), _size(size)
    { }

    /// \brief Copy constructor
    template<typename S>
    WeakVector(const WeakVector<S>& weakvector)
      : _first(&weakvector.front()), _end(&weakvector.back() + 1),
        _size(weakvector.size())
    { }

    /// \param vector Vector providing the data memory
    WeakVector(std::vector<typename remove_const<T>::type>& vector)
      : _first(&vector.front()), _end(&vector.back() + 1),
        _size(vector.size())
    { }

    /// \brief Assignment operator
    /// \attention Assigns the contents and doesn't change the
    ///            internal pointers!
    /// \remark Doesn't work if \c T is \c const !
    template<typename S>
    const WeakVector& operator= (const WeakVector<S>& y)
    {
      const typename WeakVector::size_type min_size =
        std::min(size(), y.size());

      Blas<T>::copy(min_size, &y[0], 1, _first, 1);
      return *this;
    }

    /// \brief Assignment operator
    /// \attention Assigns the contents and doesn't change the
    ///            internal pointers!
    /// \remark Doesn't work if \c T is \c const !
    const WeakVector&
    operator= (const std::vector<typename remove_const<T>::type>& y)
    {
      const typename WeakVector<T>::size_type min_size =
        std::min(size(), y.size());
      Blas<T>::copy(min_size, &y[0], 1, _first, 1);
      return *this;
    }

    // iterator begin() { return _first; }
    iterator begin() const { return _first; }
    // const_iterator begin() const { return _first; }
    // iterator end() { return _end; }
    iterator end() const { return _end; }
    // const_iterator end() const { return _end; }

    // reverse_iterator rbegin() { return _end - 1; }
    reverse_iterator rbegin() const { return _end - 1; }
    // const_reverse_iterator rbegin() const { return _end - 1; }
    // reverse_iterator rend() { return _first - 1; }
    reverse_iterator rend() const { return _first - 1; }
    // const_reverse_iterator rend() const { return _first - 1; }

    size_type size() const { return _size; }
    size_type max_size() const { return std::numeric_limits<size_type>::max(); }

    void resize(size_type s)
    {
      _size = s;
      _end = _first + _size;
    }

    bool empty() const { return _size == 0; }

    // reference at(size_type n) {
    //   range_check(n);
    //   return _first[n];
    // }

    reference at(size_type n) const
    {
      range_check(n);
      return _first[n];
    }

    // const_reference at(size_type n) const {
    //   range_check(n);
    //   return _first[n];
    // }

    reference operator[] (size_type n) const { return _first[n]; }
    // const_reference operator[] (size_type n) const { return _first[n]; }
    // reference operator[] (size_type n) { return _first[n]; }

    // reference front() { return *_first; }
    reference front() const { return *_first; }
    // const_reference front() const { return *_first; }
    // reference back() { return *(_end - 1); }
    reference back() const { return *(_end - 1); }
    // const_reference back() const { return *(_end - 1); }

    template<typename InputIterator>
    void assign(InputIterator a, InputIterator b)
    {
      _size = std::distance(a, b);
      _end = _first + _size;
      std::copy(a, b, _first);
    }

    template<typename InputIterator>
    void assign(InputIterator a)
    {
      InputIterator b = a;
      std::advance(b, _size);
      std::copy(a, b, _first);
    }

    void assign(size_type n, T x = T())
    {
      _size = n;
      _end = _first + _size;
      std::fill_n(_first, n, x);
    }

    void push_back(const T& x)
    {
      ++_size;
      *_end = x;
      ++_end;
    }

    void clear()
    {
      _size = 0;
      _end = _first;
    }

    operator WeakVector<const T> ()
    {
      return WeakVector<const T>(*this);
    }

  private:

    pointer _first;
    pointer _end;
    size_type _size;

    void range_check(size_type n) const
    {
      if ((n < 0) || (n >= _size))
      {
        throw std::out_of_range("WeakVector::range_check");
      }
    }

  };

} // namespace TensorCalculus

#endif // __WEAK_VECTOR_HPP
