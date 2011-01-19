/*
 * Copyright (C) 2010 Philipp Wähnert
 *               2011 Stefan Handschuh
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
 * \file Utilities.hpp
 * \brief Some generally usefull templates
 * \author Philipp Wähnert
 */

#ifndef __UTILITIES_HPP
#define __UTILITIES_HPP

#include <algorithm> // T std::max<T>(T, T), std::generate
// #include <cmath>     // std::fabs, std::max
#include <stdexcept> // std::logic_error
#include <complex>

namespace TensorCalculus {
  template<typename T>
  struct Utilities;

  template<typename InputIterator, typename Random>
  void fill_randomly(InputIterator begin, InputIterator end, Random rand)
  {
    std::generate(begin, end, rand);
  }

  /// \brief Removes the \c const before in the type \c T
  ///
  /// Actually in general nothing will be removed
  /// \see remove_const<const T>
  template<typename T>
  struct remove_const
  {
    typedef T type; // nothing to remove
  };

  /// \brief Removes the \c const before in the type \c T
  ///
  /// Partial template specialization for types containing a \c const
  /// \see remove_const
  template<typename T>           // partial template specialization of
  struct remove_const<const T>   // remove_const for types like const T
  {
    typedef T type; // removes the const
  };

  template<typename T>
  struct remove_ref
  {
    typedef T type;
  };
  
  template<typename T>
  struct remove_ref<T&>
  {
    typedef T type;
  };

  /// \brief Class to avoid memory leaks by using arrays
  ///        caused by thrown exceptions (RAII concept!)
  ///
  /// \remark The scoped_array holds the resource exclusivly
  ///         thus it isn't standard nor copy constructable
  ///         nor copyable
  template<typename T>
  class scoped_array
  {
  private:
    T* ptr;
    scoped_array();
    scoped_array(const scoped_array& x);
    scoped_array& operator= (const scoped_array& x);
  public:
    explicit scoped_array(T* p) // explicit!
      : ptr(p)
    { }

    ~scoped_array()
    {
      delete[] ptr; // Attention: Usage of delete[]
    }

    T& operator[] (int i) const { return ptr[i]; }
    T* get() const { return ptr; }
  };

  /// \brief Class to avoid memory leaks by using arrays
  ///        caused by thrown exceptions (RAII concept!)
  ///
  /// \remark The scoped_array holds the resource exclusivly
  ///         thus it isn't standard nor copy constructable
  ///         nor copyable
  template<typename T>
  class scoped_ptr
  {
  private:
    T* ptr;

    scoped_ptr(const scoped_ptr& x);
    scoped_ptr& operator= (const scoped_ptr& x);
  public:
    explicit scoped_ptr(T* p = 0) // explicit!
      : ptr(p)
    { }
    ~scoped_ptr()
    {
      delete ptr; // Attention: Usage of delete only for non-arrays
    }
    
    void reset(T* p)
    {
      delete ptr;
      ptr = p;
    }

    T& operator[] (int i) const { return ptr[i]; }
    T* get() const { return ptr; }
    
    T& operator * () const { return *ptr; }
    T* operator -> () const { return ptr; }
  };

  /// \brief Commonly used class to report an error while interpreting
  ///        stream data
  class StreamFormatError : public std::logic_error
  {
  public:
    StreamFormatError(const std::string& what) : std::logic_error(what) { }

    StreamFormatError()
      : std::logic_error("Stream has wrong format")
    { }
  };

  /// \brief Exception class reporting intervall errors in several iterations
  class intervall_error : public std::logic_error {
  public:
    intervall_error(const std::string& what)
      : std::logic_error(what) { }
  };

  template<typename T, typename F, typename I>
  T solve(T a, T b, T epsilon, F func, I& i)
  {
    i.initialize(a, b, epsilon, func);
    while(!i.iterate()) { }
    return i.value();
  }

  /*
    has been defined elsewhere (in some other header) and since this method is not
    used anywhere, I encommented it

    Stefan Handschuh, 08.02.2011 11:20
  void trim_left(std::string& string)
  {
    string.erase(0, string.find_first_not_of(' '));
  }
  */
  
  template<typename F>
  struct function_type;
  
  template<typename Res, typename Arg>
  struct function_type<Res (*)(Arg)> {
    typedef Res result_type;
    typedef Arg argument_type;
  };
  
  template<>
  struct Utilities<double> {
	static double rand(){
	  return (double) std::rand() / (double) RAND_MAX;
	}

	static double zero(){
	  return 0.0;
	}

	static double one(){
	  return 1.0;
	}
  };

  template<>
  struct Utilities< std::complex<double> > {
	static std::complex<double> rand() {
	  return std::complex<double>(Utilities<double>::rand(), Utilities<double>::rand());
	}

	static std::complex<double> zero() {
	  return std::complex<double>(0.0, 0.0);
	}
  };


} // namespace TensorCalculus

#endif // __UTILITIES_HPP
