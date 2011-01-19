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
 * \file Function.hpp
 * \brief General function functors
 * \author Philipp Wähnert
 */

#ifndef __FUNCTION_HPP
#define __FUNCTION_HPP

#include <stdexcept>
#include <functional>

namespace TensorCalculus {

  template<typename Arg, typename Res>
  struct FunctionInterface : public std::unary_function<Arg, Res> {

    virtual Res operator() (Arg a)
    {
      throw std::runtime_error("Invokation of abstract "
                               "'FunctionInterface' functor");
    }

  };

  template<typename Param, typename Arg, typename Res>
  struct ParameterFunctionInterface
    : public std::binary_function<Param, Arg, Res> {

    virtual Res operator() (Param p, Arg a)
    {
      throw std::runtime_error("Invokation of abstract "
                               "'ParameterFunctionInterface' functor");
    }

  };

  template<typename F, typename Arg, typename Res>
  struct Function : public FunctionInterface<Arg, Res> {

  private:
    F f;
  public:
    Function(F f) : f(f) { }

    virtual Res operator() (Arg a)
    {
      return f(a);
    }

  };

  template<typename F, typename Param, typename Arg, typename Res>
  struct ParameterFunction
    : public ParameterFunctionInterface<Param, Arg, Res> {

  private:
    F f;
  public:
    ParameterFunction(F f) : f(f) { }

    virtual Res operator() (Param p, Arg a)
    {
      return f(p, a);
    }

  };

} // namespace TensorCalculus

#endif // __FUNCTION_HPP
