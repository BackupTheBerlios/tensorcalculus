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
 * \file TensorBasisOperators.hpp
 * \brief Declares function to transform a CPTensor into an
 *        other CPTensor by a change of basis defined by a TensorBasis
 * \author Philipp Wähnert
 */

#ifndef __TENSOR_BASIS_OPERATORS_HPP
#define __TENSOR_BASIS_OPERATORS_HPP

#include <ostream>             // std::ostream
#include <algorithm>           // std::max_element
#include <iterator>            // std::ostream_iterator<T>
#include <iomanip>             // std::setprecision
#include "Tensor/CPTensor.hpp"        // CPTensor<T>
#include "Tensor/TensorBasis.hpp"     // TensorBasis<T>
#include "Constants.hpp"       // Constants<T>::epsilon

namespace TensorCalculus {

  /// \brief Generates
  template<typename T>
  void generate(const CPTensor<T>& x, const TensorBasis<T>& b,
                CPTensor<T>& y)
  {
    const int spat_dirs_x = x.getSpatialDirections();
    const int rang_x = x.getSeparationRang();
    const std::vector<int>& dims_x = x.getDimensions();

    const int spat_dirs_b = b.getSpatDirs();
    const std::vector<int>& dims_b = b.getDimsOfSubspaces();
    const std::vector<int>& bases_b = b.getDimsOfBases();

#ifdef RANGE_CHECKS_ON
    if ((spat_dirs_x != spat_dirs_b) ||
        (dims_x != bases_b))
    {
      throw std::out_of_range("generate(const CPTensor<T>&, "
                              "const TensorBasis<T>&, CPTensor<T>&)"
                              " : Dimensions mismatch");
    }
#endif

    y.resize(rang_x, spat_dirs_x, dims_b);
    for (int mu = 0; mu < spat_dirs_x; mu++) {
      Blas<T>::gemm(NoTrans, NoTrans, dims_b[mu], rang_x, dims_x[mu], 1,
                    &b.getBasis(mu)[0], dims_b[mu], &x(mu)[0], dims_x[mu],
                    0.0, &y(mu)[0], dims_b[mu]);
    }
  }

  template<typename T>
  bool isOrthogonal(const TensorBasis<T>& b,
                    T epsilon = Constants<T>::epsilon)
  {
    const int d = b.get_d();
    const std::vector<T>& n = b.get_n();
    const std::vector<T>& m = b.get_m();

    for (int mu = 0; mu < d; ++mu)
    {
      for (int k = 0; k < m[mu]; ++k)
      {
        const T* ptr = &b(mu, k)[0];
        for (int l = k+1; l < m[mu]; ++l)
        {
          const T result = Blas<T>::dot(n[mu], ptr, 1, &b(mu, l)[0], 1);
          if ( abs(result) >= epsilon) return false;
        }
      }
    }

    return true;
  }

  template<typename T>
  bool isOrthonormal(const TensorBasis<T>& b,
                     T epsilon = Constants<T>::epsilon)
  {
    const int d = b.get_d();
    const std::vector<T>& n = b.get_n();
    const std::vector<T>& m = b.get_m();

    for (int mu = 0; mu < d; ++mu)
    {
      for (int k = 0; k < m[mu]; ++k)
      {
        const T* ptr = &b(mu, k)[0];
        T result = nrm2(n[mu], ptr, 1);

        if (abs(result - 1) >= epsilon) return false;

        for (int l = k+1; l < m[mu]; ++l)
        {
          const T result = Blas<T>::dot(n[mu], ptr, 1, &b(mu, l)[0], 1);
          if ( abs(result) >= epsilon) return false;
        }
      }
    }
    return true;
  }

  template<typename T>
  void coefficients(const CPTensor<T>& y, const TensorBasis<T>& b,
                    CPTensor<T>& x)
  {
    const int d = y.get_d();
    const std::vector<int>& n = y.get_n();

#ifdef RANGE_CHECKS_ON
    if ((d != b.get_d()) || (n != b.get_n()))
    {
      throw std::out_of_range("coefficients(const CPTensor<T>&, "
                              "const TensorBasis<T>&, CPTensor<T>&)"
                              " : Dimensions mismatch");
    }
#endif

    const std::vector<int>& m = b.get_m();
    const int r = y.get_r();

    x.resize(r, d, m);
    for (int mu = 0; mu < d; ++mu)
    {
      // Blas<T>::gemm(true, false, true, r, m[mu], n[mu], 1.0,
      //               &y(mu)[0], n[mu],
      //               &b(mu)[0], n[mu], 0.0,
      //               &x(mu)[0], m[mu]);
    }
  }

  template<typename T>
  std::ostream& operator<< (std::ostream& out, const TensorBasis<T>& b)
  {
    const int d = b.getSpatDirs();
    out << "(" << d << ", [ ";

    const std::vector<int>& dims = b.getDimsOfSubspaces();
    std::copy(dims.begin(), dims.end(), std::ostream_iterator<int>(out, " "));
    out << "], [ ";

    const std::vector<int>& bases = b.getDimsOfBases();
    std::copy(bases.begin(), bases.end(), std::ostream_iterator<int>(out, " "));
    out << "])" << std::endl;

    out << std::setprecision(3) << std::scientific
        << std::setiosflags(std::ios::left) << std::setw(8);

    const int max_dim = *std::max_element(dims.begin(), dims.end());

    for (int i = 0; i < max_dim; i++)
    {
      for (int mu = 0; mu < d; mu++)
      {
        for (int k = 0; k < bases[mu]; k++)
        {
          if (i < dims[mu])
          {
            out << b.getBasisVector(mu, k)[i] << " ";
          } else {
            out << "          ";
          }
        }
        out << "|";
      }
      out << std::endl;
    }

    return out;
  }

} // namespace TensorCalculus

#endif // __TENSOR_BASIS_OPERATORS_HPP
