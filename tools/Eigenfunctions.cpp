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

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

// #include <LinAlg/Vector.hpp>
// #include <LinAlg/Matrix.hpp>
// #include <LinAlg/Factorizations/Cholesky.hpp>
// #include <LinAlg/Factorizations/LU.hpp>
// #include <LinAlg/Eigenvalues/Symmetric.hpp>

#include "Mesh.hpp"
#include "BlasInterface.hpp"
#include "LapackInterface2.hpp"
// #include "ScopedPtr.hpp"

using namespace TensorCalculus;

// extern "C" void dsygvx_(const int* itype, const char* jobz, const char* range, const char* uplo,
//                     const int* n, double* a, const int* lda, double* b, const int* ldb,
//                     const double* vl, const double* vu,
//                     const int* il, const int* iu,
//                     const double* abstol,
//                     int* m, double* w, double* z, const int* ldz,
//                     double* work, const int* lwork, int* iwork,
//                     int* ifail, int* info);

const double ALPHA = 0.25;
// const int NUM_EIGENVALUES = 3;

double covariance(const std::vector<double>& x, const std::vector<double>& y) {
  // double arg = 1.0;
  // for (int i = 0; i < x.size(); ++i) {
  //   arg *= std::min(x[i],y[i]);
  // }
  // return arg;
  
  double arg = 0.0;
  for (int i = 0; i < x.size(); ++i) {
    arg += (x[i]-y[i])*(x[i]-y[i]);
  }
  return std::exp(-ALPHA * arg);

  // double arg = 0.0;
  // for (int i = 0; i < x.size(); ++i) {
  //   arg += (x[i]-y[i])*(x[i]-y[i]);
  // }
  // arg = 1 - std::sqrt(arg);
  // return arg > 0 ? arg : 0;  
}

// const std::string filename_base("data/rect");

void print_matrix(std::ostream& stream, const double *ptr, int height, int width)
{
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      stream << ptr[i + j*height] << ' ';
    }
    stream << '\n';
  }
}

int main(int argc, const char* argv[]) {
  if (argc < 3) {
    std::cout << "Usage: Eigenfunctions <meshfile> <num-evals>" << std::endl;
    return 1;
  }
  
  const std::string filename_base(argv[1]);
  const std::string num_evals(argv[2]);
  
  int NUM_EIGENVALUES;
  {
    std::istringstream iss(num_evals);
    iss >> NUM_EIGENVALUES;
    if (!iss) {
      std::cout << "Can't read number of eigenvalues" << std::endl;
      return 1;
    }
  }
  
  const std::string mesh_filename(filename_base + ".txt");
  std::ifstream mesh_file(mesh_filename.c_str());
  if (!mesh_file) {
    std::cout << "Can't open \"" << mesh_filename << '"' << std::endl;
    return 1;
  }
  
  std::cout << "Loading mesh...";
  std::cout.flush();
  Mesh<double> mesh(mesh_file);
  std::cout << "OK" << std::endl;
  
  const unsigned int POINTS = mesh.countPoints();
  const int dimension = mesh.getDimension();
  
  std::cout << "Creating mass matrix...";
  std::cout.flush();
  std::vector<double> M(POINTS*POINTS, 0.0);
  for (unsigned int i = 0; i < POINTS; ++i) {
    for (unsigned int j = 0; j < POINTS; ++j) {
      const double val = mesh.get_mass_matrix(i, j);
      M[i + j*POINTS] = val;
    }
  }
  // {
  //   std::ofstream mass_matrix("mass.mat");
  //   print_matrix(mass_matrix, &M[0], POINTS, POINTS);
  // }
  std::cout << "OK" << std::endl;

  std::cout << "Create covariance matrix...";
  std::cout.flush();
  std::vector<double> C(POINTS * POINTS);
  for (unsigned int i = 0; i < POINTS; ++i) {
    for (unsigned int j = 0; j < POINTS; ++j) {
      const double val = covariance(mesh.getPoint(i).getCoordinates(), mesh.getPoint(j).getCoordinates());
      C[i + j*POINTS] = val;
    }
  }
  // {
  //   std::ofstream covariance_matrix("covariance.mat");
  //   print_matrix(covariance_matrix, &C[0], POINTS, POINTS);
  // }
  std::cout << "OK" << std::endl;

  std::cout << "Calculate eigenvalues and -vectors...";
  std::cout.flush();
  std::vector<double> eval(POINTS);
  std::vector<double> evec(POINTS * NUM_EIGENVALUES);
  {
    int m;
    double dbl_lwork;
    std::vector<int> iwork(5 * POINTS);
    std::vector<int> ifail(POINTS);
    int result =
      Lapack<double>::sygvx(2, 'V', 'I', 'U',
                            POINTS, &C[0], POINTS, &M[0], POINTS,
                            0.0, 0.0, POINTS - NUM_EIGENVALUES + 1, POINTS, -1.0,
                            m, &eval[0], &evec[0], POINTS,
                            &dbl_lwork, -1, &iwork[0], &ifail[0]);
    if (result != 0) {
      std::cout << "Query failed! Error: " << result << std::endl;
      return 1;
    }
    int lwork = static_cast<int>(dbl_lwork);
    std::vector<double> work(lwork);
    result =
      Lapack<double>::sygvx(2, 'V', 'I', 'L',
                            POINTS, &C[0], POINTS, &M[0], POINTS,
                            0.0, 0.0, POINTS - NUM_EIGENVALUES + 1, POINTS, -1.0,
                            m, &eval[0], &evec[0], POINTS,
                            &work[0], lwork, &iwork[0], &ifail[0]);
    if (result != 0) {
      std::cout << "Failed! Error: " << result << std::endl;
      return 1;
    }
  }
  std::cout << "OK" << std::endl;
  
  // std::cout << "Scaling eigenvalues...";
  // std::cout.flush();
  // {
  //   std::vector<double> temp(POINTS, 0.0);
  //   for (int i = 0; i < NUM_EIGENVALUES; ++i) {
  //       Blas<double>::gemv('N', POINTS, POINTS, 1.0, &M[0], POINTS, &evec[(NUM_EIGENVALUES - i -1) * POINTS], 1, 0.0, &temp[0], 1);
  //       const double norm = Blas<double>::dot(POINTS, &evec[(NUM_EIGENVALUES - i -1) * POINTS], 1, &temp[0], 1);
  //       Blas<double>::scal(POINTS, 1.0/std::sqrt(norm), &evec[(NUM_EIGENVALUES - i -1) * POINTS], 1);
  //   }
  // }
  // The factor obtained here isn't our right one...
  // for (int i = 0; i < NUM_EIGENVALUES; ++i) {
  //   const double norm = Blas<double>::nrm2(POINTS, &evec[i * POINTS], 1);
  //   Blas<double>::scal(POINTS, 1.0/norm, &evec[i * POINTS], 1);
  // }
  // std::cout << "OK" << std::endl;

  std::cout << "Performing output operations...";
  std::cout.flush();
  {
    std::ofstream output((filename_base + ".evl").c_str());
    // output << NUM_EIGENVALUES << std::endl;
    for (int i = 0; i < NUM_EIGENVALUES; ++i) {
      output << eval[NUM_EIGENVALUES - i - 1] << '\n';
      // output << gsl_vector_get(eval.get(), NUM_EIGENVALUES - i - 1) << std::endl;
    }
  }
  {
    std::ofstream output((filename_base + ".evc").c_str());
    // output << POINTS << " " << NUM_EIGENVALUES << std::endl;
    for (unsigned int point_index = 0; point_index < POINTS; ++point_index) {
      // for (unsigned int i = 0; i < dimension; ++i) {
      //   output << mesh.getPoint(point_index)[i] << " ";
      // }
      // output << point_index << " ";
      for (unsigned int j = 0; j < NUM_EIGENVALUES; ++j) {
        output << evec[point_index + POINTS*(NUM_EIGENVALUES - j - 1)] << " ";
        // output << gsl_matrix_get(E.get(), point_index, NUM_EIGENVALUES - j - 1) << " ";
      }
      output << '\n';
    }
  }
  {
    for (int i = 0; i < NUM_EIGENVALUES; ++i) {
      Blas<double>::scal(POINTS, std::sqrt(eval[i]), &evec[i * POINTS], 1);
    }
    std::ofstream output((filename_base + ".func").c_str());
    output << POINTS << ' ' << NUM_EIGENVALUES << '\n';
    for (int i = 0; i < POINTS; ++i) {
      for (int j = 0; j < NUM_EIGENVALUES; ++j) {
        output << evec[i + (NUM_EIGENVALUES - j - 1)*POINTS] << " ";
      }
      output << '\n';
    }
  }
  std::cout << "OK" << std::endl;

  return 0;
}
