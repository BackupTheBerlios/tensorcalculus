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
#include <vector>
#include <cstdlib>
#include <sstream>

#include "BlasInterface.hpp"
#include "LapackInterface2.hpp"

using namespace TensorCalculus;

// const int N = 100;
// const int M = 20;

void print_matrix(std::ostream& stream, const double *ptr, int height, int width)
{
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      stream << ptr[i + j*height] << ' ';
    }
    stream << '\n';
  }
}
                         
int main(int argc, const char* argv[])
{
  if (argc < 3) {
    std::cout << "Usage: LapackTest <matrix-width> <num-evals>" << std::endl;
    return 0;
  }
  int N;
  int M;
  {
    std::string args(argv[1]);
    std::istringstream iss(args);
    iss >> N;
    if (!iss) {
      std::cout << "Error: Can't read matrix width" << std::endl;
      return -1;
    }
  }
  {
    std::string args(argv[2]);
    std::istringstream iss(args);
    iss >> M;
    if (!iss) {
      std::cout << "Error: Can't read number of eigenvalues" << std::endl;
      return -1;
    }
  }

  std::cout << "Fill matrix A...";
  std::cout.flush();
  std::vector<double> A(N * N, 0.0);
  {
    std::srand(0);
    std::vector<double> x(N);
    for (int j = 1; j <= N; ++j) {
      for (int i = 0; i < N; ++i) {
        x[i] = static_cast<double>(std::rand()) / RAND_MAX;
      }
      Blas<double>::ger(N, N, j, &x[0], 1, &x[0], 1, &A[0], N);
    }
  }
  std::cout << "OK" << std::endl;
  
  std::cout << "Fill matrix B...";
  std::cout.flush();
  std::vector<double> B(N * N);
  {
    std::vector<double> pre_B(N * N, 0.0);
    for (int i = 0; i < N * N; ++i) {
      pre_B[i] = static_cast<double>(std::rand()) / RAND_MAX;
    }
    Blas<double>::gemm('T', 'N', N, N, N, 1.0, &pre_B[0], N, &pre_B[0], N, 0.0, &B[0], N);
  }
  std::cout << "OK" << std::endl;
  
  std::cout << "Write A and B to the files...";
  std::cout.flush();
  {
    std::ofstream A_file("A.mat");
    print_matrix(A_file, &A[0], N, N);
  }
  {
    std::ofstream B_file("B.mat");
    print_matrix(B_file, &B[0], N, N);
  }
  std::cout << "OK" << std::endl;

  std::cout << "Calculate eigenvalues and -vectors...";
  std::cout.flush();
  int m;
  std::vector<double> eval(N);
  std::vector<double> evec(N * M);
  {
    double dbl_lwork;
    std::vector<int> iwork(5 * N);
    std::vector<int> ifail(N);
    int result = 
      Lapack<double>::sygvx(2, 'V', 'I', 'U',
                            N, &A[0], N, &B[0], N,
                            0.0, 0.0, N - M + 1, N,
                            -1.0,
                            m, &eval[0], &evec[0], N,
                            &dbl_lwork, -1, &iwork[0], &ifail[0]);
    if (result != 0) {
      std::cout << "Query failed! Error: " << result << std::endl;
      return -1;
    }
    int lwork = static_cast<int>(dbl_lwork);
    std::vector<double> work(lwork);
    result =
      Lapack<double>::sygvx(2, 'V', 'I', 'U',
                            N, &A[0], N, &B[0], N,
                            0.0, 0.0, N - M + 1, N,
                            -1.0,
                            m, &eval[0], &evec[0], N,
                            &work[0], lwork, &iwork[0], &ifail[0]);
    if (result != 0) {
      std::cout << "Failed! Error: " << result << std::endl;
      return -1;
    }
  }
  std::cout << "OK" << std::endl;
  
  std::cout << "Write eigenvalues and -vectors to the files...";
  std::cout.flush();
  {
    std::ofstream eval_file("eval.mat");
    print_matrix(eval_file, &eval[0], N, 1);
  }
  {
    std::ofstream evec_file("evec.mat");
    print_matrix(evec_file, &evec[0], N, M);
  }
  std::cout << "OK" << std::endl;
  
  return 0;
}
