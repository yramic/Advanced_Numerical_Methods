#include <stdlib.h>

#include <fstream>
#include <iomanip>
#include <iostream>

#include "lowtriangtoeplitz.h"

int main() {
  // Compare runtime of matrix multiplication
  std::cout << std::left << std::setw(20) << "\nMatrix size " << std::left
            << std::setw(20) << "Runtime for M-M" << std::left << std::setw(20)
            << "Runtime for M-V" << std::left << std::setw(20)
            << "Runtime for V-V" << std::endl;
  std::ofstream f_mult(CURRENT_SOURCE_DIR "/runtime_mult.txt");
  for (int p = 3; p < 10; p++)  // take too long for p > 10
  {
    const int N = pow(2, p);
    auto [dense, mv, ltp] = LowTriangToeplitz::runtimes_ltpMult(N);
    std::cout << std::left << std::setw(20) << N << std::left << std::setw(20)
              << dense << std::left << std::setw(20) << mv << std::left
              << std::setw(20) << ltp << std::endl;
    f_mult << N << " " << dense << " " << mv << " " << ltp << std::endl;
  }
  f_mult.close();

  // Compare runtime of solving linear system of equations
  std::cout << std::left << std::setw(20) << "\nMatrix size" << std::left
            << std::setw(20) << "Runtime for tria" << std::left << std::setw(20)
            << "Runtime for ltp" << std::endl;
  std::ofstream f_solve(CURRENT_SOURCE_DIR "/runtime_solve.txt");
  for (int p = 5; p < 14; p++)  // take too long for p > 10
  {
    int N = pow(2, p);
    auto [tria, ltp] = LowTriangToeplitz::runtimes_ltpSolve(N);
    std::cout << std::left << std::setw(20) << N << std::left << std::setw(20)
              << tria << std::left << std::setw(20) << ltp << std::endl;
    f_solve << N << " " << tria << " " << ltp << std::endl;
  }
  f_solve.close();

  // Call python script
  std::system("python3 " CURRENT_SOURCE_DIR "/plot.py " CURRENT_SOURCE_DIR
              "/runtime_mult.txt " CURRENT_SOURCE_DIR
              "/runtime_solve.txt " CURRENT_SOURCE_DIR "/runtime.png");

  return 0;
}