#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "lowtriangtoeplitz.h"

int main() {

#if SOLUTION
  // Compare runtime of matrix multiplication
  std::cout << "\nMatrix size \t";
  std::cout << "Runtime for M-M \t";
  std::cout << "Runtime for M-V \t";
  std::cout << "Runtime for V-V" << std::endl;
  std::ofstream f_mult(CURRENT_SOURCE_DIR "/runtime_mult.txt");
  for (int p = 3; p < 10; p++) // take too long for p > 10
  {
    int N = pow(2, p);
    auto [dense, mv, ltp] = LowTriangToeplitz::runtimes_ltpMult(N);
    std::cout << N << "\t\t" << dense << "\t\t" << mv << "\t\t" << ltp << std::endl;
    f_mult << N << " " << dense << " " << mv << " " << ltp << std::endl;
  }
  f_mult.close();
  
  // Compare runtime of solving linear system of equations
  std::cout << "\nMatrix size \t";
  std::cout << "Runtime for tria \t";
  std::cout << "Runtime for ltp"<< std::endl;
  std::ofstream f_solve(CURRENT_SOURCE_DIR "/runtime_solve.txt");
  for (int p = 3; p < 10; p++) // take too long for p > 10
  {
    int N = pow(2, p);
    auto [tria, ltp] = LowTriangToeplitz::runtimes_ltpSolve(N);
    std::cout << N << "\t\t" << tria << "\t\t" << ltp << std::endl;
    f_solve << N << " " << tria << " " << ltp << std::endl;
  }
  f_solve.close();

  // Call python script
  std::system("python3 " CURRENT_SOURCE_DIR "/plot.py " CURRENT_SOURCE_DIR "/runtime_mult.txt " 
                CURRENT_SOURCE_DIR "/runtime_solve.txt " CURRENT_SOURCE_DIR "/runtime.png");
#else
  // **********************************************************************
  // Code to be supplemented
  // **********************************************************************
#endif

  return 0;
}