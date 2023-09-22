#include "lowtriangtoeplitz.h"

int main() {
  LowTriangToeplitz::test_accuracy_ltpMult();
  LowTriangToeplitz::time_measure_ltpMult();

  LowTriangToeplitz::test_accuracy_ltpSolve();
  LowTriangToeplitz::time_measure_ltpSolve();
  return 0;
}