#include "s21_matrix_oop.hpp"

#include <iostream>

int main() {
  try {
    double *m1_raw = new double[4]{1, -2, 3, -4};
    double *m2_raw = new double[4]{1, 2, 3, 4};
    S21Matrix m1(m1_raw, 2, 2);
    S21Matrix m2(m2_raw, 2, 2);
    S21Matrix m3 = m1;

    m1.SetRows(3);
    std::cout << "m1 == m3: " << (m1 == m3) << "\n";

  } catch (const std::bad_alloc &ex) {
    std::cout << "ERROR: S21Matrix: " << ex.what() << "\n";
  }
  return 0;
}
