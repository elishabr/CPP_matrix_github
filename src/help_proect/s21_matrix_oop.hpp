#ifndef S21_MATRIX_OOP_H
#define S21_MATRIX_OOP_H

#include <vector>

class S21Matrix {
 private:
  // For debug purposes
  static const bool log_lifetimes = true;
  static inline int new_id = 0;
  const int id = ++new_id;

  int _rows = 4;
  int _cols = 4;
  
  // Could have used std::vector, but that would be too easy!
  // NOTE: storing in row-major format
  double* _matrix = new double[16]{};

  // Some helper computation functions
  void remove_row_col(const double* mat, int mat_size, double* tmp,
                      int row_to_skip, int col_to_skip) const;
  double det(const double* mat, int n) const;

 public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const double matrix[], int rows, int cols);

  // Copy constructor
  S21Matrix(const S21Matrix& other);
  // Move constructor
  S21Matrix(S21Matrix&& other) noexcept;
  // Destructor
  ~S21Matrix() noexcept;

  S21Matrix& operator+=(const S21Matrix& rhs);
  S21Matrix& operator-=(const S21Matrix& rhs);
  S21Matrix& operator*=(const S21Matrix& rhs);
  S21Matrix& operator*=(double rhs);
  S21Matrix operator+(const S21Matrix& rhs) const;
  S21Matrix operator-(const S21Matrix& rhs) const;
  S21Matrix operator*(const S21Matrix& rhs) const;
  S21Matrix operator*(double rhs) const;
  bool operator==(const S21Matrix& rhs) const;
  bool operator!=(const S21Matrix& rhs) const;
  double & operator()(int row, int col) const;

  double & at(int row, int col) const;

  bool EqMatrix(const S21Matrix& other) const;
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(double num) noexcept;
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose() const;
  double Determinant() const;
  double Minor(int row, int col) const;
  S21Matrix CalcComplements() const;
  S21Matrix InverseMatrix() const;

  int GetRows() const;
  int GetCols() const;
  double* GetRaw() const;
  void SetRows(int rows);
  void SetCols(int cols);

  void print() const;
};

#endif
