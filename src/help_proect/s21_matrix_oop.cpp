#include "s21_matrix_oop.hpp"

#include <iostream>
#include <algorithm>
#include <cmath>

// TODO more comprehensive exception messages?
// I thought the call stack is printed to stderr, but that's not the case,
// which is weird.

// TODO fix wasteful exceptions?
// Take a look at CalcComplements(), for example. It throws. It calls to,
// Minor(), which also throws. Minor() calls to Determinant(), which also
// throws. Does the compiler optimize all this stuff? Otherwise, it's kinda
// wasteful.

S21Matrix::S21Matrix() {
  _rows = 4;
  _cols = 4;
  _matrix = new double[16]{};

  for (int i = 0; i != 16; ++i) {
    _matrix[i * _cols + i] = 1;
  }
}

S21Matrix::S21Matrix(int rows, int cols) {
  if (rows < 0 || cols < 0) {
    throw std::invalid_argument("S21Matrix: rows or cols less than zero");
  }
  if (rows == 0 || cols == 0) {
    throw std::invalid_argument("S21Matrix: empty matrix");
  }

  _rows = rows;
  _cols = cols;

  _matrix = new double[static_cast<long>(_rows) * _cols]{};

  int min_dim = std::min(_rows, _cols);
  for (int i = 0; i != min_dim; ++i) {
    _matrix[i * _cols + i] = 1;
  }

  if (log_lifetimes) std::cerr << id << ": created!\n";
}

// The caller MUST provide correct matrix, rows and cols. Because we have no way
// to tell if matrix has exactly rows * cols * sizeof(double) allocated memory.
S21Matrix::S21Matrix(const double* matrix, int rows, int cols) {
  if (matrix == nullptr) {
    throw std::invalid_argument("S21Matrix: can't create matrix from nullptr");
  }
  if (rows < 0 || cols < 0) {
    throw std::invalid_argument("S21Matrix: rows or cols less than zero");
  }
  if (rows == 0 || cols == 0) {
    throw std::invalid_argument("S21Matrix: empty matrix");
  }
  _rows = rows;
  _cols = cols;
  _matrix = (double*)matrix;

  if (log_lifetimes) std::cerr << id << ": created!\n";
}

// Copy constructor
S21Matrix::S21Matrix(const S21Matrix& other) {
  _rows = other._rows;
  _cols = other._cols;

  _matrix = new double[static_cast<long>(_rows) * _cols]{};
  std::copy_n(other._matrix, _rows * _cols, _matrix);

  if (log_lifetimes) std::cerr << id << ": copied!\n";
}

// Move constructor
S21Matrix::S21Matrix(S21Matrix&& other) noexcept {
  _rows = other._rows;
  _cols = other._cols;
  _matrix = other._matrix;

  // other._rows = 0;
  // other._cols = 0;
  other._matrix = nullptr;

  if (log_lifetimes) std::cerr << id << ": moved!\n";
}

// Destructor
S21Matrix::~S21Matrix() noexcept {
  delete[] _matrix;
  if (log_lifetimes) std::cerr << id << ": destroyed!\n";
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& rhs) {
  SumMatrix(rhs);
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& rhs) {
  SubMatrix(rhs);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& rhs) {
  MulMatrix(rhs);
  return *this;
}

S21Matrix& S21Matrix::operator*=(double rhs) {
  MulNumber(rhs);
  return *this;
}

S21Matrix S21Matrix::operator+(const S21Matrix& rhs) const {
  S21Matrix sum = *this;
  sum += rhs;
  return sum;
}

S21Matrix S21Matrix::operator-(const S21Matrix& rhs) const {
  S21Matrix sum = *this;
  sum -= rhs;
  return sum;
}

S21Matrix S21Matrix::operator*(const S21Matrix& rhs) const {
  S21Matrix sum = *this;
  sum *= rhs;
  return sum;
}

S21Matrix S21Matrix::operator*(double rhs) const {
  S21Matrix sum = *this;
  sum *= rhs;
  return sum;
}

bool S21Matrix::operator==(const S21Matrix& rhs) const { return EqMatrix(rhs); }

bool S21Matrix::operator!=(const S21Matrix& rhs) const {
  return !EqMatrix(rhs);
}

double& S21Matrix::operator()(int row, int col) const {
  if (row < 0 || col < 0 || row > _rows - 1 || col > _cols - 1) {
    throw std::out_of_range("S21Matrix: matrix indices out of range");
  }

  return _matrix[row * _cols + col];
}

// Because using operator() for indexing is stupid
// https://stackoverflow.com/q/317450
double& S21Matrix::at(int row, int col) const { return (*this)(row, col); }

bool S21Matrix::EqMatrix(const S21Matrix& other) const {
  if (_rows != other._rows || _cols != other._cols) {
    return false;
  }

  long n_elements = static_cast<long>(_rows) * _cols;

  // TODO maybe try to use iterators instead of raw loops?
  for (long i = 0; i < n_elements; ++i) {
    if (_matrix[i] != other._matrix[i]) {
      return false;
    }
  }

  return true;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (_rows != other._rows || _cols != other._cols) {
    throw std::invalid_argument("S21Matrix: different matrix dimensions");
  }

  long n_elements = static_cast<long>(_rows) * _cols;

  for (long i = 0; i < n_elements; ++i) {
    _matrix[i] += other._matrix[i];
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (_rows != other._rows || _cols != other._cols) {
    throw std::invalid_argument("S21Matrix: different matrix dimensions");
  }

  long n_elements = static_cast<long>(_rows) * _cols;

  for (long i = 0; i < n_elements; ++i) {
    _matrix[i] -= other._matrix[i];
  }
}

void S21Matrix::MulNumber(const double num) noexcept {
  long n_elements = static_cast<long>(_rows) * _cols;

  for (long i = 0; i < n_elements; ++i) {
    _matrix[i] *= num;
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (_cols != other._rows) {
    throw std::invalid_argument(
        "S21Matrix: can't multiply matrices with such dimensions");
  }

  double* result = new double[static_cast<long>(_rows) * other._cols]{};

  for (int i = 0; i < _rows; ++i) {
    for (int j = 0; j < other._cols; ++j) {
      for (int k = 0; k < _cols; ++k) {
        result[i * other._cols + j] +=
            _matrix[i * _cols + k] * other._matrix[k * other._cols + j];
      }
    }
  }

  _cols = other._cols;
  _matrix = result;
}

S21Matrix S21Matrix::Transpose() const {
  S21Matrix transposed(_cols, _rows);

  for (int i = 0; i < _rows; ++i) {
    for (int j = 0; j < _cols; ++j) {
      transposed.at(j, i) = _matrix[i * _cols + j];
    }
  }

  return transposed;
}

// don't ask
void S21Matrix::remove_row_col(const double* mat, int mat_size, double* target,
                               int row_to_skip, int col_to_skip) const {
  int i = 0;

  for (int row = 0; row < mat_size; ++row) {
    if (row == row_to_skip) continue;
    for (int col = 0; col < mat_size; ++col) {
      if (col == col_to_skip) continue;
      // Row-major, remember. We're just filling target[] row by row.
      target[i] = mat[row * mat_size + col];
      i++;
    }
  }
}

// Copied from:
// https://www.tutorialspoint.com/determinant-of-a-matrix-in-cplusplus-program
// I needed it to "just work".
double S21Matrix::det(const double* mat, int n) const {
  double determinant = 0;
  if (n == 1) return mat[0];
  if (n == 2) return mat[0] * mat[3] - mat[1] * mat[2];

  double* tmp = new double[static_cast<unsigned long>(n - 1) * (n - 1)]{};
  int sign = 1;
  for (int i = 0; i < n; ++i) {
    remove_row_col(mat, n, tmp, 0, i);
    determinant += sign * mat[i] * det(tmp, n - 1);
    sign = -sign;
  }

  delete[] tmp;
  return determinant;
}

double S21Matrix::Determinant() const {
  if (_rows != _cols) {
    throw std::logic_error("S21Matrix: not a square matrix");
  }

  return det(_matrix, _rows);
}

double S21Matrix::Minor(int row, int col) const {
  if (row < 0 || col < 0 || row > _rows - 1 || col > _cols - 1) {
    throw std::out_of_range("S21Matrix: matrix indices out of range");
  }
  if (_rows != _cols) {
    throw std::logic_error("S21Matrix: not a square matrix");
  }

  int n = _rows;

  double* tmp = new double[static_cast<unsigned long>(n - 1) * (n - 1)]{};
  remove_row_col(_matrix, n, tmp, row, col);

  double minor = det(tmp, n - 1);
  delete[] tmp;
  return minor;
}

S21Matrix S21Matrix::CalcComplements() const {
  if (_rows != _cols) {
    throw std::logic_error("S21Matrix: not a square matrix");
  }

  int n = _rows;
  double* complements = new double[static_cast<unsigned long>(n) * n];

  for (int row = 0; row < n; ++row) {
    for (int col = 0; col < n; ++col) {
      int sign = (row + col) % 2 == 0 ? 1 : -1;
      complements[row * n + col] = sign * Minor(row, col);
    }
  }

  return S21Matrix(complements, n, n);
}

S21Matrix S21Matrix::InverseMatrix() const {
  if (_rows != _cols) {
    throw std::logic_error("S21Matrix: not a square matrix");
  }

  double determinant = Determinant();
  if (std::fabs(determinant) < 1e-7) {
    throw std::logic_error("S21Matrix: a singular matrix");
  }

  S21Matrix inverse = CalcComplements().Transpose();
  inverse.MulNumber(1 / determinant);
  return inverse;
}

int S21Matrix::GetRows() const { return _rows; }

int S21Matrix::GetCols() const { return _cols; }

double* S21Matrix::GetRaw() const { return _matrix; }

void S21Matrix::SetRows(int new_rows) {
  if (new_rows < 1) {
    throw std::out_of_range("S21Matrix: rows can't be set to 0 or less");
  }
  // The compiler will probably optimize this anyway, but just in case
  if (new_rows == _rows) {
    return;
  }

  double* new_matrix = new double[static_cast<long>(new_rows) * _cols]{};
  std::copy_n(_matrix, std::min(new_rows, _rows) * _cols, new_matrix);

  _rows = new_rows;
  delete[] _matrix;
  _matrix = new_matrix;
}

void S21Matrix::SetCols(int new_cols) {
  if (new_cols < 1) {
    throw std::out_of_range("S21Matrix: cols can't be set to 0 or less");
  }
  // The compiler will probably optimize this anyway, but just in case
  if (new_cols == _cols) {
    return;
  }

  double* new_matrix = new double[static_cast<long>(_rows) * new_cols]{};

  int cols_to_copy = std::min(new_cols, _cols);

  for (int i = 0; i != _rows; ++i) {
    for (int j = 0; j != cols_to_copy; ++j) {
      new_matrix[i * new_cols + j] = _matrix[i * _cols + j];
    }
  }

  _cols = new_cols;
  delete[] _matrix;
  _matrix = new_matrix;
}

// TODO probably remove this
void S21Matrix::print() const {
  std::cout << "[";
  for (int i = 0; i != _rows; ++i) {
    std::cout << "[";
    for (int j = 0; j != _cols; ++j) {
      std::cout << _matrix[i * _cols + j] << ", ";
    }
    std::cout << "], ";
  }
  std::cout << "]\n";
}

