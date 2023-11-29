#include "s21_matrix_oop.hpp"

S21Matrix::S21Matrix(): rows_(3), cols_(3) {
    matrix_ = new double * [rows_];
    matrix_[0] = new double [rows_ * cols_]{};
    for(int i = 1; i !=rows_; i++){
        matrix_[i] = matrix_[i-1]+ cols_;
    }
    if (log_lifetimes) std::cerr << "ConstructorsDefault\n";
}

S21Matrix::S21Matrix(int rows, int cols): rows_(rows), cols_(cols) {
    if(rows <=0 || cols<=0){
        throw std::invalid_argument(ERR_MAT);
    }
    matrix_ = new double * [rows_];
    matrix_[0] = new double [rows_ * cols_]{};
    for(int i = 1; i !=rows_; i++){
        matrix_[i] = matrix_[i-1]+ cols_;
    }
    if (log_lifetimes) std::cerr << "ConstructorsRowsCols\n";
}

S21Matrix::~S21Matrix(){
    if(matrix_){
        delete [] matrix_[0];
        delete [] matrix_;
    }
    if (log_lifetimes) std::cerr << "Destructor\n";
}

S21Matrix::S21Matrix(const S21Matrix& other): rows_(other.rows_), cols_(other.cols_) {
    matrix_ = new double * [rows_];
    matrix_[0] = new double [rows_ * cols_]{};
    for(int i = 1; i !=rows_; i++){
        matrix_[i] = matrix_[i-1]+ cols_;
    }
    for(int i=0; i != other.rows_; i++){
        for(int j=0; j!= other.cols_; j++){
            matrix_[i][j] = other.matrix_[i][j];
        }
    }
    if (log_lifetimes) std::cerr << "ConstructorsCopy\n";
}

S21Matrix::S21Matrix(S21Matrix&& other) noexcept {
    rows_ = other.rows_;
    cols_ = other.cols_;
    matrix_ = other.matrix_;

    other.matrix_ = nullptr;
    if (log_lifetimes) std::cerr << "ConstructorsMove\n";
}

S21Matrix & S21Matrix::operator=(const S21Matrix & other){
    if(this != &other)
        S21Matrix(other).swap_(*this);
    if (log_lifetimes) std::cerr << "Operator=\n";
    return *this;
}

void S21Matrix::swap_(S21Matrix &other){
	std::swap(rows_, other.rows_);
	std::swap(cols_, other.cols_);
    std::swap(matrix_, other.matrix_);
}

double & S21Matrix::operator()(int row, int col){
    if(row >= this->rows_ || col >= this->cols_ || row < 0 || col < 0){
        throw std::out_of_range("S21Matrix: matrix indices out of range");
    }
    if (log_lifetimes) std::cerr << "Operator()\n";
    return this->matrix_[row][col];
}

bool S21Matrix::EqMatrix(const S21Matrix& other){
    if(rows_ != other.rows_ || cols_ != other.cols_){
        return false;
    }
    for(int i=0; i != rows_; i++){
        for(int j=0; j != cols_; j++){
            if(matrix_[i][j] != other.matrix_[i][j]) return false;
        }
    }
    return true;
}

void S21Matrix::SumMatrix(const S21Matrix& other){
    if(){

    }
}

// void S21Matrix::SumMatrix(const S21Matrix& other) {
//   if (_rows != other._rows || _cols != other._cols) {
//     throw std::invalid_argument("S21Matrix: different matrix dimensions");
//   }

//   long n_elements = static_cast<long>(_rows) * _cols;

//   for (long i = 0; i < n_elements; ++i) {
//     _matrix[i] += other._matrix[i];
//   }
// }

int S21Matrix::getRows() const { return this->rows_; }

int S21Matrix::getCols() const { return this->cols_;}

void S21Matrix::setRows(int newRows){this->rows_ = newRows;}

void S21Matrix::setCols(int newCols){ this->cols_ = newCols;}