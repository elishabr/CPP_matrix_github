#ifndef __S21MATRIX_H__
#define __S21MATRIX_H__

#include <iostream>
#include <cstddef>

#define ERR_MAT "S21Matrix: rows or cols less than zero"

class S21Matrix {

public:
    S21Matrix();
    S21Matrix(int rows, int cols);
    S21Matrix(S21Matrix&& other) noexcept;
    S21Matrix(const S21Matrix& other);
    ~S21Matrix();

    S21Matrix& operator=(const S21Matrix& other);
    double& operator()(int row, int col);

    bool EqMatrix(const S21Matrix& other);
    void SumMatrix(const S21Matrix& other);

    //Accessor 
    int getRows() const;
    int getCols() const;
    //Mutator
    void setRows(int newRows);
    void setCols(int newCols);
private:
    // attributes
    int rows_;
    int cols_;
    double **matrix_;
    // private metod
    void swap_(S21Matrix& other);

    //Debug
    static const bool log_lifetimes = true;

};

#endif