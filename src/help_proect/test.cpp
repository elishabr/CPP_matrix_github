// By no way these tests are exhaustive. I just had to create *some* tests.

#include "s21_matrix_oop.hpp"

#include <cmath>
#include <gtest/gtest.h>

const double EPSILON = 1e-7;

bool almost_equal(double a, double b) { return std::fabs(a - b) < EPSILON; }

bool is_matrix_equal(const S21Matrix &m, const double *reference) {
  const double *m_raw = m.GetRaw();
  long n_elements = static_cast<long>(m.GetRows()) * m.GetCols();

  for (long i = 0; i < n_elements; ++i) {
    if (!almost_equal(m_raw[i], reference[i])) {
      return false;
    }
  }

  return true;
}

const double identity_4_by_4[16] = {1, 0, 0, 0, 0, 1, 0, 0,
                                    0, 0, 1, 0, 0, 0, 0, 1};

const double identity_2_by_3[6] = {1, 0, 0, 0, 1, 0};

TEST(Constructor, Default) {
  S21Matrix m;

  EXPECT_EQ(m.GetRows(), 4);
  EXPECT_EQ(m.GetCols(), 4);

  EXPECT_TRUE(is_matrix_equal(m, identity_4_by_4));
}

TEST(Constructor, RowsCols) {
  S21Matrix m(2, 3);

  EXPECT_EQ(m.GetRows(), 2);
  EXPECT_EQ(m.GetCols(), 3);

  EXPECT_TRUE(is_matrix_equal(m, identity_2_by_3));
}

TEST(Constructor, PreAllocatedMatrix) {
  double *m_raw =
      new double[16]{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  S21Matrix m(m_raw, 4, 4);

  EXPECT_EQ(m.GetRows(), 4);
  EXPECT_EQ(m.GetCols(), 4);

  EXPECT_TRUE(is_matrix_equal(m, m_raw));
}

TEST(Constructor, Copy) {
  S21Matrix m1(2, 3);
  S21Matrix m2 = m1;

  EXPECT_EQ(m2.GetRows(), 2);
  EXPECT_EQ(m2.GetCols(), 3);

  m1.at(0, 0) = 123;

  EXPECT_TRUE(is_matrix_equal(m2, identity_2_by_3));

  m2.at(0, 0) = 42;
  const double reference[6] = {42, 0, 0, 0, 1, 0};
  EXPECT_TRUE(is_matrix_equal(m2, reference));
}

TEST(Constructor, Move) {
  S21Matrix m1(2, 3);
  S21Matrix m2 = m1;
  S21Matrix m3(std::move(m1));

  EXPECT_EQ(m1.GetRaw(), nullptr);
  EXPECT_NE(m2.GetRaw(), m3.GetRaw());

  EXPECT_EQ(m2.GetRows(), 2);
  EXPECT_EQ(m2.GetCols(), 3);

  EXPECT_EQ(m3.GetRows(), 2);
  EXPECT_EQ(m3.GetCols(), 3);
}

TEST(Operator, Parenthesis) {
  S21Matrix m(2, 3);

  EXPECT_EQ(m(1, 1), 1);
  EXPECT_EQ(m(1, 0), 0);

  m(0, 0) = 123;

  EXPECT_EQ(m(0, 0), 123);
}

TEST(Function, At) {
  S21Matrix m(2, 3);

  EXPECT_EQ(m.at(1, 1), 1);
  EXPECT_EQ(m.at(1, 0), 0);

  m.at(0, 0) = 123;

  EXPECT_EQ(m.at(0, 0), 123);
}

TEST(Function, EqMatrix) {
  S21Matrix m1(2, 3);
  S21Matrix m2(2, 3);
  S21Matrix m3(4, 3);
  EXPECT_TRUE(m1.EqMatrix(m2));
  EXPECT_FALSE(m1.EqMatrix(m3));

  m2.at(0, 1) = 123;
  EXPECT_FALSE(m1.EqMatrix(m2));
}

TEST(Function, SumMatrix) {
  S21Matrix m1(2, 3);
  S21Matrix m2(2, 3);

  m1.SumMatrix(m2);
  double reference[6] = {2, 0, 0, 0, 2, 0};
  EXPECT_TRUE(is_matrix_equal(m1, reference));

  S21Matrix m3(3, 3);
  EXPECT_THROW(m1.SumMatrix(m3), std::invalid_argument);
}

TEST(Function, SubMatrix) {
  S21Matrix m1(2, 3);
  S21Matrix m2(2, 3);

  m1.SubMatrix(m2);
  double reference[6] = {0, 0, 0, 0, 0, 0};
  EXPECT_TRUE(is_matrix_equal(m1, reference));

  S21Matrix m3(3, 3);
  EXPECT_THROW(m1.SubMatrix(m3), std::invalid_argument);
}

TEST(Function, MulNumber) {
  S21Matrix m(2, 3);

  m.MulNumber(42);
  double reference[6] = {42, 0, 0, 0, 42, 0};
  EXPECT_TRUE(is_matrix_equal(m, reference));
}

TEST(Function, MulMatrix) {
  double *m1_raw = new double[6]{2, 3, 4, 5, 6, 7};
  double *m2_raw = new double[6]{3, 7, 11, 13, 23, 0};

  S21Matrix m1(m1_raw, 2, 3);
  S21Matrix m2(m2_raw, 3, 2);

  m1.MulMatrix(m2);
  EXPECT_EQ(m1.GetRows(), 2);
  EXPECT_EQ(m1.GetCols(), 2);

  double reference[4] = {131, 53, 242, 113};
  EXPECT_TRUE(is_matrix_equal(m1, reference));

  EXPECT_THROW(m1.MulMatrix(m2), std::invalid_argument);
}

TEST(Function, Transpose) {
  double *m_raw = new double[6]{2, 3, 4, 5, 6, 7};

  S21Matrix m(m_raw, 2, 3);

  S21Matrix transposed = m.Transpose();
  EXPECT_EQ(transposed.GetRows(), 3);
  EXPECT_EQ(transposed.GetCols(), 2);

  double reference[6] = {2, 5, 3, 6, 4, 7};
  EXPECT_TRUE(is_matrix_equal(transposed, reference));
}

TEST(Function, Determinant) {
  double *m_raw = new double[9]{1, 2, -3, 4, -5, 6, 7, 8, -9};

  S21Matrix m(m_raw, 3, 3);
  EXPECT_TRUE(almost_equal(m.Determinant(), -48));

  m.SetRows(2);
  EXPECT_THROW(m.Determinant(), std::logic_error);
}

TEST(Function, Minor) {
  double *m_raw = new double[9]{2, 2, 1, 7, 8, 2, 5, 3, 2};

  S21Matrix m(m_raw, 3, 3);
  EXPECT_TRUE(almost_equal(m.Minor(0, 1), 4));
  EXPECT_THROW(m.Minor(42, 0), std::out_of_range);

  m.SetRows(2);
  EXPECT_THROW(m.Minor(0, 1), std::logic_error);
}

TEST(Function, CalcComplements) {
  double *m_raw = new double[9]{2, 2, 1, 7, 8, 2, 5, 3, 2};

  S21Matrix m(m_raw, 3, 3);
  S21Matrix complements = m.CalcComplements();
  double reference[9] = {10, -4, -19, -1, -1, 4, -4, 3, 2};
  EXPECT_TRUE(is_matrix_equal(complements, reference));

  m.SetRows(2);
  EXPECT_THROW(m.CalcComplements(), std::logic_error);
}

TEST(Function, Inverse) {
  double *m_raw = new double[9]{2, 2, 1, 7, 8, 2, 5, 3, 2};

  S21Matrix m1(m_raw, 3, 3);
  S21Matrix inverse = m1.InverseMatrix();
  double reference[9] = {-10.0 / 7, 1.0 / 7,  4.0 / 7,  4.0 / 7, 1.0 / 7,
                         -3.0 / 7,  19.0 / 7, -4.0 / 7, -2.0 / 7};
  EXPECT_TRUE(is_matrix_equal(inverse, reference));

  m1.SetRows(2);
  EXPECT_THROW(m1.InverseMatrix(), std::logic_error);

  S21Matrix m2(3, 3);
  EXPECT_THROW(m1.InverseMatrix(), std::logic_error);
}

TEST(Operator, Plus) {
  double *m_raw1 = new double[4]{1, 2, 3, 4};
  double *m_raw2 = new double[4]{41, 42, 43, 44};

  S21Matrix m1(m_raw1, 2, 2);
  S21Matrix m2(m_raw2, 2, 2);

  S21Matrix m3 = m1 + m2;
  double reference[4] = {42, 44, 46, 48};
  EXPECT_TRUE(is_matrix_equal(m3, reference));

  m1 += m2;
  EXPECT_TRUE(is_matrix_equal(m1, reference));

  m1.SetRows(3);
  EXPECT_THROW({ m1 += m2; }, std::invalid_argument);
  EXPECT_THROW({ S21Matrix m4 = m1 + m2; }, std::invalid_argument);
}

TEST(Operator, Minus) {
  double *m_raw1 = new double[4]{41, 42, 43, 44};
  double *m_raw2 = new double[4]{1, 5, 10, 15};

  S21Matrix m1(m_raw1, 2, 2);
  S21Matrix m2(m_raw2, 2, 2);

  S21Matrix m3 = m1 - m2;
  double reference[4] = {40, 37, 33, 29};
  EXPECT_TRUE(is_matrix_equal(m3, reference));

  m1 -= m2;
  EXPECT_TRUE(is_matrix_equal(m1, reference));

  m1.SetRows(3);
  EXPECT_THROW({ m1 -= m2; }, std::invalid_argument);
  EXPECT_THROW({ S21Matrix m4 = m1 - m2; }, std::invalid_argument);
}

TEST(Operator, MultiplyByNumber) {
  // By a matrix
  double *m_raw = new double[4]{1, -2, 3, -4};

  S21Matrix m1(m_raw, 2, 2);

  S21Matrix m2 = m1 * 2;
  double reference[4] = {2, -4, 6, -8};
  EXPECT_TRUE(is_matrix_equal(m2, reference));

  m1 *= 2;
  EXPECT_TRUE(is_matrix_equal(m1, reference));
}

TEST(Operator, MultiplyByMatrix) {
  // By a matrix
  double *m_raw1 = new double[4]{1, 2, 3, 4};
  double *m_raw2 = new double[4]{-3, 4, -5, 6};

  S21Matrix m1(m_raw1, 2, 2);
  S21Matrix m2(m_raw2, 2, 2);

  S21Matrix m3 = m1 * m2;
  double reference[4] = {-13, 16, -29, 36};
  EXPECT_TRUE(is_matrix_equal(m3, reference));

  m1 *= m2;
  EXPECT_TRUE(is_matrix_equal(m1, reference));

  m1.SetCols(3);
  EXPECT_THROW({ m1 *= m2; }, std::invalid_argument);
  EXPECT_THROW({ S21Matrix m4 = m1 * m2; }, std::invalid_argument);
}

TEST(Logic, Equal) {
  double *m1_raw = new double[4]{1, -2, 3, -4};
  double *m2_raw = new double[4]{1, 2, 3, 4};
  S21Matrix m1(m1_raw, 2, 2);
  S21Matrix m2(m2_raw, 2, 2);
  S21Matrix m3 = m1;

  EXPECT_FALSE(m1 == m2);
  EXPECT_TRUE(m1 == m3);
  EXPECT_TRUE(m1 == m1);

  m1.SetRows(3);
  EXPECT_FALSE(m1 == m3);
}

TEST(Logic, NotEqual) {
  double *m1_raw = new double[4]{1, -2, 3, -4};
  double *m2_raw = new double[4]{1, 2, 3, 4};
  S21Matrix m1(m1_raw, 2, 2);
  S21Matrix m2(m2_raw, 2, 2);
  S21Matrix m3 = m1;

  EXPECT_TRUE(m1 != m2);
  EXPECT_FALSE(m1 != m3);
  EXPECT_FALSE(m1 != m1);

  m1.SetRows(3);
  EXPECT_TRUE(m1 != m3);
}

TEST(Setter, SetRow) {
  S21Matrix m(2, 3);

  m.SetRows(3);
  EXPECT_EQ(m.GetRows(), 3);
  EXPECT_EQ(m.GetCols(), 3);

  double reference1[9] = {1, 0, 0, 0, 1, 0, 0, 0, 0};
  EXPECT_TRUE(is_matrix_equal(m, reference1));

  m.at(2, 2) = 123;
  reference1[8] = 123;
  EXPECT_TRUE(is_matrix_equal(m, reference1));

  m.SetRows(1);
  EXPECT_EQ(m.GetRows(), 1);
  EXPECT_EQ(m.GetCols(), 3);

  double reference2[3] = {1, 0, 0};
  EXPECT_TRUE(is_matrix_equal(m, reference2));
}

TEST(Setter, SetCol) {
  S21Matrix m(3, 2);

  m.SetCols(3);
  EXPECT_EQ(m.GetRows(), 3);
  EXPECT_EQ(m.GetCols(), 3);

  double reference1[9] = {1, 0, 0, 0, 1, 0, 0, 0, 0};
  EXPECT_TRUE(is_matrix_equal(m, reference1));

  m.at(1, 2) = 123;
  reference1[5] = 123;
  EXPECT_TRUE(is_matrix_equal(m, reference1));

  m.SetCols(1);
  EXPECT_EQ(m.GetRows(), 3);
  EXPECT_EQ(m.GetCols(), 1);

  double reference2[3] = {1, 0, 0};
  EXPECT_TRUE(is_matrix_equal(m, reference2));
}


int main(int argc, char ** argv){
    testing::InitGoogleTest( &argc, argv);
    return RUN_ALL_TESTS();
}