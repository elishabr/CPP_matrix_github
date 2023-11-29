#include "../googletest/googletest/include/gtest/gtest.h"
#include "../googletest/googlemock/include/gmock/gmock.h"
#include "../My_s21_matrix_oop/s21_matrix_oop.hpp"

struct MyClassTest : public testing::Test {
  S21Matrix *m;

  void SetUp() { m = new S21Matrix; } // аналог конструктора
  void TearDown() { delete m; } // аналог деструктора
};

TEST_F(MyClassTest, ConstructorsDefault){
     //Arrange
     //Act
     int rows = m->GetRows();
     int cols = m->GetCols();
     //Assert
     EXPECT_EQ(rows, 3);
     EXPECT_EQ(cols, 3);
     
}

TEST(Constructors, RowsCols){
     //Arrange
     S21Matrix matrix1(2, 2);
     S21Matrix matrix2(2, 3);
     //Act
     //Assert
     EXPECT_EQ(matrix1.GetRows(), 2);
     EXPECT_EQ(matrix1.GetCols(), 2);
     EXPECT_EQ(matrix2.GetRows(), 2);
     EXPECT_EQ(matrix2.GetCols(), 3);

}

TEST_F(MyClassTest, ConstructorsCopy){
     //Arrange
     S21Matrix matrix1(*m);
     //Act
     //Assert
     EXPECT_EQ(matrix1.GetRows(), 3);
     EXPECT_EQ(matrix1.GetCols(), 3);
}

// TEST_F(MyClassTest, ConstructorsMove){
//      //Arrange
//      S21Matrix matrix1(m);
//      //Act
//      //Assert
//      EXPECT_EQ(matrix1.GetRows(), 3);
//      EXPECT_EQ(matrix1.GetCols(), 3);
// }

int main(int argc, char ** argv){
     ::testing::InitGoogleTest( &argc, argv);
     return RUN_ALL_TESTS();
}
