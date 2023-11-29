
all: clean gcov_report

install_gtest:
	git init
	git clone https://github.com/google/googletest.git
test: 
	cmake -DCMAKE_C_FLAGS="-O0 -g -fprofile-arcs -ftest-coverage" -S . -B build
	cmake --build build
	./build/s21_test

gcov_report: test
	lcov --ignore-errors mismatch -t "test" -o test.info -c -d .
	genhtml -o report test.info
	rm -rf gcovreport gcovreport.info *.gcda *.gcno math_test.gcda math_test.gcno
	open ./report/index.html
clean:
	rm -rf build/*
	rm -rf *.o test *.a a.out *.gcno *.gcda test.dSYM *.info report

