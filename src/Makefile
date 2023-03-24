CFLAGS = -Wall -Werror -Wextra -std=c11

all: clean s21_matrix.a test gcov_report

ubuntu: clean s21_matrix.a test_ubuntu gcov_report

test: *.c
	gcc -o test $(CFLAGS) --coverage *.c -lcheck -lpthread -lm  -I /usr/local/include -L/usr/local/lib
	
test_ubuntu: *.c
	gcc -o test $(CFLAGS) --coverage *.c -lcheck -lm -lpthread -lsubunit

s21_matrix.a: s21_matrix.o
	ar rc s21_matrix.a s21_matrix.o
	ranlib s21_matrix.a

s21_matrix.o: s21_matrix.c
	gcc $(CFLAGS) -c s21_matrix.c

gcov_report: test
	./test
	gcovr --html-details rep.html

check:
	cp ../materials/linters/.clang-format ../src/.clang-format
	clang-format -i *.c
	clang-format -i *.h

clean:
	rm -rf *.o *.gcda *.a *.html *.gcno *.css s21_matrix
	rm -f test
	rm -f test_ubuntu

rebuild: clean all
