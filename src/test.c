#include <check.h>

#include "s21_matrix.h"

//----------DELETE--------------//
// void print_matrix(matrix_t matrix) {
//   for (int i = 0; i < matrix.rows; i++) {
//     for (int j = 0; j < matrix.columns; j++) {
//       printf("%lf ", matrix.matrix[i][j]);
//     }
//     printf("\n");
//   }
// }
//----------DELETE--------------//

START_TEST(create) {
  matrix_t matrix;
  int err = 0;
  err = s21_create_matrix(5, 5, &matrix);
  ck_assert_int_eq(0, err);
  s21_remove_matrix(&matrix);

  err = s21_create_matrix(5, 10, &matrix);
  ck_assert_int_eq(0, err);
  s21_remove_matrix(&matrix);

  err = s21_create_matrix(10, 5, &matrix);
  ck_assert_int_eq(0, err);
  s21_remove_matrix(&matrix);

  err = s21_create_matrix(0, 0, &matrix);
  ck_assert_int_eq(1, err);

  err = s21_create_matrix(0, 3, &matrix);
  ck_assert_int_eq(1, err);

  err = s21_create_matrix(-4, -4, &matrix);
  ck_assert_int_eq(1, err);

  err = s21_create_matrix(-4, 5, &matrix);
  ck_assert_int_eq(1, err);

  err = s21_create_matrix(4, -5, &matrix);
  ck_assert_int_eq(1, err);

  err = s21_create_matrix(2, 2, NULL);
  ck_assert_int_eq(1, err);
}
END_TEST

START_TEST(remove_m) {
  matrix_t matrix_r;
  int err = 0;

  err = s21_create_matrix(3, 3, &matrix_r);
  s21_remove_matrix(&matrix_r);
  if (matrix_r.matrix != NULL) {
    err = INCORRECT;
  }
  ck_assert_int_eq(err, OK);
}
END_TEST

START_TEST(equal) {
  matrix_t matrix_1;
  matrix_t matrix_2;
  s21_create_matrix(5, 5, &matrix_1);
  s21_create_matrix(5, 5, &matrix_2);
  int err = 0;
  err = s21_eq_matrix(&matrix_1, &matrix_2);
  ck_assert_int_eq(SUCCESS, err);
  err = s21_eq_matrix(&matrix_1, NULL);
  ck_assert_int_eq(FAILURE, err);

  matrix_1.matrix[0][0] = 5.0;
  matrix_1.matrix[0][1] = 4.0;
  matrix_1.matrix[0][2] = 3.0;
  matrix_1.matrix[0][3] = 2.0;
  matrix_1.matrix[0][4] = 1.0;
  matrix_1.matrix[1][0] = 1.0;
  matrix_1.matrix[1][1] = 2.0;
  matrix_1.matrix[1][2] = 3.0;
  matrix_1.matrix[1][3] = 4.0;
  matrix_1.matrix[1][4] = 5.0;
  matrix_1.matrix[2][0] = 2.0;
  matrix_1.matrix[2][1] = 1.0;
  matrix_1.matrix[2][2] = 3.0;
  matrix_1.matrix[2][3] = 3.0;
  matrix_1.matrix[2][4] = 5.0;
  matrix_1.matrix[3][0] = 3.0;
  matrix_1.matrix[3][1] = 2.0;
  matrix_1.matrix[3][2] = 4.0;
  matrix_1.matrix[3][3] = 1.0;
  matrix_1.matrix[3][4] = 5.0;
  matrix_1.matrix[4][0] = 5.0;
  matrix_1.matrix[4][1] = 3.0;
  matrix_1.matrix[4][2] = 1.0;
  matrix_1.matrix[4][3] = 2.0;
  matrix_1.matrix[4][4] = 4.0;

  matrix_2.matrix[0][0] = 5.0;
  matrix_2.matrix[0][1] = 4.0;
  matrix_2.matrix[0][2] = 3.0;
  matrix_2.matrix[0][3] = 2.0;
  matrix_2.matrix[0][4] = 1.0;
  matrix_2.matrix[1][0] = 1.0;
  matrix_2.matrix[1][1] = 2.0;
  matrix_2.matrix[1][2] = 3.0;
  matrix_2.matrix[1][3] = 4.0;
  matrix_2.matrix[1][4] = 5.0;
  matrix_2.matrix[2][0] = 2.0;
  matrix_2.matrix[2][1] = 1.0;
  matrix_2.matrix[2][2] = 3.0;
  matrix_2.matrix[2][3] = 3.0;
  matrix_2.matrix[2][4] = 5.0;
  matrix_2.matrix[3][0] = 3.0;
  matrix_2.matrix[3][1] = 2.0;
  matrix_2.matrix[3][2] = 4.0;
  matrix_2.matrix[3][3] = 1.0;
  matrix_2.matrix[3][4] = 5.0;
  matrix_2.matrix[4][0] = 5.0;
  matrix_2.matrix[4][1] = 3.0;
  matrix_2.matrix[4][2] = 1.0;
  matrix_2.matrix[4][3] = 2.0;
  matrix_2.matrix[4][4] = 4.0;

  err = s21_eq_matrix(&matrix_1, &matrix_2);
  ck_assert_int_eq(SUCCESS, err);
  matrix_1.matrix[1][1] = 8.0;
  err = s21_eq_matrix(&matrix_1, &matrix_2);
  ck_assert_int_eq(FAILURE, err);

  matrix_1.matrix[0][0] = 5;
  matrix_1.matrix[0][1] = 4;
  matrix_1.matrix[0][2] = 3;
  matrix_1.matrix[0][3] = 2;
  matrix_1.matrix[0][4] = 1;
  matrix_1.matrix[1][0] = 1;
  matrix_1.matrix[1][1] = 2;
  matrix_1.matrix[1][2] = 3;
  matrix_1.matrix[1][3] = 4;
  matrix_1.matrix[1][4] = 5;
  matrix_1.matrix[2][0] = 2;
  matrix_1.matrix[2][1] = 1;
  matrix_1.matrix[2][2] = 3;
  matrix_1.matrix[2][3] = 3;
  matrix_1.matrix[2][4] = 5;
  matrix_1.matrix[3][0] = 3;
  matrix_1.matrix[3][1] = 2;
  matrix_1.matrix[3][2] = 4;
  matrix_1.matrix[3][3] = 1;
  matrix_1.matrix[3][4] = 5;
  matrix_1.matrix[4][0] = 5;
  matrix_1.matrix[4][1] = 3;
  matrix_1.matrix[4][2] = 1;
  matrix_1.matrix[4][3] = 2;
  matrix_1.matrix[4][4] = 4;

  matrix_2.matrix[0][0] = 5;
  matrix_2.matrix[0][1] = 4;
  matrix_2.matrix[0][2] = 3;
  matrix_2.matrix[0][3] = 2;
  matrix_2.matrix[0][4] = 1;
  matrix_2.matrix[1][0] = 1;
  matrix_2.matrix[1][1] = 2;
  matrix_2.matrix[1][2] = 3;
  matrix_2.matrix[1][3] = 4;
  matrix_2.matrix[1][4] = 5;
  matrix_2.matrix[2][0] = 2;
  matrix_2.matrix[2][1] = 1;
  matrix_2.matrix[2][2] = 3;
  matrix_2.matrix[2][3] = 3;
  matrix_2.matrix[2][4] = 5;
  matrix_2.matrix[3][0] = 3;
  matrix_2.matrix[3][1] = 2;
  matrix_2.matrix[3][2] = 4;
  matrix_2.matrix[3][3] = 1;
  matrix_2.matrix[3][4] = 5;
  matrix_2.matrix[4][0] = 5;
  matrix_2.matrix[4][1] = 3;
  matrix_2.matrix[4][2] = 1;
  matrix_2.matrix[4][3] = 2;
  matrix_2.matrix[4][4] = 4;

  matrix_2.matrix[1][1] = 5;
  err = s21_eq_matrix(&matrix_1, &matrix_2);
  ck_assert_int_eq(FAILURE, err);
  s21_remove_matrix(&matrix_2);

  s21_create_matrix(5, 5, &matrix_2);
  err = s21_eq_matrix(&matrix_1, &matrix_2);
  ck_assert_int_eq(FAILURE, err);

  matrix_1.columns = 0;
  err = s21_eq_matrix(&matrix_1, &matrix_2);
  ck_assert_int_eq(FAILURE, err);

  s21_remove_matrix(&matrix_1);
  s21_remove_matrix(&matrix_2);

  s21_create_matrix(2, 2, &matrix_1);
  s21_create_matrix(2, 2, &matrix_2);

  err = s21_eq_matrix(&matrix_1, &matrix_2);
  ck_assert_int_eq(SUCCESS, err);

  err = s21_eq_matrix(&matrix_1, NULL);
  ck_assert_int_eq(FAILURE, err);

  matrix_1.matrix[0][0] = 5.00000011;
  matrix_1.matrix[0][1] = 4.0;
  matrix_1.matrix[1][1] = 3.0;
  matrix_1.matrix[1][2] = 2.0;

  matrix_2.matrix[0][0] = 5.000000012;
  matrix_2.matrix[0][1] = 4.0;
  matrix_2.matrix[1][1] = 3.0;
  matrix_2.matrix[1][2] = 2.0;

  err = s21_eq_matrix(&matrix_1, &matrix_2);
  ck_assert_int_eq(SUCCESS, err);

  s21_remove_matrix(&matrix_1);
  s21_remove_matrix(&matrix_2);
}
END_TEST

START_TEST(sum) {
  matrix_t matrix_1, matrix_2, result, correct_res;
  s21_create_matrix(2, 2, &matrix_1);
  s21_create_matrix(2, 2, &matrix_2);
  s21_create_matrix(2, 2, &correct_res);
  int err = 0;
  err = s21_sum_matrix(&matrix_1, &matrix_2, NULL);
  ck_assert_int_eq(INCORRECT, err);
  err = s21_sum_matrix(&matrix_1, NULL, &result);
  ck_assert_int_eq(INCORRECT, err);
  err = s21_sum_matrix(NULL, &matrix_2, &result);
  ck_assert_int_eq(INCORRECT, err);
  err = s21_sum_matrix(&matrix_1, &matrix_2, &result);
  ck_assert_int_eq(OK, err);

  matrix_1.matrix[0][0] = 1.2;
  matrix_1.matrix[0][1] = 2;
  matrix_1.matrix[1][0] = 3;
  matrix_1.matrix[1][1] = 4;
  matrix_2.matrix[0][0] = 1;
  matrix_2.matrix[0][1] = 2;
  matrix_2.matrix[1][0] = 3;
  matrix_2.matrix[1][1] = 4;
  correct_res.matrix[0][0] = 2.2;
  correct_res.matrix[0][1] = 4;
  correct_res.matrix[1][0] = 6;
  correct_res.matrix[1][1] = 8;

  err = s21_sum_matrix(&matrix_1, &matrix_2, &result);
  ck_assert_int_eq(OK, err);
  err = s21_eq_matrix(&result, &correct_res);
  ck_assert_int_eq(SUCCESS, err);

  matrix_1.matrix[0][0] = 1;
  matrix_1.matrix[0][1] = 2;
  matrix_1.matrix[1][0] = 3;
  matrix_1.matrix[1][1] = 4;
  matrix_2.matrix[0][0] = 1;
  matrix_2.matrix[0][1] = 2;
  matrix_2.matrix[1][0] = 3;
  matrix_2.matrix[1][1] = 4;
  correct_res.matrix[0][0] = 2;
  correct_res.matrix[0][1] = 4;
  correct_res.matrix[1][0] = 22;
  correct_res.matrix[1][1] = 8;

  err = s21_sum_matrix(&matrix_1, &matrix_2, &result);
  ck_assert_int_eq(OK, err);
  err = s21_eq_matrix(&result, &correct_res);
  ck_assert_int_eq(FAILURE, err);

  matrix_1.columns = -10;
  err = s21_sum_matrix(&matrix_1, &matrix_2, &result);
  ck_assert_int_eq(INCORRECT, err);

  s21_create_matrix(2, 6, &matrix_1);
  err = s21_sum_matrix(&matrix_1, &matrix_2, &result);
  ck_assert_int_eq(CALCULATION_ERR, err);

  s21_remove_matrix(&matrix_1);
  s21_remove_matrix(&matrix_2);
  s21_remove_matrix(&result);
  s21_remove_matrix(&correct_res);
}
END_TEST

START_TEST(sub) {
  matrix_t matrix_1, matrix_2, result, correct_res;
  s21_create_matrix(2, 2, &matrix_1);
  s21_create_matrix(2, 2, &matrix_2);
  s21_create_matrix(2, 2, &correct_res);
  int err = 0;
  err = s21_sum_matrix(&matrix_1, &matrix_2, NULL);
  ck_assert_int_eq(INCORRECT, err);
  err = s21_sum_matrix(&matrix_1, NULL, &result);
  ck_assert_int_eq(INCORRECT, err);
  err = s21_sum_matrix(NULL, &matrix_2, &result);
  ck_assert_int_eq(INCORRECT, err);
  err = s21_sum_matrix(&matrix_1, &matrix_2, &result);
  ck_assert_int_eq(OK, err);

  matrix_1.matrix[0][0] = 2.2;
  matrix_1.matrix[0][1] = 10;
  matrix_1.matrix[1][0] = 7;
  matrix_1.matrix[1][1] = 14;
  matrix_2.matrix[0][0] = 1;
  matrix_2.matrix[0][1] = 2;
  matrix_2.matrix[1][0] = 3;
  matrix_2.matrix[1][1] = 4;
  correct_res.matrix[0][0] = 1.2;
  correct_res.matrix[0][1] = 8;
  correct_res.matrix[1][0] = 4;
  correct_res.matrix[1][1] = 10;

  err = s21_sub_matrix(&matrix_1, &matrix_2, &result);
  ck_assert_int_eq(OK, err);
  err = s21_eq_matrix(&result, &correct_res);
  ck_assert_int_eq(SUCCESS, err);

  matrix_1.matrix[0][0] = 2;
  matrix_1.matrix[0][1] = 4;
  matrix_1.matrix[1][0] = 6;
  matrix_1.matrix[1][1] = 8;
  matrix_2.matrix[0][0] = 1;
  matrix_2.matrix[0][1] = 2;
  matrix_2.matrix[1][0] = 3;
  matrix_2.matrix[1][1] = 4;
  correct_res.matrix[0][0] = 1;
  correct_res.matrix[0][1] = 2;
  correct_res.matrix[1][0] = 22;
  correct_res.matrix[1][1] = 4;

  err = s21_sub_matrix(&matrix_1, &matrix_2, &result);
  ck_assert_int_eq(OK, err);
  err = s21_eq_matrix(&result, &correct_res);
  ck_assert_int_eq(FAILURE, err);

  matrix_1.columns = -10;
  err = s21_sub_matrix(&matrix_1, &matrix_2, &result);
  ck_assert_int_eq(INCORRECT, err);

  s21_create_matrix(2, 6, &matrix_1);
  err = s21_sub_matrix(&matrix_1, &matrix_2, &result);
  ck_assert_int_eq(CALCULATION_ERR, err);

  s21_remove_matrix(&matrix_1);
  s21_remove_matrix(&matrix_2);
  s21_remove_matrix(&result);
  s21_remove_matrix(&correct_res);
}
END_TEST

START_TEST(mult_number) {
  matrix_t matrix_1, result, correct_res;
  s21_create_matrix(2, 2, &matrix_1);
  s21_create_matrix(2, 2, &correct_res);
  int err = 0;

  err = s21_mult_number(NULL, 0, &result);
  ck_assert_int_eq(INCORRECT, err);
  err = s21_mult_number(&matrix_1, 0, NULL);
  ck_assert_int_eq(INCORRECT, err);

  matrix_1.matrix[0][0] = 1.1;
  matrix_1.matrix[0][1] = 2.0;
  matrix_1.matrix[1][0] = 3.0;
  matrix_1.matrix[1][1] = 4.0;
  correct_res.matrix[0][0] = 2.2;
  correct_res.matrix[0][1] = 4.0;
  correct_res.matrix[1][0] = 6.0;
  correct_res.matrix[1][1] = 8.0;

  err = s21_mult_number(&matrix_1, 2, &result);
  ck_assert_int_eq(OK, err);
  err = s21_eq_matrix(&result, &correct_res);
  ck_assert_int_eq(SUCCESS, err);

  correct_res.matrix[0][0] = 2.2;
  correct_res.matrix[0][1] = 4.0;
  correct_res.matrix[1][0] = 50.0;
  correct_res.matrix[1][1] = 8.0;

  err = s21_mult_number(&matrix_1, 2, &result);
  ck_assert_int_eq(OK, err);
  err = s21_eq_matrix(&result, &correct_res);
  ck_assert_int_eq(FAILURE, err);

  matrix_1.columns = -5;
  err = s21_mult_number(&matrix_1, 2, &result);
  ck_assert_int_eq(INCORRECT, err);

  s21_remove_matrix(&matrix_1);
  s21_remove_matrix(&result);
  s21_remove_matrix(&correct_res);
}
END_TEST

START_TEST(mult_matrix) {
  matrix_t matrix_1, matrix_2, result, correct_res;
  s21_create_matrix(3, 1, &matrix_1);
  s21_create_matrix(1, 3, &matrix_2);
  s21_create_matrix(3, 3, &correct_res);
  int err = 0;

  err = s21_mult_matrix(NULL, &matrix_2, &result);
  ck_assert_int_eq(INCORRECT, err);
  err = s21_mult_matrix(&matrix_1, NULL, &result);
  ck_assert_int_eq(INCORRECT, err);
  err = s21_mult_matrix(&matrix_1, &matrix_2, NULL);
  ck_assert_int_eq(INCORRECT, err);

  matrix_1.matrix[0][0] = 2;
  matrix_1.matrix[1][0] = 5;
  matrix_1.matrix[2][0] = 3;
  matrix_2.matrix[0][0] = 2;
  matrix_2.matrix[0][1] = 5;
  matrix_2.matrix[0][2] = 3;
  correct_res.matrix[0][0] = 4.0;
  correct_res.matrix[0][1] = 10.0;
  correct_res.matrix[0][2] = 6.0;
  correct_res.matrix[1][0] = 10;
  correct_res.matrix[1][1] = 25;
  correct_res.matrix[1][2] = 15;
  correct_res.matrix[2][0] = 6;
  correct_res.matrix[2][1] = 15;
  correct_res.matrix[2][2] = 9;

  err = s21_mult_matrix(&matrix_1, &matrix_2, &result);
  ck_assert_int_eq(OK, err);
  err = s21_eq_matrix(&result, &correct_res);
  ck_assert_int_eq(SUCCESS, err);

  correct_res.matrix[2][2] = 100;

  err = s21_mult_matrix(&matrix_1, &matrix_2, &result);
  ck_assert_int_eq(OK, err);
  err = s21_eq_matrix(&result, &correct_res);
  ck_assert_int_eq(FAILURE, err);

  matrix_2.rows = -1;
  err = s21_mult_matrix(&matrix_1, &matrix_2, &result);
  ck_assert_int_eq(INCORRECT, err);

  matrix_2.rows = 5;
  err = s21_mult_matrix(&matrix_1, &matrix_2, &result);
  ck_assert_int_eq(CALCULATION_ERR, err);

  s21_create_matrix(4, 4, &matrix_2);
  err = s21_mult_matrix(&matrix_1, &matrix_2, &result);
  ck_assert_int_eq(CALCULATION_ERR, err);

  s21_remove_matrix(&matrix_1);
  s21_remove_matrix(&matrix_2);
  s21_remove_matrix(&result);
  s21_remove_matrix(&correct_res);
}
END_TEST

START_TEST(transpose) {
  matrix_t matrix_1, result, correct_res;
  s21_create_matrix(3, 2, &matrix_1);
  s21_create_matrix(2, 3, &correct_res);
  int err = 0;

  err = s21_transpose(NULL, &result);
  ck_assert_int_eq(INCORRECT, err);
  err = s21_transpose(&matrix_1, NULL);
  ck_assert_int_eq(INCORRECT, err);

  matrix_1.matrix[0][0] = 1;
  matrix_1.matrix[0][1] = 4;
  matrix_1.matrix[1][0] = 2;
  matrix_1.matrix[1][1] = 5;
  matrix_1.matrix[2][0] = 3;
  matrix_1.matrix[2][1] = 6;
  correct_res.matrix[0][0] = 1;
  correct_res.matrix[0][1] = 2;
  correct_res.matrix[0][2] = 3;
  correct_res.matrix[1][0] = 4;
  correct_res.matrix[1][1] = 5;
  correct_res.matrix[1][2] = 6;

  err = s21_transpose(&matrix_1, &result);
  ck_assert_int_eq(OK, err);
  err = s21_eq_matrix(&result, &correct_res);
  ck_assert_int_eq(SUCCESS, err);

  correct_res.matrix[1][2] = 10;
  err = s21_transpose(&matrix_1, &result);
  ck_assert_int_eq(OK, err);
  err = s21_eq_matrix(&result, &correct_res);
  ck_assert_int_eq(FAILURE, err);

  matrix_1.columns = 0;
  err = s21_transpose(&matrix_1, &result);
  ck_assert_int_eq(INCORRECT, err);

  s21_remove_matrix(&matrix_1);
  s21_remove_matrix(&result);
  s21_remove_matrix(&correct_res);
}
END_TEST

START_TEST(determinant) {
  matrix_t matrix_1;
  s21_create_matrix(3, 2, &matrix_1);
  double result;
  int err = 0;

  err = s21_determinant(NULL, &result);
  ck_assert_int_eq(INCORRECT, err);
  err = s21_determinant(&matrix_1, &result);
  ck_assert_int_eq(CALCULATION_ERR, err);

  s21_create_matrix(3, 3, &matrix_1);

  matrix_1.matrix[0][0] = 2;
  matrix_1.matrix[1][0] = 5;
  matrix_1.matrix[2][0] = 3;
  matrix_1.matrix[0][1] = 1;
  matrix_1.matrix[1][1] = -5;
  matrix_1.matrix[2][1] = 3.3;
  matrix_1.matrix[0][2] = -1;
  matrix_1.matrix[1][2] = 0;
  matrix_1.matrix[2][2] = 0.3;

  err = s21_determinant(&matrix_1, &result);
  ck_assert_int_eq(OK, err);
  ck_assert_int_eq(1, result == -36);

  s21_create_matrix(2, 2, &matrix_1);
  matrix_1.matrix[0][0] = 2;
  matrix_1.matrix[1][0] = 5;
  matrix_1.matrix[0][1] = 3;
  matrix_1.matrix[1][1] = 1;

  err = s21_determinant(&matrix_1, &result);
  ck_assert_int_eq(OK, err);
  ck_assert_int_eq(1, result == -13);

  matrix_1.rows = -1;
  err = s21_determinant(&matrix_1, &result);
  ck_assert_int_eq(INCORRECT, err);

  s21_remove_matrix(&matrix_1);
}
END_TEST

START_TEST(calc_complements) {
  matrix_t matrix_1, result, correct_res;
  s21_create_matrix(3, 2, &matrix_1);
  s21_create_matrix(3, 3, &correct_res);
  int err = 0;

  err = s21_calc_complements(NULL, &result);
  ck_assert_int_eq(INCORRECT, err);
  err = s21_calc_complements(&matrix_1, NULL);
  ck_assert_int_eq(INCORRECT, err);
  err = s21_calc_complements(&matrix_1, &result);
  ck_assert_int_eq(CALCULATION_ERR, err);

  s21_create_matrix(3, 3, &matrix_1);

  matrix_1.matrix[0][0] = 2;
  matrix_1.matrix[1][0] = 5;
  matrix_1.matrix[2][0] = 3;
  matrix_1.matrix[0][1] = 1;
  matrix_1.matrix[1][1] = -5;
  matrix_1.matrix[2][1] = 3.3;
  matrix_1.matrix[0][2] = -1;
  matrix_1.matrix[1][2] = 0;
  matrix_1.matrix[2][2] = 0.3;
  correct_res.matrix[0][0] = -1.5;
  correct_res.matrix[0][1] = -1.5;
  correct_res.matrix[0][2] = 31.5;
  correct_res.matrix[1][0] = -3.6;
  correct_res.matrix[1][1] = 3.6;
  correct_res.matrix[1][2] = -3.6;
  correct_res.matrix[2][0] = -5;
  correct_res.matrix[2][1] = -5;
  correct_res.matrix[2][2] = -15;

  err = s21_calc_complements(&matrix_1, &result);
  ck_assert_int_eq(OK, err);
  err = s21_eq_matrix(&result, &correct_res);
  ck_assert_int_eq(SUCCESS, err);

  correct_res.matrix[1][2] = 100;
  err = s21_calc_complements(&matrix_1, &result);
  ck_assert_int_eq(OK, err);
  err = s21_eq_matrix(&result, &correct_res);
  ck_assert_int_eq(FAILURE, err);

  matrix_1.rows = -3;
  err = s21_calc_complements(&matrix_1, &result);
  ck_assert_int_eq(INCORRECT, err);

  s21_remove_matrix(&matrix_1);
  s21_remove_matrix(&result);
  s21_remove_matrix(&correct_res);
}
END_TEST

START_TEST(inverse_matrix) {
  matrix_t matrix_1, result, correct_res;
  s21_create_matrix(3, 2, &matrix_1);
  s21_create_matrix(3, 3, &correct_res);
  int err = 0;

  err = s21_inverse_matrix(NULL, &result);
  ck_assert_int_eq(INCORRECT, err);
  err = s21_inverse_matrix(&matrix_1, NULL);
  ck_assert_int_eq(INCORRECT, err);
  err = s21_inverse_matrix(&matrix_1, &result);
  ck_assert_int_eq(CALCULATION_ERR, err);

  s21_create_matrix(3, 3, &matrix_1);

  matrix_1.matrix[0][0] = 2;
  matrix_1.matrix[1][0] = 5;
  matrix_1.matrix[2][0] = 3;
  matrix_1.matrix[0][1] = 1;
  matrix_1.matrix[1][1] = -5;
  matrix_1.matrix[2][1] = 3.3;
  matrix_1.matrix[0][2] = -1;
  matrix_1.matrix[1][2] = 0;
  matrix_1.matrix[2][2] = 0.3;
  correct_res.matrix[0][0] = -1.5 / -36;
  correct_res.matrix[1][0] = -1.5 / -36;
  correct_res.matrix[2][0] = 31.5 / -36;
  correct_res.matrix[0][1] = -3.6 / -36;
  correct_res.matrix[1][1] = 3.6 / -36;
  correct_res.matrix[2][1] = -3.6 / -36;
  correct_res.matrix[0][2] = -5. / -36;
  correct_res.matrix[1][2] = -5. / -36;
  correct_res.matrix[2][2] = -15. / -36;

  err = s21_inverse_matrix(&matrix_1, &result);
  ck_assert_int_eq(OK, err);
  err = s21_eq_matrix(&result, &correct_res);
  ck_assert_int_eq(SUCCESS, err);

  correct_res.matrix[2][2] = -14. / -36;
  err = s21_inverse_matrix(&matrix_1, &result);
  ck_assert_int_eq(OK, err);
  err = s21_eq_matrix(&result, &correct_res);
  ck_assert_int_eq(FAILURE, err);

  matrix_1.matrix[0][2] = 0;
  matrix_1.matrix[1][2] = 0;
  matrix_1.matrix[2][2] = 0;
  err = s21_inverse_matrix(&matrix_1, &result);
  ck_assert_int_eq(CALCULATION_ERR, err);

  matrix_1.rows = -2;
  err = s21_inverse_matrix(&matrix_1, &result);
  ck_assert_int_eq(INCORRECT, err);

  s21_remove_matrix(&matrix_1);
  s21_remove_matrix(&result);
  s21_remove_matrix(&correct_res);
}
END_TEST

Suite *s21_matrix_tests_create() {
  Suite *s21_matrix;
  TCase *s21_matrix_tests;
  s21_matrix = suite_create("s21_matrix");
  s21_matrix_tests = tcase_create("tests");
  tcase_add_test(s21_matrix_tests, create);
  tcase_add_test(s21_matrix_tests, remove_m);
  tcase_add_test(s21_matrix_tests, equal);
  tcase_add_test(s21_matrix_tests, sum);
  tcase_add_test(s21_matrix_tests, sub);
  tcase_add_test(s21_matrix_tests, mult_number);
  tcase_add_test(s21_matrix_tests, mult_matrix);
  tcase_add_test(s21_matrix_tests, transpose);
  tcase_add_test(s21_matrix_tests, determinant);
  tcase_add_test(s21_matrix_tests, calc_complements);
  tcase_add_test(s21_matrix_tests, inverse_matrix);
  suite_add_tcase(s21_matrix, s21_matrix_tests);
  return s21_matrix;
}
int main() {
  Suite *s21_matrix = s21_matrix_tests_create();
  SRunner *s21_matrix_runner = srunner_create(s21_matrix);
  int number_failed;
  srunner_run_all(s21_matrix_runner, CK_NORMAL);
  number_failed = srunner_ntests_failed(s21_matrix_runner);
  srunner_free(s21_matrix_runner);

  return number_failed == 0 ? 0 : 1;
}
