#ifndef s21_matrix_h
#define s21_matrix_h

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define OK 0
#define INCORRECT 1
#define CALCULATION_ERR 2

#define SUCCESS 1
#define FAILURE 0
#define EPS 1e-7

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

int s21_create_matrix(int rows, int columns, matrix_t *result);
void s21_remove_matrix(matrix_t *A);
int s21_eq_matrix(matrix_t *A, matrix_t *B);
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_transpose(matrix_t *A, matrix_t *result);
int s21_determinant(matrix_t *A, double *result);
int s21_calc_complements(matrix_t *A, matrix_t *result);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);

int check_start_add_sub(matrix_t *A, matrix_t *B);
int check_start_one_matrix(matrix_t *A);
int check_start_mult_matrix(matrix_t *A, matrix_t *B);
void create_minor(matrix_t A, int im, int jm, matrix_t *matrix_minor);

#endif /* s21_matrix_h */
