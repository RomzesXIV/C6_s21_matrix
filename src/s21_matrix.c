#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int err = OK;
  if (rows <= 0 || columns <= 0) {
    err = INCORRECT;
  } else if (result == NULL) {
    err = INCORRECT;
  } else {
    result->rows = rows;
    result->columns = columns;
    result->matrix = (double **)calloc(rows, sizeof(double *));
    if (result->matrix) {
      for (int i = 0; i < rows; i++) {
        result->matrix[i] = (double *)calloc(columns, sizeof(double));
      }
    }
  }
  return err;
}
void s21_remove_matrix(matrix_t *A) {
  if (A && A->rows > 0 && A->columns > 0) {
    for (int i = 0; i < A->rows; i++) {
      free(A->matrix[i]);
    }
    free(A->matrix);
  }
  A->matrix = NULL;
  A->columns = 0;
  A->rows = 0;
}
int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int err = SUCCESS;
  if (A == NULL || B == NULL) {
    err = FAILURE;
  } else if (A->rows <= 0 || A->columns <= 0 || B->rows <= 0 ||
             B->columns <= 0) {
    err = FAILURE;
  } else {
    if (A->rows == B->rows && A->columns == B->columns) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          if (fabs(A->matrix[i][j] - B->matrix[i][j]) > EPS) {
            err = FAILURE;
          }
        }
      }
    } else {
      err = FAILURE;
    }
  }
  return err;
}
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int err = check_start_add_sub(A, B);
  if (err == OK) {
    err = s21_create_matrix(A->rows, A->columns, result);
    if (err == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
        }
      }
    }
  }
  return err;
}
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int err = OK;
  err = check_start_add_sub(A, B);
  if (err == OK) {
    err = s21_create_matrix(A->rows, A->columns, result);
    if (err == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
        }
      }
    }
  }
  return err;
}
int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int err = OK;
  err = check_start_one_matrix(A);
  if (err == 0) {
    err = s21_create_matrix(A->rows, A->columns, result);
    if (err == 0) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] * number;
        }
      }
    }
  }
  return err;
}
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int err = OK;
  err = check_start_mult_matrix(A, B);
  if (err == 0) {
    err = s21_create_matrix(A->rows, B->columns, result);
    if (err == 0) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < B->columns; j++) {
          result->matrix[i][j] = 0;
          for (int z = 0; z < A->columns; z++) {
            result->matrix[i][j] =
                result->matrix[i][j] + A->matrix[i][z] * B->matrix[z][j];
          }
        }
      }
    }
  }
  return err;
}
int s21_transpose(matrix_t *A, matrix_t *result) {
  int err = check_start_one_matrix(A);
  if (err == 0) {
    err = s21_create_matrix(A->columns, A->rows, result);
    if (err == 0) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[j][i] = A->matrix[i][j];
        }
      }
    }
  }
  return err;
}
int s21_determinant(matrix_t *A, double *result) {
  int err = check_start_one_matrix(A);
  *result = 0;
  double det = 0;
  matrix_t matrix_det;
  if (result == NULL) err = INCORRECT;
  if (err == 0) {
    if (A->rows != A->columns) {
      err = CALCULATION_ERR;
    } else {
      if (A->rows == 1) {
        *result = A->matrix[0][0];
      } else if (A->rows == 2) {
        *result = A->matrix[0][0] * A->matrix[1][1] -
                  A->matrix[0][1] * A->matrix[1][0];
      } else {
        for (int j = 0; j < A->rows; j++) {
          create_minor(*A, 0, j, &matrix_det);
          err = s21_determinant(&matrix_det, &det);
          *result = *result + A->matrix[0][j] * pow(-1, j) * det;
          s21_remove_matrix(&matrix_det);
        }
      }
    }
  }
  return err;
}
int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int err = check_start_one_matrix(A);
  matrix_t matrix_det;
  double det = 0;
  if (result == NULL) err = INCORRECT;
  if (err == 0) {
    if (A->rows != A->columns) {
      err = CALCULATION_ERR;
    } else {
      err = s21_create_matrix(A->rows, A->columns, result);
      if (err == 0) {
        for (int i = 0; i < A->rows; i++) {
          for (int j = 0; j < A->columns; j++) {
            create_minor(*A, i, j, &matrix_det);
            err = s21_determinant(&matrix_det, &det);
            result->matrix[i][j] = pow(-1, i + j) * det;
            s21_remove_matrix(&matrix_det);
          }
        }
      }
    }
  }
  return err;
}
int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int err = check_start_one_matrix(A);
  double det = 0;
  matrix_t matrix_det;
  matrix_t matrix_tmp;
  if (result == NULL) err = INCORRECT;
  if (err == 0) {
    if (A->rows != A->columns) {
      err = CALCULATION_ERR;
    } else {
      err = s21_determinant(A, &det);
      if ((fabs(det) > EPS) && err == 0) {
        det = 1 / det;
        err = s21_calc_complements(A, &matrix_det);
        err = s21_transpose(&matrix_det, &matrix_tmp);
        err = s21_mult_number(&matrix_tmp, det, result);
        s21_remove_matrix(&matrix_det);
        s21_remove_matrix(&matrix_tmp);
      } else {
        err = CALCULATION_ERR;
      }
    }
  }
  return err;
}

// help
int check_start_add_sub(matrix_t *A, matrix_t *B) {
  int err = OK;
  if (A == NULL || B == NULL) {
    err = INCORRECT;
  } else if (A->rows <= 0 || A->columns <= 0 || B->rows <= 0 ||
             B->columns <= 0) {
    err = INCORRECT;
  } else if (A->rows != B->rows || A->columns != B->columns) {
    err = CALCULATION_ERR;
  }
  return err;
}
int check_start_one_matrix(matrix_t *A) {
  int err = OK;
  if (A == NULL) {
    err = INCORRECT;
  } else if (A->rows <= 0 || A->columns <= 0) {
    err = INCORRECT;
  }
  return err;
}
int check_start_mult_matrix(matrix_t *A, matrix_t *B) {
  int err = OK;
  if (A == NULL || B == NULL) {
    err = INCORRECT;
  } else if (A->rows <= 0 || A->columns <= 0 || B->rows <= 0 ||
             B->columns <= 0) {
    err = INCORRECT;
  } else if (A->rows != B->columns || A->columns != B->rows) {
    err = CALCULATION_ERR;
  }
  return err;
}
void create_minor(matrix_t A, int im, int jm, matrix_t *matrix_minor) {
  s21_create_matrix(A.rows - 1, A.columns - 1, matrix_minor);
  int index_i = 0, index_j = 0;
  for (int i = 0; i < A.rows; i++) {
    if (i != im) {
      index_j = 0;
      for (int j = 0; j < A.columns; j++) {
        if (j != jm) {
          matrix_minor->matrix[index_i][index_j] = A.matrix[i][j];
          index_j++;
        }
      }
      index_i++;
    }
  }
}
