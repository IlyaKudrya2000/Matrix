#include "s21_matrix.h"

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int sing = 1;
  int flag = s21_sub_or_sum_matrix(A, B, result, sing);
  return flag;
}
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int sing = -1;
  int flag = s21_sub_or_sum_matrix(A, B, result, sing);
  return flag;
}

void s21_remove_matrix(matrix_t *A) {
  for (int i = 0; i < A->rows; i++) {
    free(A->matrix[i]);
    A->matrix[i] = NULL;
  }
  free(A->matrix);
  A->matrix = NULL;
  A->columns = 0;
  A->rows = 0;
}

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  // printf("%d\n", result->rows);
  int rez = 0;
  if (result == NULL || rows <= 0 || columns <= 0) {
    rez = 1;
    result = NULL;
  } else {
    result->rows = rows;
    result->columns = columns;
    // if(rows <= 0 || columns <= 0){
    // rez = 1;
    // } else {
    result->matrix = malloc(sizeof(double *) * result->rows);
    for (int i = 0; i < result->rows; i++) {
      result->matrix[i] = malloc(sizeof(double) * result->columns);
    }
    for (int i = 0; i < result->rows; i++) {
      for (int j = 0; j < result->columns; j++) {
        result->matrix[i][j] = 0;
      }
    }
    // }
  }
  return rez;
}
int s21_sub_or_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result,
                          int sing) {
  int flag = 0;
  if (A == NULL || B == NULL || result == NULL) {
    flag = 1;
  } else {
    if (A->columns == B->columns && A->rows == B->rows) {
      s21_create_matrix(A->rows, A->columns, result);
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] + sing * B->matrix[i][j];
        }
      }
    } else {
      flag = 2;
    }
  }
  return flag;
}

void s21_write_matrix(matrix_t dataB) {
  for (int i = 0; i < dataB.rows; i++) {
    for (int j = 0; j < dataB.columns; j++) {
      if (j != dataB.columns - 1) {
        printf("%f ", dataB.matrix[i][j]);
      } else if (i != dataB.columns - 1) {
        printf("%f\n", dataB.matrix[i][j]);
      } else {
        printf("%f\n", dataB.matrix[i][j]);
      }
    }
  }
}

void s21_write_matrix_qu(matrix_t dataB) {
  for (int i = 0; i < dataB.rows; i++) {
    for (int j = 0; j < dataB.rows; j++) {
      if (j != dataB.rows - 1) {
        printf("%f ", dataB.matrix[i][j]);
      } else if (i != dataB.rows - 1) {
        printf("%f\n", dataB.matrix[i][j]);
      } else {
        printf("%f\n", dataB.matrix[i][j]);
      }
    }
  }
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int rez = 0;
  if (A == NULL || A->columns == 0 || A->rows == 0 || result == NULL) {
    rez = 1;
  } else {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] * number;
      }
    }
  }
  return rez;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int rez = 0;
  if (A == NULL || B == NULL || result == NULL || A->columns == 0 ||
      A->rows == 0 || B->columns == 0 || B->rows == 0) {
    rez = 1;
  } else {
    if (A->columns == B->rows) {
      s21_create_matrix(A->rows, B->columns, result);
      double tmp = 0;
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < B->columns; j++) {
          for (int k = 0; k < A->columns; k++) {
            tmp = A->matrix[i][k] * B->matrix[k][j] + tmp;
          }
          result->matrix[i][j] = tmp;
          tmp = 0;
        }
      }
    } else {
      rez = 2;
    }
  }
  return rez;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int rez = 0;
  if (A == NULL || result == NULL || A->columns == 0 || A->rows == 0) {
    rez = 1;
  } else {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[j][i] = A->matrix[i][j];
      }
    }
  }
  return rez;
}

int s21_determinant(matrix_t *A, double *result) {
  int res = 0;
  *result = 0;
  if (A != NULL && result != NULL && A->columns > 0) {
    if (A->columns == A->rows) {
      // if (A->columns == 2) {
      //     *result = A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] *
      //     A->matrix[1][0];
      if (A->columns == 1 || A->rows == 1) {
        *result = A->matrix[0][0];
      } else {
        for (int j = 0; j < A->rows; j++) {
          int i = 0;
          double pp = 0;
          matrix_t ppp;
          s21_create_matrix(A->rows - 1, A->columns - 1, &ppp);
          s21_mimor(A, &ppp, i, j);
          s21_determinant(&ppp, &pp);
          // printf("%.0f\n", pp);

          s21_remove_matrix(&ppp);
          *result = *result + pow(-1, i + j) * A->matrix[i][j] * pp;
        }
      }
    } else
      res = 2;
  } else
    res = 1;
  return res;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int rez = 0;
  int sign = 0;
  if (A == NULL || result == NULL || A->columns < 2 || A->rows < 2) {
    rez = 1;
  } else {
    if (A->columns == A->rows) {
      double tmp_res = 0;
      matrix_t tmp;

      s21_create_matrix(A->columns, A->columns, result);
      s21_create_matrix(A->columns - 1, A->columns - 1, &tmp);
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          s21_mimor(A, &tmp, i, j);
          s21_determinant(&tmp, &tmp_res);
          if ((i + j) % 2 == 1) {
            sign = -1;
          } else {
            sign = 1;
          }
          result->matrix[i][j] = tmp_res * sign;
        }
      }
      s21_remove_matrix(&tmp);

    } else {
      rez = 2;
    }
  }
  return rez;
}

void s21_mimor(matrix_t *A, matrix_t *tmp, int i, int j) {
  for (int k = 0, k1 = 0; k < A->rows; k++) {
    if (k != i) {
      for (int q = 0, q1 = 0; q < A->columns; q++) {
        if (j != q) {
          tmp->matrix[k1][q1] = A->matrix[k][q];
          q1++;
        }
      }
      k1++;
    }
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int rez = 1;
  if (A == NULL || B == NULL) {
    rez = 0;
  } else {
    if (A->columns == B->columns && A->rows == B->rows) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          if (fabs(A->matrix[i][j] - B->matrix[i][j]) < pow(10, -7)) {
            //  rez = 1;
          } else {
            rez = 0;
            break;
          }
        }
      }
    } else {
      rez = 0;
    }
  }
  return rez;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int rez = 0;
  if (A == NULL || result == NULL || A->columns < 1 || A->rows < 1) {
    rez = 1;
  } else if (A->columns != A->rows) {
    rez = 2;
  } else {
    double determimant = 0;
    if (A->columns == A->rows) {
      s21_determinant(A, &determimant);
      if (determimant != 0) {
        s21_calc_complements(A, result);
        matrix_t tmp = {0};
        s21_transpose(result, &tmp);
        for (int i = 0; i < A->rows; i++) {
          for (int j = 0; j < A->columns; j++) {
            result->matrix[i][j] = tmp.matrix[i][j] / determimant;
          }
        }
        s21_remove_matrix(&tmp);
      } else {
        rez = 2;
      }
    } else {
      rez = 2;
    }
  }
  return rez;
}