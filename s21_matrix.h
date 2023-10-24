#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define SUCCESS 1
#define FAILURE 0

typedef struct matrix_struct {
  double **matrix;
  int rows;     // количество рядов size y rows
  int columns;  // количество столбцов size x columns
} matrix_t;

void s21_write_matrix(matrix_t dataB);
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_or_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result, int sing);
int s21_create_matrix(int rows, int columns, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_transpose(matrix_t *A, matrix_t *result);
int s21_determinant(matrix_t *A, double *result);
void s21_write_matrix_qu(matrix_t dataB);
int s21_calc_complements(matrix_t *A, matrix_t *result);
void s21_remove_matrix(matrix_t *A);
void s21_mimor(matrix_t *A, matrix_t *tmp, int i, int j);
int s21_eq_matrix(matrix_t *A, matrix_t *B);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
