# s21_matrix

## Contents

1. [Структура матрицы на языке C](#структура-матрицы-на-языке-c) \
2. [Операции над матрицами](#операции-над-матрицами) \
3. [Реализация функции библиотеки matrix.h](#реализация-функциий-библиотеки-matrix) \

## Структура матрицы на языке C

```c
typedef struct matrix_struct {
    double** matrix;
    int rows;
    int columns;
} matrix_t;
```

## Операции над матрицами

Все операции (кроме сравнения матриц) возвращают результирующий код:  
- 0 - OK
- 1 - Ошибка, некорректная матрица   
- 2 - Ошибка вычисления (несовпадающие размеры матриц; матрица, для которой нельзя провести вычисления и т.д.)

## Реализация функциий библиотеки matrix

Реализованы основные действия с матрицами : s21_create_matrix (создание), s21_remove_matrix (очистка и уничтожение), s21_eq_matrix (сравнение), s21_sum_matrix (сложение), s21_sub_matrix (вычитание), s21_mult_matrix (умножение), s21_mult_number (умножение на число), s21_transpose (транспонирование), s21_determinant (вычисление определителя), s21_calc_complements (вычисление матрицы алгебраических дополнений), s21_inverse_matrix (поиск обратной матрицы). 

- Библиотека разработана на языке Си стандарта C11 с использованием компилятора gcc  
- Решение оформлено как статическая библиотека (с заголовочным файлом s21_matrix.h)
- Библиотека разработана в соответствии с принципами структурного программирования
- Подготовлено полное покрытие unit-тестами функций библиотеки c помощью библиотеки Check
- Предусмотрен Makefile для сборки библиотеки и тестов (с целями all, clean, test, s21_matrix.a, gcov_report)
- В цели gcov_report формируется отчёт gcov в виде html страницы. Для этого unit-тесты должны запускаться с флагами gcov 
- Проверяемая точность дробной части - максимум 6 знаков после запятой.

