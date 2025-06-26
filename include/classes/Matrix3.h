#ifndef MATRIX3_H
#define MATRIX3_H

#include <ostream>
#include <stdexcept>

#include "Vector3.h"

/**
 * @brief Класс для хранения матрицы 3x3.
 *
 * Класс для хранения матрицы 3x3.
 */
class Matrix3 {
private:
  /**
   * @brief Матрица 3x3.
   *
   * Матрица 3x3.
   */
  double m[3][3];

public:
  /**
   * @brief Конструктор по умолчанию.
   *
   * Конструктор по умолчанию. Все элементы матрицы инициализируются нулями.
   */
  Matrix3() {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        m[i][j] = 0.0;
  }

  /**
   * @brief Конструктор с инициализацией всех элементов матрицы.
   *
   * Конструктор с инициализацией всех элементов матрицы.
   * @param init_value - значение для инициализации всех элементов матрицы.
   */
  Matrix3(double init_value) {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        m[i][j] = init_value;
  }

  
  inline double &xx() { return m[0][0]; }
  inline double &xy() { return m[0][1]; }
  inline double &xz() { return m[0][2]; }
  inline double &yx() { return m[1][0]; }
  inline double &yy() { return m[1][1]; }
  inline double &yz() { return m[1][2]; }
  inline double &zx() { return m[2][0]; }
  inline double &zy() { return m[2][1]; }
  inline double &zz() { return m[2][2]; }

  inline const double &xx() const { return m[0][0]; }
  inline const double &xy() const { return m[0][1]; }
  inline const double &xz() const { return m[0][2]; }
  inline const double &yx() const { return m[1][0]; }
  inline const double &yy() const { return m[1][1]; }
  inline const double &yz() const { return m[1][2]; }
  inline const double &zx() const { return m[2][0]; }
  inline const double &zy() const { return m[2][1]; }
  inline const double &zz() const { return m[2][2]; }

  /**
   * @brief Перегрузка оператора () для доступа к элементам матрицы.
   *
   * Перегрузка оператора () для доступа к элементам матрицы.
   * @param row - индекс строки.
   * @param col - индекс столбца.
   * @return значение элемента матрицы.
   */
  inline double &operator()(int row, int col) {
    if (row < 0 || row > 2 || col < 0 || col > 2)
      throw std::out_of_range("Matrix3 indices out of range");
    return m[row][col];
  }

  /**
   * @brief Перегрузка оператора () для доступа к элементам матрицы.
   *
   * Перегрузка оператора () для доступа к элементам матрицы.
   * @param row - индекс строки.
   * @param col - индекс столбца.
   * @return значение элемента матрицы.
   */
  inline const double &operator()(int row, int col) const {
    if (row < 0 || row > 2 || col < 0 || col > 2)
      throw std::out_of_range("Matrix3 indices out of range");
    return m[row][col];
  }

  /**
   * @brief Перегрузка оператора + для сложения двух матриц.
   *
   * Перегрузка оператора + для сложения двух матриц.
   * @param other - матрица для сложения.
   * @return результат сложения двух матриц.
   */
  inline Matrix3 operator+(const Matrix3 &other) const {
    Matrix3 result;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        result(i, j) = m[i][j] + other(i, j);
    return result;
  }

  /**
   * @brief Перегрузка оператора += для сложения двух матриц.
   *
   * Перегрузка оператора += для сложения двух матриц.
   * @param other - матрица для сложения.
   * @return результат сложения двух матриц.
   */
  inline Matrix3 &operator+=(const Matrix3 &other) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        m[i][j] += other.m[i][j];
      }
    }
    return *this;
  }

  /**
   * @brief Перегрузка оператора * для умножения матрицы на скаляр.
   *
   * Перегрузка оператора * для умножения матрицы на скаляр.
   * @param scalar - скаляр для умножения.
   * @return результат умножения матрицы на скаляр.
   */
  inline Matrix3 operator*(const double scalar) const {
    Matrix3 result;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        result(i, j) = m[i][j] * scalar;
    return result;
  }

  /**
   * @brief Перегрузка оператора * для умножения скаляра на матрицу.
   *
   * Перегрузка оператора * для умножения скаляра на матрицу.
   * @param scalar - скаляр для умножения.
   * @return результат умножения матрицы на скаляр.
   */
  friend Matrix3 operator*(double scalar, const Matrix3 &matrix) {
    return matrix * scalar;
  }

  /**
   * @brief Перегрузка оператора / для деления матрицы на скаляр.
   *
   * Перегрузка оператора / для деления матрицы на скаляр.
   * @param scalar - скаляр для деления.
   * @return результат деления матрицы на скаляр.
   */
  inline Matrix3 operator/(const double scalar) const {
    Matrix3 result;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        result(i, j) = m[i][j] / scalar;
    return result;
  }

  /**
   * @brief Перегрузка оператора / для деления скаляра на матрицу.
   *
   * Перегрузка оператора / для деления скаляра на матрицу.
   * @param scalar - скаляр для деления.
   * @return результат деления матрицы на скаляр.
   */
  friend Matrix3 operator/(double scalar, const Matrix3 &matrix) {
    return matrix / scalar;
  }

  /**
   * @brief Умножение вектора на вектор.
   *
   * Умножение вектора на вектор.
   * @param v1 - вектор для умножения.
   * @param v2 - вектор для умножения.
   * @return результат умножения вектора на вектор.
   */
  static inline Matrix3 outerProduct(const Vector3<double> &v1,
                                     const Vector3<double> &v2) {
    Matrix3 result;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        result.m[i][j] = v1[i] * v2[j];
      }
    return result;
  }

  /**
   * @brief Перегрузка оператора << для вывода матрицы.
   *
   * Перегрузка оператора << для вывода матрицы.
   */
  friend std::ostream &operator<<(std::ostream &os, const Matrix3 &mat) {
    os << "[ " << mat.xx() << " " << mat.xy() << " " << mat.xz() << " ],"
       << "[ " << mat.yx() << " " << mat.yy() << " " << mat.yz() << " ],"
       << "[ " << mat.zx() << " " << mat.zy() << " " << mat.zz() << " ]";
    return os;
  }
};

#endif
