#ifndef MATRIX3_H
#define MATRIX3_H

class Matrix3;
#include "Vector3.h"

class Matrix3 {
private:
  double m[3][3];

public:
  // Default constructor: initialize all elements to zero.
  Matrix3() {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        m[i][j] = 0.0;
  }
  Matrix3(double init_value) {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        m[i][j] = init_value;
  }

  // Named accessors (non-const)
  inline double &xx() { return m[0][0]; }
  inline double &xy() { return m[0][1]; }
  inline double &xz() { return m[0][2]; }
  inline double &yx() { return m[1][0]; }
  inline double &yy() { return m[1][1]; }
  inline double &yz() { return m[1][2]; }
  inline double &zx() { return m[2][0]; }
  inline double &zy() { return m[2][1]; }
  inline double &zz() { return m[2][2]; }

  // Named accessors (const)
  inline const double &xx() const { return m[0][0]; }
  inline const double &xy() const { return m[0][1]; }
  inline const double &xz() const { return m[0][2]; }
  inline const double &yx() const { return m[1][0]; }
  inline const double &yy() const { return m[1][1]; }
  inline const double &yz() const { return m[1][2]; }
  inline const double &zx() const { return m[2][0]; }
  inline const double &zy() const { return m[2][1]; }
  inline const double &zz() const { return m[2][2]; }

  // Provide operator() overloads for row/column access with bounds checking.
  inline double &operator()(int row, int col) {
    if (row < 0 || row > 2 || col < 0 || col > 2)
      throw std::out_of_range("Matrix3 indices out of range");
    return m[row][col];
  }

  inline const double &operator()(int row, int col) const {
    if (row < 0 || row > 2 || col < 0 || col > 2)
      throw std::out_of_range("Matrix3 indices out of range");
    return m[row][col];
  }

  inline Matrix3 operator+(const Matrix3 &other) const {
    Matrix3 result;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        result(i, j) = m[i][j] + other(i, j);
    return result;
  }

  inline Matrix3 &operator+=(const Matrix3 &other) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        m[i][j] += other.m[i][j];
      }
    }
    return *this;
  }

  inline Matrix3 operator*(const double scalar) const {
    Matrix3 result;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        result(i, j) = m[i][j] * scalar;
    return result;
  }
  friend Matrix3 operator*(double scalar, const Matrix3 &matrix) {
    return matrix * scalar;
  }

  inline Matrix3 operator/(const double scalar) const {
    Matrix3 result;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        result(i, j) = m[i][j] / scalar;
    return result;
  }
  friend Matrix3 operator/(double scalar, const Matrix3 &matrix) {
    return matrix / scalar;
  }

  inline Matrix3 outerProduct(const Vector3<double> &v1,
                              const Vector3<double> &v2) {
    Matrix3 result;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        result.m[i][j] = v1[i] * v2[j];
    return result;
  }

  // Override << operator to output the matrix.
  friend std::ostream &operator<<(std::ostream &os, const Matrix3 &mat) {
    os << "[ " << mat.xx() << " " << mat.xy() << " " << mat.xz() << " ],"
       << "[ " << mat.yx() << " " << mat.yy() << " " << mat.yz() << " ],"
       << "[ " << mat.zx() << " " << mat.zy() << " " << mat.zz() << " ]";
    return os;
  }
};

#endif