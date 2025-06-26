#ifndef VECTOR3_H
#define VECTOR3_H

#include <cmath>
#include <cstddef>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <vector>

/**
 * @brief Класс для работы с векторами в трехмерном пространстве.
 * @details Класс предоставляет интерфейс для работы с векторами в трехмерном
 * пространстве. Он позволяет создавать вектора, получать их компоненты,
 * выполнять операции сложения, вычитания, умножения и деления векторов,
 * а также выполнять другие операции.
 */
template <typename T> class Vector3 {
private:
  /**
   * @brief Компоненты вектора.
   */
  T components[3];

public:
  /**
   * @brief Конструктор по умолчанию.
   */
  Vector3() : components{T(0), T(0), T(0)} {}

  /**
   * @brief Конструктор с заданными значениями.
   * @param vec - вектор.
   * @details Конструктор создает вектор с заданными значениями.
   */
  Vector3(const std::vector<T> &vec) {
    if (vec.size() != 3) {
      throw std::invalid_argument("Vector3 can only be initialized from a "
                                  "std::vector with exactly 3 elements.");
    }
    components[0] = vec[0];
    components[1] = vec[1];
    components[2] = vec[2];
  }

  /**
   * @brief Конструктор с заданными значениями.
   * @param x - значение x.
   * @param y - значение y.
   * @param z - значение z.
   * @details Конструктор создает вектор с заданными значениями.
   */
  Vector3(T x, T y, T z) : components{x, y, z} {}

  /**
   * @brief Получение значения x.
   * @return Значение x.
   */
  T &x() { return components[0]; }

  /**
   * @brief Получение значения y.
   * @return Значение y.
   */
  T &y() { return components[1]; }

  /**
   * @brief Получение значения z.
   * @return Значение z.
   */
  T &z() { return components[2]; }

  /**
   * @brief Получение значения x.
   * @return Значение x.
   */
  const T &x() const { return components[0]; }

  /**
   * @brief Получение значения y.
   * @return Значение y.
   */
  const T &y() const { return components[1]; }

  /**
   * @brief Получение значения z.
   * @return Значение z.
   */
  const T &z() const { return components[2]; }

  /**
   * @brief Конвертирующий конструктор.
   * @param other - другой вектор.
   * @details Конвертирующий конструктор создает вектор с заданными значениями.
   */
  template <typename U> Vector3(const Vector3<U> &other) {
    components[0] = static_cast<T>(other.x());
    components[1] = static_cast<T>(other.y());
    components[2] = static_cast<T>(other.z());
  }

  /**
   * @brief Перегрузка оператора [].
   * @param index - индекс.
   * @return Значение вектора.
   */
  T &operator[](size_t index) {
    if (index >= 3)
      throw std::out_of_range("Index out of range in Vector3");
    return components[index];
  }

  /**
   * @brief Перегрузка оператора [].
   * @param index - индекс.
   * @return Значение вектора.
   */
  const T &operator[](size_t index) const {
    if (index >= 3)
      throw std::out_of_range("Index out of range in Vector3");
    return components[index];
  }

  /**
   * @brief Получение квадрата длины вектора.
   * @return Квадрат длины вектора.
   */
  T lengthSquared() const {
    return components[0] * components[0] + components[1] * components[1] +
           components[2] * components[2];
  }

  /**
   * @brief Получение длины вектора.
   * @return Длина вектора.
   */
  T length() const { return std::sqrt(lengthSquared()); }

  /**
   * @brief Перегрузка оператора *.
   * @param scalar - скаляр.
   * @return Результат умножения вектора на скаляр.
   */
  Vector3<T> operator*(T scalar) const {
    return Vector3<T>(components[0] * scalar, components[1] * scalar,
                      components[2] * scalar);
  }

  /**
   * @brief Перегрузка оператора *.
   * @param scalar - скаляр.
   * @param vec - вектор.
   * @return Результат умножения вектора на скаляр.
   */
  friend Vector3<T> operator*(T scalar, const Vector3<T> &vec) {
    return vec * scalar;
  }

  /**
   * @brief Перегрузка оператора +.
   * @param scalar - скаляр.
   * @return Результат сложения вектора с скаляром.
   */
  Vector3<T> operator+(T scalar) const {
    return Vector3<T>(components[0] + scalar, components[1] + scalar,
                      components[2] + scalar);
  }

  /**
   * @brief Перегрузка оператора +.
   * @param scalar - скаляр.
   * @param vec - вектор.
   * @return Результат сложения вектора с скаляром.
   */
  friend Vector3<T> operator+(T scalar, const Vector3<T> &vec) {
    return vec + scalar;
  }

  /**
   * @brief Перегрузка оператора -.
   * @param scalar - скаляр.
   * @return Результат вычитания скаляра из вектора.
   */
  Vector3<T> operator-(T scalar) const {
    return Vector3<T>(components[0] - scalar, components[1] - scalar,
                      components[2] - scalar);
  }

  /**
   * @brief Перегрузка оператора -.
   * @param scalar - скаляр.
   * @param vec - вектор.
   * @return Результат вычитания скаляра из вектора.
   */
  friend Vector3<T> operator-(T scalar, const Vector3<T> &vec) {
    return vec - scalar;
  }

  /**
   * @brief Перегрузка оператора /.
   * @param scalar - скаляр.
   * @return Результат деления вектора на скаляр.
   */
  Vector3<T> operator/(T scalar) const {
    return Vector3<T>(components[0] / scalar, components[1] / scalar,
                      components[2] / scalar);
  }

  /**
   * @brief Перегрузка оператора /.
   * @param scalar - скаляр.
   * @param vec - вектор.
   * @return Результат деления вектора на скаляр.
   */
  friend Vector3<T> operator/(T scalar, const Vector3<T> &vec) {
    return vec / scalar;
  }

  /**
   * @brief Перегрузка оператора *=.
   * @param scalar - скаляр.
   * @return Результат умножения вектора на скаляр.
   */
  Vector3<T> &operator*=(T scalar) {
    components[0] *= scalar;
    components[1] *= scalar;
    components[2] *= scalar;
    return *this;
  }

  /**
   * @brief Перегрузка оператора /=.
   * @param scalar - скаляр.
   * @return Результат деления вектора на скаляр.
   */
  Vector3<T> &operator/=(T scalar) {
    components[0] /= scalar;
    components[1] /= scalar;
    components[2] /= scalar;
    return *this;
  }

  /**
   * @brief Перегрузка оператора +=.
   * @param scalar - скаляр.
   * @return Результат сложения вектора с скаляром.
   */
  Vector3<T> &operator+=(T scalar) {
    components[0] += scalar;
    components[1] += scalar;
    components[2] += scalar;
    return *this;
  }

  /**
   * @brief Перегрузка оператора -=.
   * @param scalar - скаляр.
   * @return Результат вычитания скаляра из вектора.
   */
  Vector3<T> &operator-=(T scalar) {
    components[0] -= scalar;
    components[1] -= scalar;
    components[2] -= scalar;
    return *this;
  }

  /**
   * @brief Перегрузка оператора +.
   * @param other - другой вектор.
   * @return Результат сложения вектора с другим вектором.
   */
  Vector3<T> operator+(const Vector3<T> &other) const {
    return Vector3<T>(components[0] + other.components[0],
                      components[1] + other.components[1],
                      components[2] + other.components[2]);
  }

  /**
   * @brief Перегрузка оператора +=.
   * @param other - другой вектор.
   * @return Результат сложения вектора с другим вектором.
   */
  Vector3<T> operator+=(const Vector3<T> &other) {
    components[0] += other.components[0];
    components[1] += other.components[1];
    components[2] += other.components[2];
    return *this;
  }

  /**
   * @brief Перегрузка оператора -.
   * @param other - другой вектор.
   * @return Результат вычитания вектора из другого вектора.
   */
  Vector3<T> operator-(const Vector3<T> &other) const {
    return Vector3<T>(components[0] - other.components[0],
                      components[1] - other.components[1],
                      components[2] - other.components[2]);
  }

  /**
   * @brief Перегрузка оператора -=.
   * @param other - другой вектор.
   * @return Результат вычитания вектора из другого вектора.
   */
  Vector3<T> operator-=(const Vector3<T> &other) {
    components[0] -= other.components[0];
    components[1] -= other.components[1];
    components[2] -= other.components[2];
    return *this;
  }

  /**
   * @brief Перегрузка оператора <<.
   * @return Результат вывода вектора.
   */
  friend std::ostream &operator<<(std::ostream &os, const Vector3<T> &vec) {
    os << "(" << vec.components[0] << ", " << vec.components[1] << ", "
       << vec.components[2] << ")";
    return os;
  }

  /**
   * @brief Преобразование вектора в строку.
   * @return Строка, представляющая вектор.
   */
  inline const std::string to_str() const {
    std::ostringstream os;
    os.precision(16);
    os << "(" << components[0] << ", " << components[1] << ", " << components[2]
       << ")";
    return os.str();
  }
};

#endif
