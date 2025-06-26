#ifndef DIMENSIONS_H
#define DIMENSIONS_H

#include <sstream>
#include <stdexcept>

#include "classes/Vector3.h"
#include "nlohmann/json.hpp"

using json = nlohmann::json;

/**
 * @brief Класс для хранения размеров системы.
 *
 * Класс для хранения размеров системы.
 */
class Dimensions {
private:
  /**
   * @brief Длина кристалической решетки.
   */
  double crist_length_;

  /**
   * @brief Кол-во кристалических решеток.
   */
  Vector3<int> num_crist_{0, 0, 0};

  /**
   * @brief Кол-во пустых решеток (для металла).
   */
  Vector3<int> num_void_{0, 0, 0};

  /**
   * @brief Размер системы (нм).
   */
  Vector3<double> sizes_{0.0, 0.0, 0.0};

  /**
   * @brief Половинчатые размеры системы.
   */
  Vector3<double> half_sizes_{0.0, 0.0, 0.0};

  /**
   * @brief Квадрат половинчатых размеров системы.
   */
  Vector3<double> sqr_half_sizes_{0.0, 0.0, 0.0};

  /**
   * @brief Объем системы.
   */
  double volume_{0.0};

public:
  Dimensions() = default;
  //
  explicit Dimensions(const json &config)
      : crist_length_(config.value("crist_length", 1.0)),
        num_crist_{config["crist_num"].get<std::vector<int>>()},
        num_void_{config["void_num"].get<std::vector<int>>()} {
    if (crist_length_ < 0)
      throw std::invalid_argument("Crystal length must be non-negative");

    if (num_crist_.x() < 0 || num_crist_.y() < 0 || num_crist_.z() < 0)
      throw std::invalid_argument("Crystal counts must be non-negative");

    if (num_void_.x() < 0 || num_void_.y() < 0 || num_void_.z() < 0)
      throw std::invalid_argument("Void counts must be non-negative");

    sizes_ = {(num_crist_.x() + num_void_.x()) * crist_length_,
              (num_crist_.y() + num_void_.y()) * crist_length_,
              (num_crist_.z() + num_void_.z()) * crist_length_};
    recalc();
  }

  // -------------------------
  // Геттеры
  // -------------------------
  /**
   * @brief Получение длины кристалической решетки.
   *
   * Получение длины кристалической решетки.
   * @return длина кристалической решетки.
   */
  inline double cristLength() const { return crist_length_; }

  /**
   * @brief Получение количества кристалических решеток.
   *
   * Получение количества кристалических решеток.
   * @return вектор с количеством кристалических решеток.
   */
  inline Vector3<int> numCrists() const { return num_crist_; }

  /**
   * @brief Получение количества пустых решеток.
   *
   * Получение количества пустых решеток.
   * @return вектор с количеством пустых решеток.
   */
  inline Vector3<int> numVoid() const { return num_void_; }

  /**
   * @brief Получение размеров системы.
   *
   * Получение размеров системы.
   * @return вектор с размерами системы.
   */
  inline Vector3<double> sizes() const { return sizes_; }

  /**
   * @brief Получение квадрата половинчатых размеров системы.
   *
   * Получение квадрата половинчатых размеров системы.
   * @return вектор с квадратами половинчатых размеров системы.
   */
  inline Vector3<double> sqrHalfSizes() const { return sqr_half_sizes_; }

  /**
   * @brief Получение количества кристалических решеток по оси X.
   *
   * Получение количества кристалических решеток по оси X.
   * @return количество кристалических решеток по оси X.
   */
  inline int numCristX() const { return num_crist_.x(); }

  /**
   * @brief Получение количества кристалических решеток по оси Y.
   *
   * Получение количества кристалических решеток по оси Y.
   * @return количество кристалических решеток по оси Y.
   */
  inline int numCristY() const { return num_crist_.y(); }

  /**
   * @brief Получение количества кристалических решеток по оси Z.
   *
   * Получение количества кристалических решеток по оси Z.
   * @return количество кристалических решеток по оси Z.
   */
  inline int numCristZ() const { return num_crist_.z(); }

  /**
   * @brief Получение длины системы по оси X.
   *
   * Получение длины системы по оси X.
   * @return длина системы по оси X.
   */
  inline double lx() const { return sizes_.x(); }

  /**
   * @brief Получение длины системы по оси Y.
   *
   * Получение длины системы по оси Y.
   * @return длина системы по оси Y.
   */
  inline double ly() const { return sizes_.y(); }

  /**
   * @brief Получение длины системы по оси Z.
   *
   * Получение длины системы по оси Z.
   * @return длина системы по оси Z.
   */
  inline double lz() const { return sizes_.z(); }

  /**
   * @brief Получение половины длины системы по оси X.
   *
   * Получение половины длины системы по оси X.
   * @return половина длины системы по оси X.
   */
  inline double halfLX() const { return half_sizes_.x(); }

  /**
   * @brief Получение половины длины системы по оси Y.
   *
   * Получение половины длины системы по оси Y.
   * @return половина длины системы по оси Y.
   */
  inline double halfLY() const { return half_sizes_.y(); }

  /**
   * @brief Получение половины длины системы по оси Z.
   *
   * Получение половины длины системы по оси Z.
   * @return половина длины системы по оси Z.
   */
  inline double halfLZ() const { return half_sizes_.z(); }

  /**
   * @brief Получение объема системы.
   *
   * Получение объема системы.
   * @return объем системы.
   */
  inline double volume() const { return volume_; }

  // -------------------------
  // Сеттеры
  // -------------------------
  /**
   * @brief Установка длины кристалической решетки.
   *
   * Установка длины кристалической решетки. И пересчет всех зависимых переменных.
   * @param crist_length - длина кристалической решетки.
   */
  inline void setCristLength(double crist_length) {
    if (crist_length < 0) {
      throw std::invalid_argument(
          "Длина кристаллической решетки не может быть отрицательной");
    }
    crist_length_ = crist_length;
    sizes_.x() = (num_crist_.x() * crist_length_);
    sizes_.y() = (num_crist_.y() * crist_length_);
    sizes_.z() = (num_crist_.z() * crist_length_);
    recalc();
  }

  /**
   * @brief Установка количества кристалических решеток.
   *
   * Установка количества кристалических решеток. И пересчет всех зависимых переменных.
   * @param num_x - количество кристалических решеток по оси X.
   * @param num_y - количество кристалических решеток по оси Y.
   * @param num_z - количество кристалических решеток по оси Z.
   */
  inline void setNumCrists(int num_x, int num_y, int num_z) {
    if (num_x < 0 || num_y < 0 || num_z < 0)
      throw std::invalid_argument(
          "Кол-во кристаллических решеток не может быть отрицательным");

    num_crist_ = {num_x, num_y, num_z};
    sizes_ = {num_x * crist_length_, num_y * crist_length_,
              num_z * crist_length_};
    recalc();
  }

  /**
   * @brief Установка количества кристалических решеток по оси X.
   *
   * Установка количества кристалических решеток по оси X. И пересчет всех зависимых переменных.
   * @param numcrist_x - количество кристалических решеток по оси X.
   */
  inline void setNumCristX(int numcrist_x) {
    if (numcrist_x < 0)
      throw std::invalid_argument(
          "Кол-во кристалических решеток не может быть отрицательным");
    num_crist_.x() = numcrist_x;
    sizes_.x() = num_crist_.x() * crist_length_;
    recalc();
  }

  /**
   * @brief Установка количества кристалических решеток по оси Y.
   *
   * Установка количества кристалических решеток по оси Y. И пересчет всех зависимых переменных.
   * @param numcrist_y - количество кристалических решеток по оси Y.
   */
  inline void setNumCristY(int numcrist_y) {
    if (numcrist_y < 0)
      throw std::invalid_argument(
          "Кол-во кристалических решеток не может быть отрицательным");
    num_crist_.y() = numcrist_y;
    sizes_.y() = num_crist_.y() * crist_length_;
    recalc();
  }

  /**
   * @brief Установка количества кристалических решеток по оси Z.
   *
   * Установка количества кристалических решеток по оси Z. И пересчет всех зависимых переменных.
   * @param numcrist_z - количество кристалических решеток по оси Z.
   */
  inline void setNumCristZ(int numcrist_z) {
    if (numcrist_z < 0)
      throw std::invalid_argument(
          "Кол-во кристалических решеток не может быть отрицательным");
    num_crist_.z() = numcrist_z;
    sizes_.z() = num_crist_.z() * crist_length_;
    recalc();
  }

  /**
   * @brief Установка размеров системы.
   *
   * Установка размеров системы. И пересчет всех зависимых переменных.
   * @param sizes - вектор с размерами системы.
   */
  inline void setSizes(Vector3<double> sizes) {
    if (sizes.x() < 0 || sizes.y() < 0 || sizes.z() < 0)
      throw std::invalid_argument("Размеры должны быть положительными");

    sizes_ = sizes;
    recalc();
  }

  /**
   * @brief Установка размеров системы.
   *
   * Установка размеров системы. И пересчет всех зависимых переменных.
   * @param l_x - длина системы по оси X.
   * @param l_y - длина системы по оси Y.
   * @param l_z - длина системы по оси Z.
   */
  inline void setSizes(double l_x, double l_y, double l_z) {
    if (l_x < 0 || l_y < 0 || l_z < 0)
      throw std::invalid_argument("Размеры должны быть положительными");

    sizes_ = {l_x, l_y, l_z};
    recalc();
  }

  /**
   * @brief Установка длины системы по оси X.
   *
   * Установка длины системы по оси X. И пересчет всех зависимых переменных.
   * @param l_x - длина системы по оси X.
   */
  inline void setLX(double l_x) {
    if (l_x < 0)
      throw std::invalid_argument("Длина должна быть положительной");
    sizes_.x() = l_x;
    recalc();
  }

  /**
   * @brief Установка длины системы по оси Y.
   *
   * Установка длины системы по оси Y. И пересчет всех зависимых переменных.
   * @param l_y - длина системы по оси Y.
   */
  inline void setLY(double l_y) {
    if (l_y < 0)
      throw std::invalid_argument("Длина должна быть положительной");
    sizes_.y() = l_y;
    recalc();
  }

  /**
   * @brief Установка длины системы по оси Z.
   *
   * Установка длины системы по оси Z. И пересчет всех зависимых переменных.
   * @param l_z - длина системы по оси Z.
   */
  inline void setLZ(double l_z) {
    if (l_z < 0)
      throw std::invalid_argument("Длина должна быть положительной");
    sizes_.z() = l_z;
    recalc();
  }

  /**
   * @brief Пересчет всех зависимых переменных.
   *
   * Пересчет всех зависимых переменных.
   */
  inline void recalc() {
    half_sizes_ = {sizes_.x() / 2.0, sizes_.y() / 2, sizes_.z() / 2};
    sqr_half_sizes_ = {half_sizes_.x() * half_sizes_.x(),
                       half_sizes_.y() * half_sizes_.y(),
                       half_sizes_.z() * half_sizes_.z()};
    volume_ = sizes_.x() * sizes_.y() * sizes_.z();
  }

  /**
   * @brief Получение базовой информации о размерах системы.
   *
   * Получение базовой информации о размерах системы.
   * @return строка с информацией о размерах системы.
   */
  std::string getData() const {
    std::ostringstream oss;
    oss.precision(16);
    oss << "System size params:\n\t Crist lenght: " << crist_length_
        << "\n\t Crist num (x,y,z): " << num_crist_
        << "\n\t Void num (x,y,z): " << num_void_
        << "\n\t Size (lx,ly,lz): " << sizes_ << "\n\t Volume: " << volume_
        << "\n";

    return oss.str();
  }
};

#endif // DIMENSIONS_H
