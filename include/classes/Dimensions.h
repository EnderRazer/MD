#ifndef DIMENSIONS_H
#define DIMENSIONS_H

#include <sstream>
#include <stdexcept>
#include <vector>

#include "nlohmann/json.hpp"

using json = nlohmann::json;

class Dimensions {
private:
  double crist_length_{0.0}; // Длина кристалической решетки

  int num_crist_x_{0}, num_crist_y_{0},
      num_crist_z_{0}; // Кол-во кристалических решеток

  int num_void_x_{0}, num_void_y_{0},
      num_void_z_{0}; // Кол-во пустых решеток (для металла)

  double size_x_{0}, size_y_{0}, size_z_{0}; // Размер системы (нм)
  double half_size_x_{0}, half_size_y_{0},
      half_size_z_{0}; // Половинчатые размеры системы
  double sqr_half_size_x_{0}, sqr_half_size_y_{0},
      sqr_half_size_z_{0}; // Квадрат половинчатых размеров системы

  double volume_{0.0}; // Объем

public:
  Dimensions() = default;
  ~Dimensions() = default;
  //
  Dimensions(const json &config) {
    double crist_length = config.value("crist_length", 1.0);
    std::vector<int> num_crist = config["crist_num"].get<std::vector<int>>();
    std::vector<int> num_void = config["void_num"].get<std::vector<int>>();

    if (crist_length < 0)
      throw std::invalid_argument("Crystal length must be non-negative");

    if (num_crist[0] < 0 || num_crist[1] < 0 || num_crist[2] < 0)
      throw std::invalid_argument("Crystal counts must be non-negative");

    if (num_void[0] < 0 || num_void[1] < 0 || num_void[2] < 0)
      throw std::invalid_argument("Void counts must be non-negative");

    crist_length_ = crist_length;
    num_crist_x_ = num_crist[0], num_crist_y_ = num_crist[1],
    num_crist_z_ = num_crist[2];
    num_void_x_ = num_void[0], num_void_y_ = num_void[1],
    num_void_z_ = num_void[2];

    size_x_ = (num_crist_x_ + num_void_x_) * crist_length_;
    size_y_ = (num_crist_y_ + num_void_y_) * crist_length_;
    size_z_ = (num_crist_z_ + num_void_z_) * crist_length_;

    recalc();
  }
  // Функция пересчета зависимых переменных.
  inline void recalc() {
    half_size_x_ = 0.5 * size_x_, half_size_y_ = 0.5 * size_y_,
    half_size_z_ = 0.5 * size_z_;

    sqr_half_size_x_ = half_size_x_ * half_size_x_,
    sqr_half_size_y_ = half_size_y_ * half_size_y_,
    sqr_half_size_z_ = half_size_z_ * half_size_z_;

    volume_ = size_x_ * size_y_ * size_z_;
  }

  std::string getData() const {
    std::ostringstream oss;
    oss.precision(16);
    oss << "System size params:\n\t Crist lenght: " << crist_length_
        << "\n\t Crist num (x,y,z): (" << num_crist_x_ << ", " << num_crist_y_
        << ", " << num_crist_z_ << ")"
        << "\n\t Void num (x,y,z): (" << num_void_x_ << ", " << num_void_y_
        << ", " << num_void_z_ << ")"
        << "\n\t Size (lx,ly,lz): (" << size_x_ << ", " << size_y_ << ", "
        << size_z_ << ")"
        << "\n\t Volume: " << volume_ << "\n";

    return oss.str();
  }

  inline double &cristLength() { return crist_length_; }
  inline const double &cristLength() const { return crist_length_; }

  inline int &numCristX() { return num_crist_x_; }
  inline const int &numCristX() const { return num_crist_x_; }
  inline int &numCristY() { return num_crist_y_; }
  inline const int &numCristY() const { return num_crist_y_; }
  inline int &numCristZ() { return num_crist_z_; }
  inline const int &numCristZ() const { return num_crist_z_; }

  inline int &numVoidX() { return num_void_x_; }
  inline const int &numVoidX() const { return num_void_x_; }
  inline int &numVoidY() { return num_void_y_; }
  inline const int &numVoidY() const { return num_void_y_; }
  inline int &numVoidZ() { return num_void_z_; }
  inline const int &numVoidZ() const { return num_void_z_; }

  inline double &sizeX() { return size_x_; }
  inline const double &sizeX() const { return size_x_; }
  inline double &sizeY() { return size_y_; }
  inline const double &sizeY() const { return size_y_; }
  inline double &sizeZ() { return size_z_; }
  inline const double &sizeZ() const { return size_z_; }

  inline double &halfSizeX() { return half_size_x_; }
  inline const double &halfSizeX() const { return half_size_x_; }
  inline double &halfSizeY() { return half_size_y_; }
  inline const double &halfSizeY() const { return half_size_y_; }
  inline double &halfSizeZ() { return half_size_z_; }
  inline const double &halfSizeZ() const { return half_size_z_; }

  inline double &sqrHalfSizeX() { return sqr_half_size_x_; }
  inline const double &sqrHalfSizeX() const { return sqr_half_size_x_; }
  inline double &sqrHalfSizeY() { return sqr_half_size_y_; }
  inline const double &sqrHalfSizeY() const { return sqr_half_size_y_; }
  inline double &sqrHalfSizeZ() { return sqr_half_size_z_; }
  inline const double &sqrHalfSizeZ() const { return sqr_half_size_z_; }

  inline double &volume() { return volume_; }
  inline const double &volume() const { return volume_; }

  Dimensions(const Dimensions &) = delete;
  Dimensions &operator=(const Dimensions &) = delete;
};

#endif // DIMENSIONS_H
