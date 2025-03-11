#ifndef DIMENSIONS_H
#define DIMENSIONS_H

#include <stdexcept>
using json = nlohmann::json;

class Dimensions {
  private:
    double crist_length_;                           // Длина кристалической решетки
    Vector3<int> num_crist_{0, 0, 0};               // Кол-во кристалических решеток
    Vector3<double> sizes_{0.0, 0.0, 0.0};          // Размер системы (нм)
    Vector3<double> half_sizes_{0.0, 0.0, 0.0};     // Половинчатые размеры системы
    Vector3<double> sqr_half_sizes_{0.0, 0.0, 0.0}; // Квадрат половинчатых размеров системы
    double volume_{0.0};                            // Объем

  public:
    Dimensions() = default;
    //
    explicit Dimensions(const json &config)
        : crist_length_(config.value("crist_length", 1.0)),
          num_crist_{config.value("crist_num_x", 1), config.value("crist_num_y", 1), config.value("crist_num_z", 1)} {
        if (crist_length_ < 0)
            throw std::invalid_argument("Crystal length must be non-negative");

        if (num_crist_.x() < 0 || num_crist_.y() < 0 || num_crist_.z() < 0)
            throw std::invalid_argument("Crystal counts must be non-negative");

        sizes_ = {num_crist_.x() * crist_length_, num_crist_.y() * crist_length_, num_crist_.z() * crist_length_};
        recalc();
    }

    // -------------------------
    // Геттеры
    // -------------------------
    inline double cristLength() const { return crist_length_; }

    inline Vector3<int> numCrists() const { return num_crist_; }
    inline Vector3<double> sizes() const { return sizes_; }
    inline Vector3<double> halfSizes() const { return half_sizes_; }
    inline Vector3<double> sqrHalfSizes() const { return sqr_half_sizes_; }

    inline int numCristX() const { return num_crist_.x(); }
    inline int numCristY() const { return num_crist_.y(); }
    inline int numCristZ() const { return num_crist_.z(); }

    inline double lx() const { return sizes_.x(); }
    inline double ly() const { return sizes_.y(); }
    inline double lz() const { return sizes_.z(); }

    inline double halfLX() const { return half_sizes_.x(); }
    inline double halfLY() const { return half_sizes_.y(); }
    inline double halfLZ() const { return half_sizes_.z(); }

    inline double sqrHalfLX() const { return sqr_half_sizes_.x(); }
    inline double sqrHalfLY() const { return sqr_half_sizes_.y(); }
    inline double sqrHalfLZ() const { return sqr_half_sizes_.z(); }

    inline double volume() const { return volume_; }
    // -------------------------
    // Сеттеры
    // -------------------------
    inline void setCristLength(double crist_length) {
        if (crist_length < 0) {
            throw std::invalid_argument("Длина кристаллической решетки не может быть отрицательной");
        }
        crist_length_ = crist_length;
        sizes_.x() = (num_crist_.x() * crist_length_);
        sizes_.y() = (num_crist_.y() * crist_length_);
        sizes_.z() = (num_crist_.z() * crist_length_);
        recalc();
    }

    // Изменение кол-ва кристалических решеток
    inline void setNumCrists(int num_x, int num_y, int num_z) {
        if (num_x < 0 || num_y < 0 || num_z < 0)
            throw std::invalid_argument("Кол-во кристаллических решеток не может быть отрицательным");

        num_crist_ = {num_x, num_y, num_z};
        sizes_ = {num_x * crist_length_, num_y * crist_length_, num_z * crist_length_};
        recalc();
    }
    inline void setNumCristX(int numcrist_x) {
        if (numcrist_x < 0)
            throw std::invalid_argument("Кол-во кристалических решеток не может быть отрицательным");
        num_crist_.x() = numcrist_x;
        sizes_.x() = num_crist_.x() * crist_length_;
        recalc();
    }
    inline void setNumCristY(int numcrist_y) {
        if (numcrist_y < 0)
            throw std::invalid_argument("Кол-во кристалических решеток не может быть отрицательным");
        num_crist_.y() = numcrist_y;
        sizes_.y() = num_crist_.y() * crist_length_;
        recalc();
    }
    inline void setNumCristZ(int numcrist_z) {
        if (numcrist_z < 0)
            throw std::invalid_argument("Кол-во кристалических решеток не может быть отрицательным");
        num_crist_.z() = numcrist_z;
        sizes_.z() = num_crist_.z() * crist_length_;
        recalc();
    }

    inline void setSizes(Vector3<double> sizes) {
        if (sizes.x() < 0 || sizes.y() < 0 || sizes.z() < 0)
            throw std::invalid_argument("Размеры должны быть положительными");

        sizes_ = sizes;
        recalc();
    }
    inline void setSizes(double l_x, double l_y, double l_z) {
        if (l_x < 0 || l_y < 0 || l_z < 0)
            throw std::invalid_argument("Размеры должны быть положительными");

        sizes_ = {l_x, l_y, l_z};
        recalc();
    }
    // Меняем размеры напрямую
    inline void setLX(double l_x) {
        if (l_x < 0)
            throw std::invalid_argument("Длина должна быть положительной");
        sizes_.x() = l_x;
        recalc();
    }
    inline void setLY(double l_y) {
        if (l_y < 0)
            throw std::invalid_argument("Длина должна быть положительной");
        sizes_.y() = l_y;
        recalc();
    }
    inline void setLZ(double l_z) {
        if (l_z < 0)
            throw std::invalid_argument("Длина должна быть положительной");
        sizes_.z() = l_z;
        recalc();
    }

    // Функция пересчета зависимых переменных.
    inline void recalc() {
        half_sizes_ = {sizes_.x() / 2.0, sizes_.y() / 2, sizes_.z() / 2};
        sqr_half_sizes_ = {half_sizes_.x() * half_sizes_.x(), half_sizes_.y() * half_sizes_.y(),
                           half_sizes_.z() * half_sizes_.z()};
        volume_ = sizes_.x() * sizes_.y() * sizes_.z();
    }

    std::string getData() const {
        std::ostringstream oss;
        oss.precision(16);
        oss << "System size params:\n\t Crist lenght: " << crist_length_ << "\n\t Crist num (x,y,z): " << num_crist_
            << "\n\t Size (lx,ly,lz): " << sizes_ << "\n\t Volume: " << volume_ << "\n";

        return oss.str();
    }
};

#endif // DIMENSIONS_H