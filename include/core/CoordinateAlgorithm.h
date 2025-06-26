#ifndef COORDINATE_ALGORITHM_H
#define COORDINATE_ALGORITHM_H

#include "classes/Dimensions.h"
#include "classes/Particle.h"
#include "core/Settings.h"

/**
 * @brief Класс для вычисления координат частицы.
 * @details Класс предоставляет методы для вычисления координат частицы и
 * применения периодических граничных условий.
 */
class CoordinateAlgorithm {
private:
  /**
   * @brief Настройки.
   */
  Settings &settings_;

  
public:
/**
   * @brief Вычисление координат частицы.
   */
  inline void compute(Particle &p) {
    p.addCoord(p.velocity() * settings_.dt());
  }

  /**
   * @brief Применение периодических граничных условий.
   */
  inline void applyPBC(Particle &p, Dimensions dim) {
    Vector3<double> coords = p.coord();
    Vector3<double> dimSizes = dim.sizes();
    if (coords.x() < 0)
      coords.x() += dimSizes.x();
    if (coords.x() >= dimSizes.x())
      coords.x() -= dimSizes.x();

    if (coords.y() < 0)
      coords.y() += dimSizes.y();
    if (coords.y() >= dimSizes.y())
      coords.y() -= dimSizes.y();

    if (coords.z() < 0)
      coords.z() += dimSizes.z();
    if (coords.z() >= dimSizes.z())
      coords.z() -= dimSizes.z();

    p.setCoord(coords);
  }

  /**
   * @brief Конструктор.
   */
  CoordinateAlgorithm() = delete;

  /**
   * @brief Конструктор.
   */
  CoordinateAlgorithm(Settings &settings) : settings_(settings) {}

  /**
   * @brief Деструктор.
   */
  ~CoordinateAlgorithm() = default;

  CoordinateAlgorithm(const CoordinateAlgorithm &) = delete;
  CoordinateAlgorithm &operator=(const CoordinateAlgorithm &) = delete;
};
#endif // COORDINATE_ALGORITHM_H
