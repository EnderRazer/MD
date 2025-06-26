#ifndef VELOCITY_ALGORITHM_H
#define VELOCITY_ALGORITHM_H

#include <classes/Particle.h>
#include <core/Settings.h>

/**
 * @brief Класс для работы с алгоритмом скоростей.
 * @details Класс предоставляет методы для работы с алгоритмом скоростей.
 */
class VelocityAlgorithm {
private:
  /**
   * @brief Настройки.
   */
  Settings &settings_;

  /**
   * @brief Константа для расчета скоростей.
   * @details Константа для расчета скоростей. dt / (2 * mass).
   */
  const double mt_;

public:
  /**
   * @brief Конструктор.
   * @param settings - настройки.
   */
  VelocityAlgorithm(Settings &settings)
      : settings_(settings), mt_(settings_.dt() / (2 * settings_.mass())) {}

  /**
   * @brief Расчет скоростей.
   * @param p - частица.
   * @details Расчет скоростей частицы.
   */
  inline void compute(Particle &p) { p.addVelocity(p.force() * mt_); }

  VelocityAlgorithm() = delete;
  ~VelocityAlgorithm() = default;

  // Запрещаем копирование
  VelocityAlgorithm(const VelocityAlgorithm &) = delete;
  VelocityAlgorithm &operator=(const VelocityAlgorithm &) = delete;
};
#endif
