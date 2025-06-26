#ifndef SPEED_GENERATOR_H
#define SPEED_GENERATOR_H

#include "core/Settings.h"
#include "core/System.h"
#include "helpers/random_normalized.h"

/**
 * @brief Класс для генерации скоростей частиц.
 * @details Класс предоставляет методы для генерации скоростей частиц.
 */
class SpeedGenerator {
private:
  /**
   * @brief Константа для распределения скоростей.
   * @details Константа для распределения скоростей.
   */
  double sigma_maxwell_{0.0};

public:
  SpeedGenerator() = default;
  ~SpeedGenerator() = default;

  // Запрещаем копирование
  SpeedGenerator(const SpeedGenerator &) = delete;
  SpeedGenerator &operator=(const SpeedGenerator &) = delete;

  /**
   * @brief Конструктор.
   * @param settings - настройки.
   */
  SpeedGenerator(Settings &settings) {
    sigma_maxwell_ = sqrt(settings.constants().KBOLTZMAN *
                          (settings.startTemp()) / settings.mass());
  };

  /**
   * @brief Генерация начальных скоростей.
   * @param sys - система частиц.
   * @param settings - настройки.
   */
  void startingSpeeds(System &sys, Settings &settings) {
    int pn = sys.particleNumber();
    std::vector<Particle> &particles = sys.particles();
    for (int i = 0; i < pn / 2; i++) {
      Vector3<double> speed =
          Vector3<double>(sigma_maxwell_ * gasdev(settings.seed()),
                          sigma_maxwell_ * gasdev(settings.seed()),
                          sigma_maxwell_ * gasdev(settings.seed()));
      particles[i].setVelocity(speed);
      particles[i + pn / 2].setVelocity(-1 * speed);
    }
  }
};

#endif
