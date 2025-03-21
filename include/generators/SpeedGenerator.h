#ifndef SPEED_GENERATOR_H
#define SPEED_GENERATOR_H

#include "core/Settings.h"
#include "core/System.h"
#include "random_normalized.h"

class SpeedGenerator {
private:
  double sigma_maxwell_{0.0}; // Константа для распределения скоростей
public:
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
  SpeedGenerator() = default;
  SpeedGenerator(Settings &settings) {
    sigma_maxwell_ = sqrt(settings.constants().KBOLTZMAN *
                          (settings.startTemp()) / settings.mass());
  };
  ;
  ~SpeedGenerator() = default;

  // Запрещаем копирование
  SpeedGenerator(const SpeedGenerator &) = delete;
  SpeedGenerator &operator=(const SpeedGenerator &) = delete;
};

#endif
