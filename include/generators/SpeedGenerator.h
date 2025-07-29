#ifndef SPEED_GENERATOR_H
#define SPEED_GENERATOR_H

#include "classes/Particles.h"
#include "core/Settings.h"
#include "core/System.h"
#include "random_normalized.h"
#include <ostream>

class SpeedGenerator {
private:
  double sigma_maxwell_{0.0}; // Константа для распределения скоростей
public:
  void startingSpeeds(System &sys, Settings &settings) {
    Particles &particles = sys.particles();
    int pn = particles.size();
    double speed_x,speed_y,speed_z;
    for (int i = 0; i < pn / 2; i++) {
      speed_x = sigma_maxwell_ * gasdev(settings.seed());
      speed_y = sigma_maxwell_ * gasdev(settings.seed());
      speed_z = sigma_maxwell_ * gasdev(settings.seed());

      particles.velocityX(i)= speed_x;
      particles.velocityY(i)= speed_y;
      particles.velocityZ(i)= speed_z;

      particles.velocityX(i + pn / 2) = -speed_x;
      particles.velocityY(i + pn / 2) = -speed_y;
      particles.velocityZ(i + pn / 2) = -speed_z;
    }
    //std::cout << sys.getParticlesVelocityInfo() << std::endl;
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
