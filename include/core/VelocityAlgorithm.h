#ifndef VELOCITY_ALGORITHM_H
#define VELOCITY_ALGORITHM_H

#include "classes/Timer.h"
#include <classes/Particles.h>
#include <core/Settings.h>

class VelocityAlgorithm {
private:
  Settings &settings_;
  const double mt_;

public:
  inline void compute(Particles &p) { 
    Timer timer{1};
    timer.start();
    p.computeVelocitiesSIMD(mt_);
    timer.stop();
    std::cout << "Velocity time elapsed = " << timer.elapsed() << "ms"<<std::endl;
  }

  VelocityAlgorithm() = delete;
  VelocityAlgorithm(Settings &settings)
      : settings_(settings), mt_(settings_.dt() / (2 * settings_.mass())) {}
  ~VelocityAlgorithm() = default;

  // Запрещаем копирование
  VelocityAlgorithm(const VelocityAlgorithm &) = delete;
  VelocityAlgorithm &operator=(const VelocityAlgorithm &) = delete;
};
#endif
