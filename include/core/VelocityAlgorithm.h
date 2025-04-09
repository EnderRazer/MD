#ifndef VELOCITY_ALGORITHM_H
#define VELOCITY_ALGORITHM_H

#include <classes/Particle.h>
#include <core/Settings.h>

class VelocityAlgorithm {
private:
  Settings &settings_;
  const double mt_;

public:
  inline void compute(Particle &p) { p.addVelocity(p.force() * mt_); }

  VelocityAlgorithm() = delete;
  VelocityAlgorithm(Settings &settings)
      : settings_(settings), mt_(settings_.dt() / (2 * settings_.mass())) {}
  ~VelocityAlgorithm() = default;

  // Запрещаем копирование
  VelocityAlgorithm(const VelocityAlgorithm &) = delete;
  VelocityAlgorithm &operator=(const VelocityAlgorithm &) = delete;
};
#endif
