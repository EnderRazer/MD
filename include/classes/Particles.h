#ifndef PARTICLE_H
#define PARTICLE_H

#include <cstddef>
#include <iostream>
#include <vector>

#include "AlignedAllocator.h"

using AlignedVector = std::vector<double, AlignedAllocator<double, 64>>;

struct ForceCalcValues {
  int interaction_count{0};
  double rVec_x{0}, rVec_y{0}, rVec_z{0};

  double force_x{0}, force_y{0}, force_z{0};

  double e_pot{0.0};

  double virial_xx{0}, virial_xy{0}, virial_xz{0};
  double virial_yx{0}, virial_yy{0}, virial_yz{0};
  double virial_zx{0}, virial_zy{0}, virial_zz{0};
};

class Particles {
public:
  int N{};
  alignas(64) std::vector<double, AlignedAllocator<double, 64>> mass_{}; // Particle
                                                                         // mass

  // Position (coordinate)
  alignas(64) std::vector<double, AlignedAllocator<double, 64>> coord_x_,
      coord_y_, coord_z_;

  // Velocity
  alignas(64) std::vector<double, AlignedAllocator<double, 64>> velocity_x_{},
      velocity_y_{}, velocity_z_{};

  // Force acting on the particle
  alignas(64) std::vector<double, AlignedAllocator<double, 64>> force_x_{},
      force_y_{}, force_z_{};

  // Virials (3x3 matrix)
  alignas(64) std::vector<double, AlignedAllocator<double, 64>> virial_xx_{},
      virial_xy_{}, virial_xz_{};
  alignas(64) std::vector<double, AlignedAllocator<double, 64>> virial_yx_{},
      virial_yy_{}, virial_yz_{};
  alignas(64) std::vector<double, AlignedAllocator<double, 64>> virial_zx_{},
      virial_zy_{}, virial_zz_{};

  // Energies
  alignas(64) std::vector<double, AlignedAllocator<double, 64>> e_pot_{},
      e_kin_{}, e_int_{}, e_term_{}, e_full_{};

  // Pulse of the particle
  alignas(64) std::vector<double, AlignedAllocator<double, 64>> pulse_{};

  // EAM Properties
  // Electron density of the particle
  alignas(
      64) std::vector<double, AlignedAllocator<double, 64>> electron_density_{};

  // Pair potential energy of the particle
  alignas(
      64) std::vector<double, AlignedAllocator<double, 64>> pair_potential_{};

  // OM of particles
  alignas(64) std::vector<double, AlignedAllocator<double, 64>> ce_{};
  alignas(64) std::vector<double, AlignedAllocator<double, 64>> d_ce_{};

public:
  // Default constructor: Initializes everything to zero.
  Particles() = default;
  ~Particles() = default;

  Particles(int n) : N(n) {
    mass_.resize(N);
    coord_x_.resize(N);
    coord_y_.resize(N);
    coord_z_.resize(N);
    velocity_x_.resize(N);
    velocity_y_.resize(N);
    velocity_z_.resize(N);
    force_x_.resize(N);
    force_y_.resize(N);
    force_z_.resize(N);

    virial_xx_.resize(N);
    virial_xy_.resize(N);
    virial_xz_.resize(N);
    virial_yx_.resize(N);
    virial_yy_.resize(N);
    virial_yz_.resize(N);
    virial_zx_.resize(N);
    virial_zy_.resize(N);
    virial_zz_.resize(N);

    e_kin_.resize(N), e_pot_.resize(N), e_int_.resize(N), e_term_.resize(N),
        e_full_.resize(N);

    pulse_.resize(N);

    electron_density_.resize(N);
    pair_potential_.resize(N);
    ce_.resize(N);
    d_ce_.resize(N);
  };

  void resize(int n) {
    std::cout << "Calling resize with n = " << n << std::endl;
    N = n;
    mass_.resize(n);
    coord_x_.resize(n);
    coord_y_.resize(n);
    coord_z_.resize(n);
    velocity_x_.resize(n);
    velocity_y_.resize(n);
    velocity_z_.resize(n);
    force_x_.resize(n);
    force_y_.resize(n);
    force_z_.resize(n);

    virial_xx_.resize(n);
    virial_xy_.resize(n);
    virial_xz_.resize(n);
    virial_yx_.resize(n);
    virial_yy_.resize(n);
    virial_yz_.resize(n);
    virial_zx_.resize(n);
    virial_zy_.resize(n);
    virial_zz_.resize(n);

    e_kin_.resize(n), e_pot_.resize(n), e_int_.resize(n), e_term_.resize(n),
        e_full_.resize(n);
    pulse_.resize(n);

    electron_density_.resize(n);
    pair_potential_.resize(n);
    ce_.resize(N);
    d_ce_.resize(N);
  }

  inline void applyForceInteraction(const int index,
                                    const ForceCalcValues &result) {
    force_x_[index] = result.force_x;
    force_y_[index] = result.force_y;
    force_z_[index] = result.force_z;

    virial_xx_[index] = result.virial_xx;
    virial_xy_[index] = result.virial_xy;
    virial_xz_[index] = result.virial_xz;

    virial_yx_[index] = result.virial_yx;
    virial_yy_[index] = result.virial_yy;
    virial_yz_[index] = result.virial_yz;

    virial_zx_[index] = result.virial_zx;
    virial_zy_[index] = result.virial_zy;
    virial_zz_[index] = result.virial_zz;

    e_pot_[index] = result.e_pot;
  }

  inline void updateEnergy(double &vcm_x, double &vcm_y, double vcm_z) {
    double halfMass = 0.0;
    double velocitySquared = 0.0;
    double relativeVelocitySquared = 0.0;

    double e_kin = 0.0;
    double e_term = 0.0;

    for (int i = 0; i < N; i++) {
      halfMass = mass_[i] * 0.5;
      velocitySquared = (velocity_x_[i] * velocity_x_[i]) +
                        (velocity_y_[i] * velocity_y_[i]) +
                        (velocity_z_[i] * velocity_z_[i]);
      relativeVelocitySquared =
          ((velocity_x_[i] - vcm_x) * (velocity_x_[i] - vcm_x)) +
          ((velocity_y_[i] - vcm_y) * (velocity_y_[i] - vcm_y)) +
          ((velocity_z_[i] - vcm_z) * (velocity_z_[i] - vcm_z));
      e_kin = halfMass * velocitySquared;
      e_term = halfMass * relativeVelocitySquared;

      e_kin_[i] = e_kin;
      e_term_[i] = e_term;
      e_int_[i] = e_pot_[i] + e_term;
      e_full_[i] = e_pot_[i] + e_kin;
    }
  }

  inline void updatePulse() {
    double velocitySquared = 0.0;
    for (int i = 0; i < N; i++) {
      velocitySquared = (velocity_x_[i] * velocity_x_[i]) +
                        (velocity_y_[i] * velocity_y_[i]) +
                        (velocity_z_[i] * velocity_z_[i]);
      pulse_[i] = mass_[i] * velocitySquared;
    }
  }

  constexpr int size() const { return N; }

  // Запрещаем копирование
  Particles(const Particles &) = delete;
  Particles &operator=(const Particles &) = delete;
};
#endif // PARTICLE_H
