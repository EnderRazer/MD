/**
 * @file ForceAlgorithm.h
 * @brief Defines the ForceAlgorithm class for calculating forces between
 * particles.
 */

#ifndef FORCE_ALGORITHM_H
#define FORCE_ALGORITHM_H

#include "classes/Dimensions.h"
#include "classes/Particles.h"
#include "core/Settings.h"
#include "core/System.h"
#include "potentials/Potential.h"
#include <cmath>

/**
 * @class ForceAlgorithm
 * @brief Handles calculation of forces between particles using different
 * potential models.
 *
 * This class encapsulates the logic for computing forces, potential energies,
 * and other related values between particles in a molecular dynamics
 * simulation. It supports different potential types including LJ and
 * EAM.
 */
class ForceAlgorithm {
private:
  Settings &settings_; /**< Reference to simulation settings */
  Dimensions &dim_;    /**< Reference to system dimensions */
  

public:
std::unique_ptr<Potential> potential_; /**< Pointer to the potential model */
  /**
   * @brief Computes force and related values between two particles.
   *
   * Calculates the force, potential energy, and virials between two particles
   * based on the configured potential model and system settings.
   *
   * @param p1 First particle
   * @param p2 Second particle
   * @return Calculation results including force, potential energy, etc.
   */
  inline double PotentialEnergy_EAM(Particles &p, int index) {
      return potential_->getU(p.electronDensity(index), p.pairPart(index));
  }
  inline ForceCalcValues compute(Particles &p, int i1, int i2) {
    ForceCalcValues output;
    double rVec_x = p.coordX(i1) - p.coordX(i2);
    double rVec_y = p.coordY(i1) - p.coordY(i2);
    double rVec_z = p.coordZ(i1) - p.coordZ(i2);

    // Отражение частицы
    if (settings_.hasPbc()) {
      double lx = dim_.sizeX();
      double ly = dim_.sizeY();
      double lz = dim_.sizeZ();

      rVec_x -= dim_.sizeX() * std::round(rVec_x / lx);
      rVec_y -= dim_.sizeY() * std::round(rVec_y / ly);
      rVec_z -= dim_.sizeZ() * std::round(rVec_z / lz);
    }
    output.rVec_x = rVec_x;
    output.rVec_y = rVec_y;
    output.rVec_z = rVec_z;

    double lengthSqr = rVec_x*rVec_x + rVec_y*rVec_y + rVec_z*rVec_z;
    // Проверка на соседство
    if (lengthSqr > potential_->getSqrRcut())
      return output;
    output.interaction_count++;
    double length = std::sqrt(lengthSqr);
    // Расчет силы
    double U = 0.0;
    double FU = 0.0;
    switch (potential_->getPotentialType()) {
    case Potential::PotentialType::LJ: {
      PotentialResult res = potential_->getAll(lengthSqr);
      // U = potential_->getU(lengthSqr);
      // FU = potential_->getFU(lengthSqr);
      U = res.u / 2;
      FU = res.fu;
      break;
    }
    case Potential::PotentialType::EAM: {
      // U = potential_->getU(p1.electron_density(), p1.pairPotential());
      FU = potential_->getFU(p.electronDensity(i1), p.electronDensity(i2),
                             length);
      break;
    }
    default: {
      throw std::runtime_error("Unknown type of potential");
      break;
    }
    }

    output.e_pot = U;

    output.force_x = rVec_x * FU / length;
    output.force_y = rVec_y * FU / length;
    output.force_z = rVec_z * FU / length;

    output.virial_xx = rVec_x * output.force_x;
    output.virial_xy = rVec_x * output.force_y;
    output.virial_xz = rVec_x * output.force_z;

    output.virial_yx = rVec_y * output.force_x;
    output.virial_yy = rVec_y * output.force_y;
    output.virial_yz = rVec_y * output.force_z;

    output.virial_zx = rVec_z * output.force_x;
    output.virial_zy = rVec_z * output.force_y;
    output.virial_zz = rVec_z * output.force_z;

    return output;
  }

  /**
   * @brief Performs pre-computation for EAM potential.
   *
   * Calculates and stores intermediate values needed for EAM potential
   * calculations.
   *
   * @param p1 First particle
   * @param p2 Second particle
   */
  inline void preCompute(Particles &p, int i1, int i2) {
    double rVec_x = p.coordX(i1) - p.coordX(i2);
    double rVec_y = p.coordY(i1) - p.coordY(i2);
    double rVec_z = p.coordZ(i1) - p.coordZ(i2);
    
    // Отражение частицы
    if (settings_.hasPbc()) {
      double lx = dim_.sizeX();
      double ly = dim_.sizeY();
      double lz = dim_.sizeZ();

      rVec_x -= dim_.sizeX() * std::round(rVec_x / lx);
      rVec_y -= dim_.sizeY() * std::round(rVec_y / ly);
      rVec_z -= dim_.sizeZ() * std::round(rVec_z / lz);
    }

    double lengthSqr = rVec_x*rVec_x + rVec_y*rVec_y + rVec_z*rVec_z;
    // Проверка на соседство
    if (lengthSqr > potential_->getSqrRcut())
      return;
    double length = std::sqrt(lengthSqr);

    double mu = potential_->getPairPart(length);
    double rho_f = potential_->getDensityPart(length);

    p.electronDensity(i1) += rho_f;
    p.pairPart(i1) += mu;
  }
  /**
   * @brief Returns the cutoff distance for the potential.
   * @return Cutoff distance
   */
  inline const double getCutOff() const { return potential_->getRcut(); }
  inline const double getCutOffSqr() const { return potential_->getSqrRcut(); }

  /**
   * @brief Returns the type of the potential being used.
   * @return Potential type
   */
  inline const Potential::PotentialType getPotentialType() const {
    return potential_->getPotentialType();
  }

  inline const double callCloudCalc(double rho) const {
    return potential_->getCloud(rho);
  }
  inline const double callDerCloudCalc(double rho) const {
    return potential_->getDerCloud(rho);
  }
  /**
   * @brief Constructs a ForceAlgorithm with specified settings and potential.
   *
   * @param settings Reference to simulation settings
   * @param potential Unique pointer to a potential model
   * @param sys Reference to the particle system
   */
  ForceAlgorithm(Settings &settings, std::unique_ptr<Potential> &&potential,
                 System &sys)
      : settings_(settings), potential_(std::move(potential)),
        dim_(sys.dimensions()) {}

  ForceAlgorithm() = delete;
  ~ForceAlgorithm() = default;

  ForceAlgorithm(const ForceAlgorithm &) = delete;
  ForceAlgorithm &operator=(const ForceAlgorithm &) = delete;
};
#endif // FORCE_ALGORITHM_H
