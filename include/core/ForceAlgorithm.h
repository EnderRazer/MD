/**
 * @file ForceAlgorithm.h
 * @brief Defines the ForceAlgorithm class for calculating forces between
 * particles.
 */

#ifndef FORCE_ALGORITHM_H
#define FORCE_ALGORITHM_H

#include "classes/Dimensions.h"
#include "classes/Matrix3.h"
#include "classes/Particle.h"
#include "core/Settings.h"
#include "core/System.h"
#include "potentials/Potential.h"

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
  std::unique_ptr<Potential> potential_; /**< Pointer to the potential model */

  /**
   * @brief Applies mirroring conditions to a vector.
   *
   * Adjusts vector components to account for periodic boundary conditions
   * in all three dimensions.
   *
   * @param vec Vector to be adjusted
   * @param lx System length in x dimension
   * @param ly System length in y dimension
   * @param lz System length in z dimension
   * @return Adjusted vector
   */
  inline Vector3<double> mirror_vector(Vector3<double> &vec, double lx,
                                       double ly, double lz) {
    double halfLX = lx / 2;
    double halfLY = ly / 2;
    double halfLZ = lz / 2;

    if (vec.x() > halfLX)
      vec.x() -= lx;
    if (vec.x() <= -halfLX)
      vec.x() += lx;

    if (vec.y() > halfLY)
      vec.y() -= ly;
    if (vec.y() <= -halfLY)
      vec.y() += ly;

    if (vec.z() > halfLZ)
      vec.z() -= lz;
    if (vec.z() <= -halfLZ)
      vec.z() += lz;
    return vec;
  }

public:
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
  inline double PotentialEnergy_EAM(Particle &p1) {
    return potential_->getU(p1.electron_density(), p1.pairPotential());
  }
  inline ForceCalcValues compute(Particle &p1, Particle &p2) {
    ForceCalcValues output;
    Vector3<double> rVec = p1.coord() - p2.coord();
    // Отражение частицы
    if (settings_.hasPbc()) {
      double lx = dim_.lx();
      double ly = dim_.ly();
      double lz = dim_.lz();
      rVec = mirror_vector(rVec, lx, ly, lz);
    }
    output.rVec = rVec;
    double lengthSqr = rVec.lengthSquared();
    // Проверка на соседство
    if (lengthSqr > potential_->getSqrRcut())
      return output;
    output.interaction_count++;
    double length = rVec.length();
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
      FU = potential_->getFU(p1.electron_density(), p2.electron_density(),
                             length);
      break;
    }
    default: {
      throw std::runtime_error("Unknown type of potential");
      break;
    }
    }

    output.e_pot = U;
    output.force = rVec * FU / length;
    output.virials = Matrix3::outerProduct(rVec, output.force);

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
  inline void preCompute(Particle &p1, Particle &p2) {
    Vector3<double> rVec = p1.coord() - p2.coord();
    // Отражение частицы
    if (settings_.hasPbc()) {
      double lx = dim_.lx();
      double ly = dim_.ly();
      double lz = dim_.lz();
      rVec = mirror_vector(rVec, lx, ly, lz);
    }
    double lengthSqr = rVec.lengthSquared();
    // Проверка на соседство
    if (lengthSqr > potential_->getSqrRcut())
      return;
    double length = rVec.length();

    double mu = potential_->getPairPart(length);
    double rho_f = potential_->getDensityPart(length);

    p1.addElectronDensity(rho_f);
    p1.addPairPotential(mu);
  }
  /**
   * @brief Returns the cutoff distance for the potential.
   * @return Cutoff distance
   */
  inline const double getCutOff() const { return potential_->getRcut(); }

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
  /** @brief Default constructor is deleted */
  inline ForceAlgorithm() = delete;
  /** @brief Default destructor */
  ~ForceAlgorithm() = default;
  /** @brief Copy constructor is deleted */
  ForceAlgorithm(const ForceAlgorithm &) = delete;
  /** @brief Assignment operator is deleted */
  ForceAlgorithm &operator=(const ForceAlgorithm &) = delete;
};
#endif // FORCE_ALGORITHM_H
