#ifndef PARTICLE_H
#define PARTICLE_H

#include "classes/Energy.h"
#include "classes/Matrix3.h"
#include "classes/Vector3.h"

#include <iostream>
#include <ostream>

struct ForceCalcValues {
  int interaction_count{0};
  Vector3<double> rVec{0, 0, 0};
  double e_pot{0.0};
  Vector3<double> force{0, 0, 0};
  Matrix3 virials{};

  friend std::ostream &operator<<(std::ostream &os,
                                  const ForceCalcValues &fcv) {
    os << "\n\tInteractions: " << fcv.interaction_count
       << "\n\trVec: " << fcv.rVec << "\n\tPotential energy: " << fcv.e_pot
       << "\n\tForce: " << fcv.force << "\n\tVirials: " << fcv.virials;
    return os;
  }
};

class Particle {
private:
  int id_{0};        // Particle ID
  double mass_{0.0}; // Particle mass

  Vector3<double> coord_{};    // Position (coordinate)
  Vector3<double> velocity_{}; // Velocity
  Vector3<double> force_{};    // Force acting on the particle

  Matrix3 virials_{}; // Virials (3x3 matrix)

  // Energies
  Energy energies_{}; // Energy values

  double pulse_{0.0}; // Pulse of the particle

  // EAM Properties
  double electron_density_{0.0}; // Electron density of the particle
  double pair_potential_{0.0};   // Pair potential energy of the particle

public:
  // Default constructor: Initializes everything to zero.
  Particle() = default;
  ~Particle() = default;

  // Constructor with initial mass
  Particle(int id, double mass) : id_(id), mass_(mass) {}

  // Constructor with initial mass, position, and velocity
  Particle(double mass, const Vector3<double> &coord,
           const Vector3<double> &velocity)
      : mass_(mass), coord_(coord), velocity_(velocity) {}

  // Constructor to initialize all members including energy
  Particle(double mass, const Vector3<double> &coord,
           const Vector3<double> &velocity, const Vector3<double> &force,
           const Matrix3 &virials, Energy &energy, double pulse)
      : mass_(mass), coord_(coord), velocity_(velocity), force_(force),
        virials_(virials), energies_(energy), pulse_(pulse) {}

  // Getters and setters

  // ID
  inline int getId() const { return id_; }
  inline void setId(int id) { id_ = id; }

  // Mass
  inline double getMass() const { return mass_; }
  inline void setMass(double m) { mass_ = m; }

  // Coord (Position)
  inline const Vector3<double> &coord() const { return coord_; }
  inline double getCoordX() const { return coord_.x(); }
  inline double getCoordY() const { return coord_.y(); }
  inline double getCoordZ() const { return coord_.z(); }

  inline void setCoord(const Vector3<double> &pos) { coord_ = pos; }
  inline void setCoordX(double x) { coord_.x() = x; }
  inline void setCoordY(double y) { coord_.y() = y; }
  inline void setCoordZ(double z) { coord_.z() = z; }

  inline void addCoord(const Vector3<double> &pos) { coord_ += pos; }

  // Velocity
  inline const Vector3<double> &velocity() const { return velocity_; }
  inline double getVelocityX() const { return velocity_.x(); }
  inline double getVelocityY() const { return velocity_.y(); }
  inline double getVelocityZ() const { return velocity_.z(); }

  inline void setVelocity(const Vector3<double> &vel) {
    velocity_ = vel;
    updatePulse();
  }
  inline void setVelocityX(double x) {
    velocity_.x() = x;
    updatePulse();
  }
  inline void setVelocityY(double y) {
    velocity_.y() = y;
    updatePulse();
  }
  inline void setVelocityZ(double z) {
    velocity_.z() = z;
    updatePulse();
  }

  inline void addVelocity(const Vector3<double> &vel) {
    velocity_ += vel;
    updatePulse();
  }

  // Force
  inline Vector3<double> &force() { return force_; }
  inline const Vector3<double> &force() const { return force_; }
  inline double getForceX() const { return force_.x(); }
  inline double getForceY() const { return force_.y(); }
  inline double getForceZ() const { return force_.z(); }

  inline void setForce(const Vector3<double> &f) { force_ = f; }
  inline void setForceX(double x) { force_.x() = x; }
  inline void setForceY(double y) { force_.y() = y; }
  inline void setForceZ(double z) { force_.z() = z; }

  // Virials
  inline const Matrix3 &virials() const { return virials_; }
  inline void setVirials(const Matrix3 &v) { virials_ = v; }

  // Energy
  inline double energy(Energy::EnergyType type) { return energies_.get(type); }
  inline const double energy(Energy::EnergyType type) const {
    return energies_.get(type);
  }

  inline Energy &energies() { return energies_; }
  inline const Energy &energies() const { return energies_; }

  inline void setEnergy(Energy::EnergyType type, double value) {
    energies_.set(type, value);
  }
  inline void setEnergies(const Energy &e) { energies_ = e; }

  // Pulse
  inline const double pulse() const { return pulse_; }
  inline void setPulse(double p) { pulse_ = p; }
  inline void updatePulse() {
    pulse_ = (velocity_.x() + velocity_.y() + velocity_.z()) * mass_;
  }

  inline double const electron_density() const { return electron_density_; }
  inline double const pairPotential() const { return pair_potential_; }

  inline void setElectronDensity(double density) {
    electron_density_ = density;
  }
  inline void setPairPotential(double potential) {
    pair_potential_ = potential;
  }

  inline void addElectronDensity(double density) {
    electron_density_ += density;
  }
  inline void addPairPotential(double potential) {
    pair_potential_ += potential;
  }

  inline void applyForceInteraction(const ForceCalcValues &result) {
    force_ = result.force;
    virials_ = result.virials;
    energies_.set(Energy::EnergyType::Potential, result.e_pot);
  }
  inline void updateEnergy(const Vector3<double> &vcm) noexcept {
    const double velocitySquared = velocity_.lengthSquared();
    const double relativeVelocitySquared = (velocity_ - vcm).lengthSquared();

    const double kineticEnergy = 0.5 * mass_ * velocitySquared;
    const double thermoEnergy = 0.5 * mass_ * relativeVelocitySquared;
    const double potentialEnergy = energies_.get(Energy::EnergyType::Potential);

    energies_.set(Energy::EnergyType::Kinetic, kineticEnergy);
    energies_.set(Energy::EnergyType::Thermodynamic, thermoEnergy);
    energies_.set(Energy::EnergyType::Internal, potentialEnergy + thermoEnergy);
    energies_.set(Energy::EnergyType::Full, potentialEnergy + kineticEnergy);
  }
  // Overload operator<< to output the particle's data
  friend std::ostream &operator<<(std::ostream &os, const Particle &particle) {
    os << "Particle data:\n"
       << "Mass: " << particle.mass_ << "\n"
       << "Position: " << particle.coord_ << "\n"
       << "Velocity: " << particle.velocity_ << "\n"
       << "Force: " << particle.force_ << "\n"
       << "Virials: " << particle.virials_ << "\n"
       << "Energy: " << particle.energies_ << "\n"
       << "Pulse: " << particle.pulse_ << "\n";
    return os;
  }
};
#endif // PARTICLE_H
