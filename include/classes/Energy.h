#ifndef ENERGIES_H
#define ENERGIES_H

#include <array>
#include <ostream>

struct Energy {
  // Define energy types.
  enum EnergyType { Full, Potential, Kinetic, Thermodynamic, Internal, Count };
  // Arrays for current energy values and their averages.
  std::array<double, static_cast<size_t>(Count)> values{}; // Current energies

  // Setter for the energy value.
  inline void set(EnergyType type, double value) {
    values[static_cast<size_t>(type)] = value;
  }
  inline void add(EnergyType type, double value) {
    values[static_cast<size_t>(type)] += value;
  }
  // Getter for the energy value.
  inline double get(EnergyType type) const {
    return values[static_cast<size_t>(type)];
  }

  // Overload operator<< to output the energy values.
  inline friend std::ostream &operator<<(std::ostream &os,
                                         const Energy &energy) {
    os << "Energy values:"
       << "\n\tFull: " << energy.values[EnergyType::Full]
       << "\n\tPotential: " << energy.values[EnergyType::Potential]
       << "\n\tKinetic: " << energy.values[EnergyType::Kinetic]
       << "\n\tThermodynamic: " << energy.values[EnergyType::Thermodynamic]
       << "\n\tInternal: " << energy.values[EnergyType::Internal] << "\n";
    return os;
  }
};

#endif // ENERGIES_H
