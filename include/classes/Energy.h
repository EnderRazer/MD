#ifndef ENERGIES_H
#define ENERGIES_H

struct Energy {
    // Define energy types.
    enum EnergyType { Full, Potential, Kinetic, Thermodynamic, Internal, Count };
    // Arrays for current energy values and their averages.
    std::array<double, static_cast<size_t>(Count)> values{}; // Current energies

    // Setter for the energy value.
    inline void set(EnergyType type, double value) { values[static_cast<size_t>(type)] = value; }
    inline void add(EnergyType type, double value) { values[static_cast<size_t>(type)] += value; }
    // Getter for the energy value.
    inline double get(EnergyType type) const { return values[static_cast<size_t>(type)]; }

    // Overload operator<< to output the energy values.
    inline friend std::ostream &operator<<(std::ostream &os, const Energy &energy) {
        os << "Energy values:\n"
           << "Full: " << energy.values[static_cast<size_t>(EnergyType::Full)] << "\n"
           << "Potential: " << energy.values[static_cast<size_t>(EnergyType::Potential)] << "\n"
           << "Kinetic: " << energy.values[static_cast<size_t>(EnergyType::Kinetic)] << "\n"
           << "Thermodynamic: " << energy.values[static_cast<size_t>(EnergyType::Thermodynamic)] << "\n"
           << "Internal: " << energy.values[static_cast<size_t>(EnergyType::Internal)] << "\n";
        return os;
    }
};

#endif // ENERGIES_H