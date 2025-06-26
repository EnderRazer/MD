#ifndef ENERGIES_H
#define ENERGIES_H

#include <array>
#include <ostream>

/**
 * @brief Класс для хранения энергий системы.
 *
 * Класс для хранения энергий системы.
 */
struct Energy {
  /**
   * @brief Типы энергий.
   *
   * Типы энергий.
   */
  enum EnergyType { Full, Potential, Kinetic, Thermodynamic, Internal, Count };

  /**
   * @brief Значения энергий.
   *
   * Значения энергий.
   */
  std::array<double, static_cast<size_t>(Count)> values{};

  /**
   * @brief Установка значения энергии.
   *
   * Установка значения энергии.
   */
  inline void set(EnergyType type, double value) {
    values[static_cast<size_t>(type)] = value;
  }

  /**
   * @brief Добавление значения энергии.
   *
   * Добавление значения энергии.
   */
  inline void add(EnergyType type, double value) {
    values[static_cast<size_t>(type)] += value;
  }

  /**
   * @brief Получение значения энергии.
   *
   * Получение значения энергии.
   */
  inline double get(EnergyType type) const {
    return values[static_cast<size_t>(type)];
  }

  /**
   * @brief Перегрузка оператора << для вывода значений энергий.
   *
   * Перегрузка оператора << для вывода значений энергий.
   */
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
