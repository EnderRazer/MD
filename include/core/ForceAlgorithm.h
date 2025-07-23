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
#include "potentials/LJ.h"
#include "potentials/EAM.h"

/**
 * @brief Класс для вычисления сил между частицами.
 * @details Класс предоставляет методы для вычисления сил между частицами и
 * применения периодических граничных условий.
 */
class ForceAlgorithm {
private:
  /**
   * @brief Настройки.
   */
  Settings &settings_;

  /**
   * @brief Размеры системы.
   */
  Dimensions &dim_;

public:
  /**
   * @brief Вычисление LJ потенциала между двумя частицами.
   * @param p1 - первая частица.
   * @param p2 - вторая частица.
   * @return Результат вычисления силы.
   */
  inline ForceCalcValues compute(const LJ &potential,Particle &p1, Particle &p2, bool pbc=false) {
    ForceCalcValues output;
    Vector3<double> rVec = p1.coord() - p2.coord();

    // Отражение частицы
    rVec.x() -= pbc ? dim_.lx() * std::round(rVec.x() / dim_.lx()) : 0.0;
    rVec.y() -= pbc ? dim_.ly() * std::round(rVec.y() / dim_.ly()) : 0.0;
    rVec.z() -= pbc ? dim_.lz() * std::round(rVec.z() / dim_.lz()) : 0.0;
    output.rVec = rVec;

    double lengthSqr = rVec.lengthSquared();
    // Проверка на соседство
    if (lengthSqr > potential.getSqrRcut())
      return output;
    output.interaction_count++;
    double length = rVec.length();
    // Расчет силы
    double U = 0.0;
    double FU = 0.0;

    PotentialResult res = potential.getAll(lengthSqr);
    U = res.u / 2;
    FU = res.fu;

    output.e_pot = U;
    output.force = rVec * FU / length;
    output.virials = Matrix3::outerProduct(rVec, output.force);

    return output;
  }
  /**
   * @brief Вычисление EAM потенциала между двумя частицами.
   * @param p1 - первая частица.
   * @param p2 - вторая частица.
   * @return Результат вычисления потенциала.
   */
  inline ForceCalcValues compute(const EAM &potential,Particle &p1, Particle &p2, bool pbc=false) {
    ForceCalcValues output;
    Vector3<double> rVec = p1.coord() - p2.coord();

    // Отражение частицы
    rVec.x() -= pbc ? dim_.lx() * std::round(rVec.x() / dim_.lx()) : 0.0;
    rVec.y() -= pbc ? dim_.ly() * std::round(rVec.y() / dim_.ly()) : 0.0;
    rVec.z() -= pbc ? dim_.lz() * std::round(rVec.z() / dim_.lz()) : 0.0;
    output.rVec = rVec;

    double lengthSqr = rVec.lengthSquared();
    // Проверка на соседство
    if (lengthSqr > potential.getSqrRcut())
      return output;
    output.interaction_count++;
    double length = rVec.length();
    // Расчет силы
    double U = 0.0;
    double FU = 0.0;

    FU = potential.getFU(p1.electron_density(), p2.electron_density(),length);

    output.e_pot = U;
    output.force = rVec * FU / length;
    output.virials = Matrix3::outerProduct(rVec, output.force);

    return output;
  }

  /**
   * @brief Предвычисление для EAM потенциала.
   * @param p1 - первая частица.
   * @param p2 - вторая частица.
   */
  inline void preComputeEAM(const EAM &potential,Particle &p1, Particle &p2, bool pbc=false) {
    Vector3<double> rVec = p1.coord() - p2.coord();
    
    // Отражение частицы
    rVec.x() -= pbc ? dim_.lx() * std::round(rVec.x() / dim_.lx()) : 0.0;
    rVec.y() -= pbc ? dim_.ly() * std::round(rVec.y() / dim_.ly()) : 0.0;
    rVec.z() -= pbc ? dim_.lz() * std::round(rVec.z() / dim_.lz()) : 0.0;
    
    double lengthSqr = rVec.lengthSquared();
    // Проверка на соседство
    if (lengthSqr > potential.getSqrRcut())
      return;
    double length = rVec.length();

    Mu_RhoF mu_rho_f = potential.getPairDesityPart(length);

    p1.addElectronDensity(mu_rho_f.rho_f);
    p1.addPairPotential(mu_rho_f.mu);
  }

  /**
   * @brief Получение типа потенциала.
   * @return Тип потенциала.
   */
  inline const PotentialType getPotentialType(EAM potential) const {
    return potential.getPotentialType();
  }
  inline const PotentialType getPotentialType(LJ potential) const {
    return potential.getPotentialType();
  }
  /**
   * @brief Конструктор.
   * @param settings - настройки.
   * @param potential - модель потенциала.
   * @param sys - система частиц.
   */
  ForceAlgorithm(Settings &settings,System &sys)
      : settings_(settings),dim_(sys.dimensions()) {}

  /**
   * @brief Конструктор по умолчанию.
   */
  inline ForceAlgorithm() = delete;

  /**
   * @brief Деструктор.
   */
  ~ForceAlgorithm() = default;

  // Запрещаем копирование
  ForceAlgorithm(const ForceAlgorithm &) = delete;
  ForceAlgorithm &operator=(const ForceAlgorithm &) = delete;
};
#endif // FORCE_ALGORITHM_H
