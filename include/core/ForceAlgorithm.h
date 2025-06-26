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

  /**
   * @brief Модель потенциала.
   */
  std::unique_ptr<Potential> potential_;

  /**
   * @brief Отражение вектора.
   * @details Отражает вектор, учитывая периодические граничные условия.
   * @param vec - вектор.
   * @param lx - длина системы в x-ой размерности.
   * @param ly - длина системы в y-ой размерности.
   * @param lz - длина системы в z-ой размерности.
   * @return Отраженный вектор.
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
   * @brief Вычисление потенциальной энергии для EAM потенциала.
   * @param p1 - первая частица.
   * @return Потенциальная энергия.
   */
  inline double PotentialEnergy_EAM(Particle &p1) {
    return potential_->getU(p1.electron_density(), p1.pairPotential());
  }

  /**
   * @brief Вычисление силы между двумя частицами.
   * @param p1 - первая частица.
   * @param p2 - вторая частица.
   * @return Результат вычисления силы.
   */
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
    case PotentialType::LJ: {
      PotentialResult res = potential_->getAll(lengthSqr);
      // U = potential_->getU(lengthSqr);
      // FU = potential_->getFU(lengthSqr);
      U = res.u / 2;
      FU = res.fu;
      break;
    }
    case PotentialType::EAM: {
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
   * @brief Предвычисление для EAM потенциала.
   * @param p1 - первая частица.
   * @param p2 - вторая частица.
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

    Mu_RhoF mu_rho_f = potential_->getPairDesityPart(length);

    p1.addElectronDensity(mu_rho_f.rho_f);
    p1.addPairPotential(mu_rho_f.mu);
  }

  /**
   * @brief Получение радиуса обрезания.
   * @return Радиус обрезания.
   */
  inline const double getCutOff() const { return potential_->getRcut(); }

  /**
   * @brief Получение типа потенциала.
   * @return Тип потенциала.
   */
  inline const PotentialType getPotentialType() const {
    return potential_->getPotentialType();
  }

  inline const double callCloudCalc(double rho) const {
    return potential_->getCloud(rho);
  }

  inline const double callDerCloudCalc(double rho) const {
    return potential_->getDerCloud(rho);
  }

  /**
   * @brief Конструктор.
   * @param settings - настройки.
   * @param potential - модель потенциала.
   * @param sys - система частиц.
   */
  ForceAlgorithm(Settings &settings, std::unique_ptr<Potential> &&potential,
                 System &sys)
      : settings_(settings), potential_(std::move(potential)),
        dim_(sys.dimensions()) {}

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
