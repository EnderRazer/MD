#ifndef PARTICLE_H
#define PARTICLE_H

#include "classes/Energy.h"
#include "classes/Matrix3.h"
#include "classes/Vector3.h"

#include <iostream>
#include <ostream>

/**
 * @brief Перечисление типов частиц.
 *
 * Перечисление типов частиц.
 */
enum class ParticleMaterial {
  None,H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,As,Se,Br,Kr,Rb,Sr,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Sb,Te,I,Xe,Cs,Ba,La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi,Po,At,Rn,Fr,Ra,Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr
};

/**
 * @brief Структура для хранения значений силы взаимодействия.
 *
 * Структура для хранения значений силы взаимодействия.
 */
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

/**
 * @brief Класс для хранения информации о частице.
 *
 * Класс для хранения информации о частице.
 */
class Particle {
private:
  /**
   * @brief ID частицы.
   */
  int id_{0};

  /**
   * @brief Тип частицы.
   */
  ParticleMaterial material_{ParticleMaterial::None};      // Particle type

  /**
   * @brief Масса частицы.
   */
  double mass_{0.0};

  /**
   * @brief Координаты частицы.
   */
  Vector3<double> coord_{};

  /**
   * @brief Скорость частицы.
   */
  Vector3<double> velocity_{};

  /**
   * @brief Сила, действующая на частицу.
   */
  Vector3<double> force_{};

  /**
   * @brief Вириалы частицы.
   */
  Matrix3 virials_{};

  /**
   * @brief Энергии частицы.
   */
  Energy energies_{};

  /**
   * @brief Импульс частицы.  
   */
  double pulse_{0.0}; // Pulse of the particle

  /**
   * @brief Электронная плотность частицы.
   */
  double electron_density_{0.0};

  /**
   * @brief Потенциал парного взаимодействия частицы.
   */
  double pair_potential_{0.0};

public:
  /**
   * @brief Конструктор по умолчанию.
   *
   * Конструктор по умолчанию.
   */
  Particle() = default;

  /**
   * @brief Деструктор по умолчанию.
   *
   * Деструктор по умолчанию.
   */
  ~Particle() = default;

  /**
   * @brief Конструктор с начальной массой.
   *
   * Конструктор с начальной массой.
   * @param id - ID частицы.
   * @param mass - масса частицы.
   */
  Particle(int id, double mass) : id_(id), mass_(mass) {} 

  /**
   * @brief Конструктор с начальной массой, координатами и скоростью.
   *
   * Конструктор с начальной массой, координатами и скоростью.
   * @param mass - масса частицы.
   * @param coord - координаты частицы.
   * @param velocity - скорость частицы.
   */
  Particle(int id, double mass, const Vector3<double> &coord,
           const Vector3<double> &velocity)
      : id_(id), mass_(mass), coord_(coord), velocity_(velocity) {}

  /**
   * @brief Конструктор с начальной массой, координатами, скоростью, силой, вириалами, энергией и пульсом.
   *
   * Конструктор с начальной массой, координатами, скоростью, силой, вириалами, энергией и пульсом.
   * @param mass - масса частицы.
   * @param coord - координаты частицы.
   * @param velocity - скорость частицы.
   * @param force - сила, действующая на частицу.
   * @param virials - вириалы частицы.
   * @param energy - энергии частицы.
   * @param pulse - пульс частицы.
   */
  Particle(int id, double mass, const Vector3<double> &coord,
           const Vector3<double> &velocity, const Vector3<double> &force,
           const Matrix3 &virials, Energy &energy, double pulse)
      : id_(id), mass_(mass), coord_(coord), velocity_(velocity), force_(force),
        virials_(virials), energies_(energy), pulse_(pulse) {}

  // Getters and setters

  /**
   * @brief Получение ID частицы.
   * @return ID частицы.
   */
  inline int id() const { return id_; }

  /**
   * @brief Установка ID частицы.
   * @param id - ID частицы.
   */
  inline void setId(int id) { id_ = id; }

  /**
   * @brief Получение массы частицы.
   * @return Масса частицы.
   */
  inline double mass() const { return mass_; }

  /**
   * @brief Установка массы частицы.
   * @param m - масса частицы.
   */
  inline void setMass(double m) { mass_ = m; }

  /**
   * @brief Получение координат частицы.
   * @return Координаты частицы.
   */
  inline const Vector3<double> &coord() const { return coord_; }

  /**
   * @brief Получение X-координаты частицы.
   * @return X-координата частицы.
   */
  inline double getCoordX() const { return coord_.x(); }

  /**
   * @brief Получение Y-координаты частицы.
   * @return Y-координата частицы.
   */
  inline double getCoordY() const { return coord_.y(); }

  /**
   * @brief Получение Z-координаты частицы.
   * @return Z-координата частицы.
   */
  inline double getCoordZ() const { return coord_.z(); }

  /**
   * @brief Установка координат частицы.
   * @param pos - координаты частицы.
   */
  inline void setCoord(const Vector3<double> &pos) { coord_ = pos; }

  /**
   * @brief Установка X-компонента координат частицы.
   * @param x - X-компонента координат частицы.
   */
  inline void setCoordX(double x) { coord_.x() = x; }

  /**
   * @brief Установка Y-компонента координат частицы.
   * @param y - Y-компонента координат частицы.
   */
  inline void setCoordY(double y) { coord_.y() = y; }

  /**
   * @brief Установка Z-компонента координат частицы.
   * @param z - Z-компонента координат частицы.
   */
  inline void setCoordZ(double z) { coord_.z() = z; }

  /**
   * @brief Добавление к координатам частицы.
   * @param pos - координаты частицы.
   */
  inline void addCoord(const Vector3<double> &pos) { coord_ += pos; }

  /**
   * @brief Получение скорости частицы.
   * @return Скорость частицы.
   */
  inline const Vector3<double> &velocity() const { return velocity_; }

  /**
   * @brief Получение X-компоненты скорости частицы.
   * @return X-компонента скорости частицы.
   */
  inline double getVelocityX() const { return velocity_.x(); }

  /**
   * @brief Получение Y-компоненты скорости частицы.
   * @return Y-компонента скорости частицы.
   */
  inline double getVelocityY() const { return velocity_.y(); }

  /**
   * @brief Получение Z-компоненты скорости частицы.
   * @return Z-компонента скорости частицы.
   */
  inline double getVelocityZ() const { return velocity_.z(); }

  /**
   * @brief Установка скорости частицы.
   * @param vel - скорость частицы.
   */
  inline void setVelocity(const Vector3<double> &vel) {
    velocity_ = vel;
    updatePulse();
  }

  /**
   * @brief Установка X-компоненты скорости частицы.
   * @param x - X-компонента скорости частицы.
   */
  inline void setVelocityX(double x) {
    velocity_.x() = x;
    updatePulse();
  }

  /**
   * @brief Установка Y-компоненты скорости частицы.
   * @param y - Y-компонента скорости частицы.
   */
  inline void setVelocityY(double y) {
    velocity_.y() = y;
    updatePulse();
  }

  /**
   * @brief Установка Z-компоненты скорости частицы.
   * @param z - Z-компонента скорости частицы.
   */
  inline void setVelocityZ(double z) {
    velocity_.z() = z;
    updatePulse();
  }

  /**
   * @brief Добавление к скорости частицы.
   * @param vel - скорость частицы.
   */
  inline void addVelocity(const Vector3<double> &vel) {
    velocity_ += vel;
    updatePulse();
  }

  /**
   * @brief Получение силы, действующей на частицу.
   * @return Сила, действующая на частицу.
   */
  inline Vector3<double> &force() { return force_; }

  /**
   * @brief Получение X-компоненты силы, действующей на частицу.
   * @return X-компонента силы, действующей на частицу.
   */
  inline const Vector3<double> &force() const { return force_; }
  
  /**
   * @brief Получение Y-компоненты силы, действующей на частицу.
   * @return Y-компонента силы, действующей на частицу.
   */
  inline double getForceY() const { return force_.y(); }

  /**
   * @brief Получение Z-компоненты силы, действующей на частицу.
   * @return Z-компонента силы, действующей на частицу.
   */
  inline double getForceZ() const { return force_.z(); }

  /**
   * @brief Установка силы, действующей на частицу.
   * @param f - сила, действующая на частицу.
   */
  inline void setForce(const Vector3<double> &f) { force_ = f; }

  /**
   * @brief Установка X-компоненты силы, действующей на частицу.
   * @param x - X-компонента силы, действующей на частицу.
   */
  inline void setForceX(double x) { force_.x() = x; }

  /**
   * @brief Установка Y-компоненты силы, действующей на частицу.
   * @param y - Y-компонента силы, действующей на частицу.
   */
  inline void setForceY(double y) { force_.y() = y; }

  /**
   * @brief Установка Z-компоненты силы, действующей на частицу.
   * @param z - Z-компонента силы, действующей на частицу.
   */
  inline void setForceZ(double z) { force_.z() = z; }

  /**
   * @brief Получение вириалов частицы.
   * @return Вириалы частицы.
   */
  inline const Matrix3 &virials() const { return virials_; }

  /**
   * @brief Установка вириалов частицы.
   * @param v - вириалы частицы.
   */
  inline void setVirials(const Matrix3 &v) { virials_ = v; }

  /**
   * @brief Получение энергии частицы.
   * @return Энергия частицы.
   */
  inline double energy(Energy::EnergyType type) { return energies_.get(type); }

  /**
   * @brief Получение энергии частицы.
   * @param type - тип энергии.
   * @return Энергия частицы.
   */
  inline const double energy(Energy::EnergyType type) const {
    return energies_.get(type);
  }

  /**
   * @brief Получение всех энергий частицы.
   * @return Все энергии частицы.
   */
  inline Energy &energies() { return energies_; }

  /**
   * @brief Получение всех энергий частицы.
   * @return Все энергии частицы.
   */
  inline const Energy &energies() const { return energies_; }

  /**
   * @brief Установка энергии частицы.
   * @param type - тип энергии.
   * @param value - значение энергии.
   */
  inline void setEnergy(Energy::EnergyType type, double value) {
    energies_.set(type, value);
  }

  /**
   * @brief Установка всех энергий частицы.
   * @param e - все энергии частицы.
   */
  inline void setEnergies(const Energy &e) { energies_ = e; }

  /**
   * @brief Получение импульса частицы.
   * @return Импульс частицы.
   */
  inline const double pulse() const { return pulse_; }

  /**
   * @brief Установка импульса частицы.
   * @param p - импульс частицы.
   */
  inline void setPulse(double p) { pulse_ = p; }

  /**
   * @brief Обновление импульса частицы.
   */
  inline void updatePulse() {
    pulse_ = (velocity_.x() + velocity_.y() + velocity_.z()) * mass_;
  }

  /**
   * @brief Применение силы взаимодействия.
   * @param result - результат расчета силы взаимодействия.
   */
  inline void applyForceInteraction(const ForceCalcValues &result) {
    force_ = result.force;
    virials_ = result.virials;
    energies_.set(Energy::EnergyType::Potential, result.e_pot);
  }

  // Energy

  /**
   * @brief Обновление энергии частицы.
   * @param vcm - скорость центра масс.
   */
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

  /**
   * @brief Получение вклада давления частицы.
   * @return Вклад давления частицы.
   */
  inline double pressureContribution(const Vector3<double> &vcm) const {
    double sum_mv = mass_ * (velocity_ - vcm).lengthSquared();
    double sum_vir = 0.5 * (virials_.xx() + virials_.yy() + virials_.zz());

    return sum_mv + sum_vir;
  }

  /**
   * @brief Получение электронной плотности частицы.
   * @return Электронная плотность частицы.
   */
  inline double const electron_density() const { return electron_density_; }

  /**
   * @brief Получение потенциала парного взаимодействия частицы.
   * @return Потенциал парного взаимодействия частицы.
   */
  inline double const pairPotential() const { return pair_potential_; }

  /**
   * @brief Установка электронной плотности частицы.
   * @param density - электронная плотность частицы.
   */
  inline void setElectronDensity(double density) {
    electron_density_ = density;
  }

  /**
   * @brief Установка потенциала парного взаимодействия частицы.
   * @param potential - потенциал парного взаимодействия частицы.
   */
  inline void setPairPotential(double potential) {
    pair_potential_ = potential;
  }

  /**
   * @brief Добавление к электронной плотности частицы.
   * @param density - электронная плотность частицы.
   */
  inline void addElectronDensity(double density) {
    electron_density_ += density;
  }

  /**
   * @brief Добавление к потенциалу парного взаимодействия частицы.
   * @param potential - потенциал парного взаимодействия частицы.
   */
  inline void addPairPotential(double potential) {
    pair_potential_ += potential;
  }
  
  /**
   * @brief Перегрузка оператора << для вывода данных частицы.
   * @return Поток вывода.
   */
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
