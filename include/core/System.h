#ifndef SYSTEM_H
#define SYSTEM_H

#include <iostream>
#include <nlohmann/json.hpp>
#include <sstream>
#include <string>
#include <deque>

#include "classes/Dimensions.h" //Размеры куба
#include "classes/Energy.h"     //Энергии
#include "classes/Particle.h"   //Класс частицы

#include "core/Settings.h" //Класс параметров
//#include "classes/Structure.h" //Класс структуры

using json = nlohmann::json;

/**
 * @brief Класс для работы с системой частиц.
 * @details Класс предоставляет методы для работы с системой частиц.
 */
class System {
private:
  int current_step_{0}; // Текущий шаг моделирования
  int particle_number_{0}; // Количество частиц

  Dimensions dimensions_{};           // Размеры системы
  std::vector<Particle> particles_{}; // Частицы
  //std::vector<Structure> structures_{}; // Структуры

  Energy energies_{};     // Усредненные энергии системы на 1 частицу
  Energy energies_avg_{}; // Усредненные энергии системы по шагам

  int window_size_{1000}; // Размер окна температуры/давления

  double temperature_{0.0}; //Температура моментальная
  std::deque<double> temperature_window_{}; // Окно температуры
  double temperature_avg_{0.0}; // Температура усредненная по шагам
  
  double pressure_{0.0}; // Давление моментальное
  std::deque<double> pressure_window_{}; // Окно давления
  double pressure_avg_{0.0}; // Давление усредненное по шагам
  
  Matrix3 pressure_tensors_{};

  Vector3<double> vcm_{}; // Скорость центра масс
  double pulse_{0.0},
      pulse_avg_{0.0}; // Импульс системы моментальный/усредненный по шагам

  inline void initializeParticles(double mass) {
    particles_.reserve(particle_number_);
    for (int i = 0; i < particle_number_; ++i) {
      particles_.emplace_back(i, mass);
    }
  }

public:
  System() = default;
  ~System() = default;

  // Запрещаем копирование
  System(const System &) = delete;
  System &operator=(const System &) = delete;

  /**
   * @brief Конструктор.
   * @param config - конфигурация.
   * @param settings - настройки.
   */
  System(const json &config, const Settings &settings)
      : dimensions_(config.at("dimensions")) {
    if (settings.structType() == std::string("CP")) {
      if (settings.hasPbc()) {
        particle_number_ = dimensions_.numCristX() * dimensions_.numCristY() *
                           dimensions_.numCristZ();
      } else {
        particle_number_ = (dimensions_.numCristX() + 1) *
                           (dimensions_.numCristY() + 1) *
                           (dimensions_.numCristZ() + 1);
      }
    } else if (settings.structType() == std::string("FCC")) {
      if (settings.hasPbc()) {
        particle_number_ = 4 * dimensions_.numCristX() *
                           dimensions_.numCristY() * dimensions_.numCristZ();
      } else {
        particle_number_ =
            (dimensions_.numCristX() + 1) * (dimensions_.numCristY() + 1) *
                (dimensions_.numCristZ() + 1) +
            (3 * dimensions_.numCristX() * dimensions_.numCristY() *
             (dimensions_.numCristZ() + 1));
      }
    } else if (settings.structType() == std::string("Custom")) {
      particle_number_ = settings.customLayout().at("particles").size();
    } else {
      throw std::runtime_error("Unknown structure type: " +
                               settings.structType());
    }
    // Предзаполнение массива частиц
    initializeParticles(settings.mass());
  }

  // Сеттеры и геттеры
  inline int currentStep() const { return current_step_; }
  inline void setCurrentStep(int step) { current_step_ = step; }

  inline int particleNumber() const {
    return particle_number_;
  }

  inline Dimensions &dimensions() { return dimensions_; }
  inline const Dimensions &dimensions() const { return dimensions_; }
  inline void setDimensions(const Dimensions &dims) { dimensions_ = dims; }

  inline std::vector<Particle> &particles() { return particles_; }
  inline const std::vector<Particle> &particles() const { return particles_; }
  inline void setParticles(std::vector<Particle> parts) {
    particles_ = std::move(parts);
  }

  inline Energy &energies() { return energies_; }
  inline const Energy &energies() const { return energies_; }
  inline void setEnergies(const Energy &energy) { energies_ = energy; }

  inline Energy &energiesAvg() { return energies_avg_; }
  inline const Energy &energiesAvg() const { return energies_avg_; }
  inline void setEnergiesAvg(const Energy &energy_avg) {
    energies_avg_ = energy_avg;
  }

  inline double temperature() const { return temperature_; }
  inline double temperatureWindow() const {
    int size = temperature_window_.size();
    double temp_window = 0.0;
    for(int i = 0; i < size;i++){
      temp_window+=temperature_window_[i];
    }
    return temp_window/size;
  }
  inline double temperatureAvg() const { return temperature_avg_/current_step_; }

  inline void setTemperature(double temp) {
    // Моментальная температура
    temperature_ = temp;
    // Оконное усреднение температуры (за n шагов)
    temperature_window_.push_back(temp);
    if (temperature_window_.size() > window_size_) {
      temperature_window_.pop_front();
    }
    // Усредненная температура по шагам
    temperature_avg_ += temp;
  }

  inline Matrix3 pressureTensors() { return pressure_tensors_; };
  inline const Matrix3 pressureTensors() const { return pressure_tensors_; };
  inline void setPressureTensors(Matrix3 p_t) { pressure_tensors_ = p_t; };

  inline double pressure() const { return pressure_; }
  inline double pressureWindow() const {
    int size = pressure_window_.size();
    double press_window = 0.0;
    for(int i = 0; i < size;i++){
      press_window+=pressure_window_[i];
    }
    return press_window/size;}
  inline double pressureAvg() const { return pressure_avg_/current_step_; }

  inline void setPressure(double pres) { 
    // Моментальное давление
    pressure_ = pres;

     // Оконное усреднение давления (за n шагов)
    pressure_window_.push_back(pres);
    if (pressure_window_.size() > window_size_) {
      pressure_window_.pop_front();
    }
    // Усредненное давление по шагам
    pressure_avg_ += pres;
  }

  inline const Vector3<double> &vcm() const { return vcm_; }
  inline void setVcm(const Vector3<double> &v) { vcm_ = v; }

  inline double pulse() const { return pulse_; }
  inline void setPulse(double p) { pulse_ = p; }

  inline double pulseAvg() const { return pulse_avg_; }
  inline void setPulseAvg(double pulse_avg) { pulse_avg_ = pulse_avg; }
  inline void updatePulseAvg() {
    pulse_avg_ = ((pulse_avg_ * (current_step_ - 1)) + pulse_) / current_step_;
  }

  /**
   * @brief Обновление энергий.
   * @details Обновление энергий системы.
   */
  void updateEnergy() {
    double e_pot = 0, e_kin = 0, e_term = 0, e_int = 0, e_full = 0;
    for (auto &p : particles_) {
      p.updateEnergy(vcm_);
      e_pot += p.energy(Energy::EnergyType::Potential);
      e_kin += p.energy(Energy::EnergyType::Kinetic);
      e_term += p.energy(Energy::EnergyType::Thermodynamic);
      e_int += p.energy(Energy::EnergyType::Internal);
      e_full += p.energy(Energy::EnergyType::Full);
    }
    // Средняя на 1 частицу
    energies_.set(Energy::EnergyType::Potential, e_pot / particle_number_);
    energies_.set(Energy::EnergyType::Kinetic, e_kin / particle_number_);
    energies_.set(Energy::EnergyType::Thermodynamic, e_term / particle_number_);
    energies_.set(Energy::EnergyType::Internal, e_int / particle_number_);
    energies_.set(Energy::EnergyType::Full, e_full / particle_number_);
  }

  /**
   * @brief Обновление усредненных энергий.
   * @details Обновление усредненных энергий системы.
   */
  void updateEnergyAvg() {
    int prev_factor = (current_step_ - 1) / current_step_;
    double curr_factor = 1.0 / current_step_;

    for (Energy::EnergyType type :
         {Energy::EnergyType::Potential, Energy::EnergyType::Kinetic,
          Energy::EnergyType::Thermodynamic, Energy::EnergyType::Internal,
          Energy::EnergyType::Full}) {
      double prev_avg = energies_avg_.get(type);
      double curr_value = energies_.get(type);

      energies_avg_.set(type,
                        prev_avg * prev_factor + curr_value * curr_factor);
    }
  }

  /**
   * @brief Обновление скорости центра масс.
   * @details Обновление скорости центра масс системы.
   */
  void updateVCM() {
    Vector3<double> vel;
    for (Particle &p : particles_) {
      vel += p.velocity();
    }
    vcm_ = vel / particle_number_;
  }

  /**
   * @brief Обновление импульса системы.
   * @details Обновление импульса системы.
   */
  void updatePulse() {
    double pulse = 0.0;
    for (Particle &p : particles_) {
      pulse += p.pulse();
    }
    pulse_ = pulse / particle_number_;
  }

  /**
   * @brief Обновление шага.
   * @details Обновление шага системы.
   */
  inline void advanceStep() { current_step_++; }

  /**
   * @brief Получение базовой информации о системе.
   * @details Получение базовой информации о системе.
   */
  std::string getData() const {
    std::ostringstream oss;
    oss.precision(16);
    oss << dimensions_.getData() << "\nParticle number: " << particle_number_
        << "\nCurrent step: " << current_step_
        << "\nTemperature: " << temperature_ << "\nPressure: " << pressure_
        << "\n";

    return oss.str();
  }

  /**
   * @brief Получение краткой информации о системе. Для вывода в консоль.
   * @details Получение краткой информации о системе.
   */
  std::string getShortData() const {
    std::ostringstream oss;
    oss.precision(16);
    oss << "Current step: " << current_step_ << 
        "\nMomentum:" << "\n\tTemperature: " << temperature_ << "\n\tPressure: " << pressure_ << 
        "\nAverage:" << "\n\tTemperature: " << temperatureAvg() << "\n\tPressure: " << pressureAvg() <<
        "\nWindow average:" << "\n\tTemperature: " << temperatureWindow() << "\n\tPressure: " << pressureWindow()  
        << "\n";

    return oss.str();
  }

  /**
   * @brief Получение информации о частицах.
   * @details Получение информации о частицах.
   */
  std::string getParticlesInfo() const {
    std::ostringstream oss;
    oss.precision(16);
    oss << "Частицы: " << particle_number_ << std::endl;
    for (const auto &particle : particles_) {
      oss << particle << std::endl;
    }

    return oss.str();
  }

  /**
   * @brief Получение информации о координатах частиц.
   * @details Получение информации о координатах частиц.
   */
  std::string getParticlesCoordsInfo() const {
    std::ostringstream oss;
    oss.precision(16);
    oss << "Coordinates:\n" << std::endl;
    for (const auto &particle : particles_) {
      oss << particle.coord() << std::endl;
    }

    return oss.str();
  }

  /**
   * @brief Получение информации о скоростях частиц.
   * @details Получение информации о скоростях частиц.
   */
  std::string getParticlesVelocityInfo() const {
    std::ostringstream oss;
    oss.precision(16);
    oss << "Velocity:\n" << std::endl;
    for (const auto &particle : particles_) {
      oss << particle.velocity() << std::endl;
    }

    return oss.str();
  }
};

#endif // SYSTEM_H
