#ifndef SYSTEM_H
#define SYSTEM_H

#include <cstddef>
#include <deque>
#include <iostream>
#include <nlohmann/json.hpp>
#include <sstream>
#include <string>

#include "classes/Dimensions.h" //Размеры куба
#include "classes/Particles.h"   //Класс частицы

#include "core/Settings.h" //Класс параметров

using json = nlohmann::json;

class System {
private:
  int current_step_{0}; // Текущий шаг моделирования

  Dimensions dimensions_{}; // Размеры системы
  Particles particles_{}; // Частицы

  // Energies
  double e_pot_{},e_kin_{},e_int_{},e_term_{},e_full_{}; // Усредненные энергии системы на 1 частицу
  double e_pot_avg_{},e_kin_avg_{},e_int_avg_{},e_term_avg_{},e_full_avg_{}; // Усредненные энергии системы по шагам

  int window_size_ = 1000;

  //Температура
  double temperature_{0.0}; // Температура моментальная
  double temperature_avg_{0.0}; // Температура усредненная по шагам
  std::deque<double> temperature_window_{}; //Температура, усредненная скользящим окном
  
  //Давление
  double pressure_{0.0}; // Давление моментальное
  double pressure_avg_{0.0}; // Давление усредненное по шагам
  std::deque<double> pressure_window_{}; //Давление, усредненная скользящим окном
 
  //Тензоры давления
  double pressure_xx_{}, pressure_xy_{}, pressure_xz_{};
  double pressure_yx_{}, pressure_yy_{}, pressure_yz_{};
  double pressure_zx_{}, pressure_zy_{}, pressure_zz_{};

  double vcm_x_{},vcm_y_{},vcm_z_{}; // Скорость центра масс

  double pulse_{0.0},
      pulse_avg_{0.0}; // Импульс системы моментальный/усредненный по шагам

public:
  // Конструкторы / Деструкторы
  System(const json &config, const Settings &settings)
      : dimensions_(config.at("dimensions")) {
    int particle_number = 0;
    if (settings.structType() == std::string("CP")) {
      if (settings.hasPbc()) {
        particle_number = dimensions_.numCristX() * dimensions_.numCristY() *
                           dimensions_.numCristZ();
      } else {
        particle_number = (dimensions_.numCristX() + 1) *
                           (dimensions_.numCristY() + 1) *
                           (dimensions_.numCristZ() + 1);
      }
    } else if (settings.structType() == std::string("FCC")) {
      if (settings.hasPbc()) {
        particle_number = 4 * dimensions_.numCristX() *
                           dimensions_.numCristY() * dimensions_.numCristZ();
      } else {
        particle_number =
            (dimensions_.numCristX() + 1) * (dimensions_.numCristY() + 1) *
                (dimensions_.numCristZ() + 1) +
            (3 * dimensions_.numCristX() * dimensions_.numCristY() *
             (dimensions_.numCristZ() + 1));
      }
    } else if (settings.structType() == std::string("Custom")) {
      particle_number = settings.customLayout().at("particles").size();
    } else {
      throw std::runtime_error("Unknown structure type: " +
                               settings.structType());
    }
    // Предзаполнение массива частиц
    particles_.resize(particle_number);

    double mass = settings.mass();
    for(int i=0;i<particles_.size();i++){
      particles_.mass_[i] = mass;
    }
  }
  ~System() = default;

  // Запрещаем копирование
  System(const System &) = delete;
  System &operator=(const System &) = delete;

  // Сеттеры и геттеры
  inline int &currentStep() { return current_step_; }
  inline const int &currentStep() const { return current_step_; }

  inline const size_t particleNumber() const { return particles_.size(); }

  inline Dimensions &dimensions() { return dimensions_; }
  inline const Dimensions &dimensions() const { return dimensions_; }

  inline Particles &particles() { return particles_; }
  inline const Particles &particles() const { return particles_; }

  // Energies
  inline const double& ePot() const {return e_pot_; }
  inline const double& eKin() const {return e_kin_; }
  inline const double& eTerm() const {return e_term_; }
  inline const double& eInt() const {return e_int_; }
  inline const double& eFull() const {return e_full_; }

  // Energies average
  inline const double ePotAvg() const {return e_pot_avg_/(current_step_+1); }
  inline const double eKinAvg() const {return e_kin_avg_/(current_step_+1); }
  inline const double eTermAvg() const {return e_term_avg_/(current_step_+1); }
  inline const double eIntAvg() const {return e_int_avg_/(current_step_+1); }
  inline const double eFullAvg() const {return e_full_avg_/(current_step_+1); }

  inline const double& temperature() const { return temperature_; }
  inline const double temperatureAvg() const { return temperature_avg_/(current_step_+1); }
  inline const double temperatureWindow() const { 
    int size = temperature_window_.size();
    double temp_window = 0.0;
    for(int i = 0; i < size;i++){
      temp_window+=temperature_window_[i];
    }
    return temp_window/size;
  }
  inline void setTemperature(double temp) { 
    temperature_ = temp; 
    temperature_window_.push_back(temp);
    if(temperature_window_.size() > window_size_)
      temperature_window_.pop_front();
    temperature_avg_+=temp;
  }
    
  inline double& pressureXX() { return pressure_xx_; }
  inline const double& pressureXX() const { return pressure_xx_; }
  inline double& pressureXY() { return pressure_xy_; }
  inline const double& pressureXY() const { return pressure_xy_; }
  inline double& pressureXZ() { return pressure_xz_; }
  inline const double& pressureXZ() const { return pressure_xz_; }

  inline double& pressureYX() { return pressure_yx_; }
  inline const double& pressureYX() const { return pressure_yx_; }
  inline double& pressureYY() { return pressure_yy_; }
  inline const double& pressureYY() const { return pressure_yy_; }
  inline double& pressureYZ() { return pressure_yz_; }
  inline const double& pressureYZ() const { return pressure_yz_; }

  inline double& pressureZX() { return pressure_zx_; }
  inline const double& pressureZX() const { return pressure_zx_; }
  inline double& pressureZY() { return pressure_zy_; }
  inline const double& pressureZY() const { return pressure_zy_; }
  inline double& pressureZZ() { return pressure_zz_; }
  inline const double& pressureZZ() const { return pressure_zz_; }

  inline const double& pressure() const { return pressure_; }
  inline const double pressureAvg() const { return pressure_avg_/(current_step_+1); }
  inline const double pressureWindow() const { 
    int size = pressure_window_.size();
    double press_window = 0.0;
    for(int i = 0; i < size;i++){
      press_window+=pressure_window_[i];
    }
    return press_window/size;
  }
  inline void setPressure(double press) { 
    pressure_ = press; 
    pressure_avg_+=press;
    pressure_window_.push_back(press);
    if(pressure_window_.size() > window_size_)
      pressure_window_.pop_front();
  }

  inline double& vcmX() { return vcm_x_; }
  inline const double& vcmX() const { return vcm_x_; }
  inline double& vcmY() { return vcm_y_; }
  inline const double& vcmY() const { return vcm_y_; }
  inline double& vcmZ() { return vcm_z_; }
  inline const double& vcmZ() const { return vcm_z_; }

  inline double& pulse() { return pulse_; }
  inline const double& pulse() const { return pulse_; }

  inline double& pulseAvg() { return pulse_avg_; }
  inline const double& pulseAvg() const { return pulse_avg_; }

  

  // Методы
  void updateEnergy() {
    int size = particles_.size();
    particles_.updateEnergy(vcm_x_,vcm_y_,vcm_z_);
    double e_pot = 0, e_kin = 0, e_term = 0, e_int = 0, e_full = 0;
    for(int i=0;i<size;i++){
      e_pot += particles_.e_pot_[i];
      e_kin += particles_.e_kin_[i];
      e_term += particles_.e_term_[i];
      e_int += particles_.e_int_[i];
      e_full += particles_.e_full_[i];
    }
    // Средняя на 1 частицу
    e_pot_ = e_pot/size;
    e_kin_ = e_kin/size;
    e_term_ = e_term/size;
    e_int_ = e_int/size;
    e_full_ = e_full/size;

    // Средняя по шагам
    e_pot_avg_ += e_pot_;
    e_kin_avg_ += e_kin_;
    e_term_avg_ += e_term_;
    e_int_avg_ += e_int_;
    e_full_avg_ += e_full_;
  }
  void updateVCM() {
    int size = particles_.size();
    double vel_x{0}, vel_y{0}, vel_z{0};
    for(int i=0;i<size;i++){
      vel_x += particles_.velocity_x_[i];
      vel_y += particles_.velocity_y_[i];
      vel_z += particles_.velocity_z_[i];
    }
    vcm_x_ = vel_x / size;
    vcm_y_ = vel_y / size;
    vcm_z_ = vel_z / size;
  }
  void updatePulse() {
    particles_.updatePulse();
    double pulse = 0.0;
    int size = particles_.size();
    for(int i=0;i<size;i++){
      pulse += particles_.pulse_[i];
    }
    pulse_ = pulse / size;
  }

  inline void advanceStep() { current_step_++; }

  // Output methods
  std::string getData() const {
    std::ostringstream oss;
    oss.precision(16);
    oss << dimensions_.getData() << "\nParticle number: " << particles_.size()
        << "\nCurrent step: " << current_step_
        << "\nTemperature: " << temperature_ << "\nPressure: " << pressure_
        << "\n";

    return oss.str();
  }
  std::string getShortData() const {
    std::ostringstream oss;
    oss.precision(16);
    oss << "Current step: " << current_step_ << 
        "\nMomentum:"
          << "\n\tTemperature: " << temperature_ 
          << "\n\tPressure: " << pressure_
        << "\nAverage:"
          << "\n\tTemperature: " << temperatureAvg()
          << "\n\tPressure: " << pressureAvg() << "\n"
        << "Window (Last "<<window_size_<<" steps):"
          << "\n\tTemperature: " << temperatureWindow()
          << "\n\tPressure: " << pressureWindow() << "\n";

    return oss.str();
  }
  
  std::string getParticlesCoordsInfo() const {
    std::ostringstream oss;
    oss.precision(16);
    oss << "Coordinates:\n" << std::endl;
    for(int i=0;i<particles_.size();i++) {
      oss << i << ": "<< particles_.coord_x_[i] << "; " << particles_.coord_y_[i] << "; " << particles_.coord_z_[i] << std::endl;
    }

    return oss.str();
  }
  std::string getParticlesVelocityInfo() const {
    std::ostringstream oss;
    oss.precision(16);
    oss << "Velocity:\n" << std::endl;
    for(int i=0;i<particles_.size();i++) {
      oss << i << ": "<< particles_.velocity_x_[i] << "; " << particles_.velocity_y_[i] << "; " << particles_.velocity_z_[i] << std::endl;
    }

    return oss.str();
  }
  std::string getParticlesForceInfo() const {
    std::ostringstream oss;
    oss.precision(16);
    oss << "Force:\n" << std::endl;
    for(int i=0;i<particles_.size();i++) {
      oss << i << ": "<< particles_.force_x_[i] << "; " << particles_.force_y_[i] << "; " << particles_.force_z_[i] << std::endl;
    }

    return oss.str();
  }
};

#endif // SYSTEM_H
