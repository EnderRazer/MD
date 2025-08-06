//===================================================================================
// Реализация термостата Ланжевена
//===================================================================================
#ifndef THERMOSTAT_LNGVN
#define THERMOSTAT_LNGVN

#include <sstream>
#include <vector>

#include "nlohmann/json.hpp"
#include "random_normalized.h"

#include "classes/Particles.h"
#include "core/Settings.h"
#include "core/System.h"

#include "Thermostat.h"

using json = nlohmann::json;

class ThermostatLangevine : public Thermostat {
public:
  const ThermostatType type_{ThermostatType::LANGEVIN};

private:
  int seed_{0};
  bool toggle_{false};
  double tau_t_{0.0};
  double pref_temperature_{0.0};
  double Vs_x_{0.0}, Vs_y_{0.0}, Vs_z_{0.0}; // Скорость потока

  double sigma_lang_{0.0}; // Константа для термостата Ланжевена

  double eps_{0.0};                              // Случайная величина (0;1)
  std::vector<double> Fr_x_{}, Fr_y_{}, Fr_z_{}; // Случайная сила
  std::vector<double> Ft_x_{}, Ft_y_{}, Ft_z_{}; // Сила трения

public:
  ~ThermostatLangevine() noexcept override = default;
  ThermostatLangevine(json &config, Settings &settings, int particle_number) {
    toggle_ = config.value("toggle", false);
    pref_temperature_ = config.value("pref_temp", 0.0);
    tau_t_ = config.value("tau", 0.0);
    std::vector<double> Vs =
        config["stream_velocity"].get<std::vector<double>>();
    Vs_x_ = Vs[0], Vs_y_ = Vs[1], Vs_z_ = Vs[2];
    sigma_lang_ = sqrt((2 * settings.constants().KBOLTZMAN * pref_temperature_ *
                        settings.mass()) /
                       (tau_t_ * settings.dt()));
    seed_ = settings.seed();
    Fr_x_.resize(particle_number), Fr_y_.resize(particle_number),
        Fr_z_.resize(particle_number);
    Ft_x_.resize(particle_number), Ft_y_.resize(particle_number),
        Ft_z_.resize(particle_number);
  }

  // Запрещаем копирование
  ThermostatLangevine(const ThermostatLangevine &) = delete;
  ThermostatLangevine &operator=(const ThermostatLangevine &) = delete;

  void LGVN_generateFr(System &sys) {
    int pn = sys.particleNumber();
    for (int i = 0; i < pn / 2; i++) {

      eps_ = gasdev(seed_);
      Fr_x_[i] = sigma_lang_ * eps_;
      Fr_x_[pn / 2 + i] = -Fr_x_[i];

      eps_ = gasdev(seed_);
      Fr_y_[i] = sigma_lang_ * eps_;
      Fr_y_[pn / 2 + i] = -Fr_y_[i];

      eps_ = gasdev(seed_);
      Fr_z_[i] = sigma_lang_ * eps_;
      Fr_z_[pn / 2 + i] = -Fr_z_[i];
    }
  };

  void LGVN_generateFt(System &sys) {
    Particles &particles = sys.particles();
    for (int i = 0; i < particles.size(); i++) {
      Ft_x_[i] =
          -particles.mass_[i] * (particles.velocity_x_[i] - Vs_x_) / tau_t_;
      Ft_y_[i] =
          -particles.mass_[i] * (particles.velocity_y_[i] - Vs_y_) / tau_t_;
      Ft_z_[i] =
          -particles.mass_[i] * (particles.velocity_z_[i] - Vs_z_) / tau_t_;
    }
  };
  void applyTemperatureControl(System &sys) override {
    LGVN_generateFr(sys);
    LGVN_generateFt(sys);
    Particles &particles = sys.particles();
    for (int i = 0; i < particles.size(); i++) {
      particles.force_x_[i] += 2 * (Fr_x_[i] + Ft_x_[i]);
      particles.force_y_[i] += 2 * (Fr_y_[i] + Ft_y_[i]);
      particles.force_z_[i] += 2 * (Fr_z_[i] + Ft_z_[i]);
    }
  };

  std::string getData() const override {
    std::ostringstream oss;
    oss.precision(16);
    oss << "Using Langevin thermostat!";
    oss << "\n\tToggled: " << toggle_;
    oss << "\n\tPrefered temperature: " << pref_temperature_;
    oss << "\n\tTau_t: " << tau_t_ << "\n\tStream velocity: (" << Vs_x_ << ", "
        << Vs_y_ << ", " << Vs_z_ << ")";
    oss << "\n\tSIGMA_LANG: " << sigma_lang_ << "\n";
    return oss.str();
  }
  inline const ThermostatType getThermostatType() const override {
    return type_;
  }
};
#endif
