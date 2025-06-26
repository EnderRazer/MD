//===================================================================================
// Реализация термостата Ланжевена
//===================================================================================
#ifndef THERMOSTAT_LNGVN
#define THERMOSTAT_LNGVN

#include <sstream>

#include "nlohmann/json.hpp"
#include "helpers/random_normalized.h"

#include "classes/Particle.h"
#include "classes/Vector3.h"
#include "core/Settings.h"
#include "core/System.h"

#include "Thermostat.h"

using json = nlohmann::json;

/**
 * @brief Класс для термостата Ланжевена.
 * @details Класс предоставляет методы для управления температурой системы.
 */
class ThermostatLangevine : public Thermostat {
private:
/**
   * @brief Тип термостата.
   */
  const ThermostatType type_{ThermostatType::LANGEVIN};

  /**
   * @brief Сид.
   */
  int seed_{0};

  /**
   * @brief Флаг включения термостата.
   */
  bool toggle_{false};

  /**
   * @brief Время взаимодействия с резервуаром.
   */
  double tau_t_{0.0};

  /**
   * @brief Предпочтительная температура.
   */
  double pref_temperature_{0.0};

  /**
   * @brief Скорость потока.
   */
  Vector3<double> Vs_{};

  /**
   * @brief Константа для термостата Ланжевена.
   */
  double sigma_lang_{0.0}; // Константа для термостата Ланжевена

  /**
   * @brief Случайная величина.
   */
  double eps_{0.0};

  /**
   * @brief Случайная сила.
   */
  std::vector<Vector3<double>> Fr_{};

  /**
   * @brief Сила трения.
   */
  std::vector<Vector3<double>> Ft_{};

public:
  /**
   * @brief Конструктор.
   * @param config - конфигурация.
   * @param settings - настройки.
   * @param particle_number - количество частиц.
   */
  ThermostatLangevine(json &config, Settings &settings, int particle_number) {
    toggle_ = config.value("toggle", false);
    pref_temperature_ = config.value("pref_temp", 0.0);
    tau_t_ = config.value("tau", 0.0);
    Vs_ = config["stream_velocity"].get<std::vector<double>>();
    sigma_lang_ = sqrt((2 * settings.constants().KBOLTZMAN * pref_temperature_ *
                        settings.mass()) /
                       (tau_t_ * settings.dt()));
    seed_ = settings.seed();
    Fr_.resize(particle_number);
    Ft_.resize(particle_number);
  }

  /**
   * @brief Деструктор.
   */
  ~ThermostatLangevine() noexcept override = default;

/**
   * @brief Получение типа термостата.
   * @return Тип термостата.
   */
  inline const ThermostatType getThermostatType() const override {
    return type_;
  }
  /**
   * @brief Генерация случайной силы.
   * @param sys - система частиц.
   */
  void LGVN_generateFr(System &sys) {
    int pn = sys.particleNumber();
    for (int i = 0; i < pn / 2; i++) {

      eps_ = gasdev(seed_);
      Fr_[i].x() = sigma_lang_ * eps_;
      Fr_[pn / 2 + i].x() = -Fr_[i].x();

      eps_ = gasdev(seed_);
      Fr_[i].y() = sigma_lang_ * eps_;
      Fr_[pn / 2 + i].y() = -Fr_[i].y();

      eps_ = gasdev(seed_);
      Fr_[i].z() = sigma_lang_ * eps_;
      Fr_[pn / 2 + i].z() = -Fr_[i].z();
    }
  };

  /**
   * @brief Генерация силы трения.
   * @param sys - система частиц.
   */
  void LGVN_generateFt(System &sys) {
    int pn = sys.particleNumber();
    std::vector<Particle> &particles = sys.particles();
    for (int i = 0; i < pn; i++) {
      Ft_[i] =
          -particles[i].mass() * (particles[i].velocity() - Vs_) / tau_t_;
    }
  };

  /**
   * @brief Применение температурного контроля.
   * @param sys - система частиц.
   */
  void applyTemperatureControl(System &sys) override {
    LGVN_generateFr(sys);
    LGVN_generateFt(sys);
    int pn = sys.particleNumber();
    std::vector<Particle> &particles = sys.particles();
    for (int i = 0; i < pn; i++) {
      particles[i].force() += 2 * (Fr_[i] + Ft_[i]);
    }
  };

  /**
   * @brief Получение базовой информации о термостате.
   * @return Базовая информация о термостате.
   */
  std::string getData() const override {
    std::ostringstream oss;
    oss.precision(16);
    oss << "Using Langevin thermostat!"
        << "\n\tToggled: " << toggle_
        << "\n\tPrefered temperature: " << pref_temperature_
        << "\n\tTau_t: " << tau_t_ << "\n\tStream velocity: " << Vs_
        << "\n\tSIGMA_LANG: " << sigma_lang_ << "\n";
    return oss.str();
  }

  // Запрещаем копирование
  ThermostatLangevine(const ThermostatLangevine &) = delete;
  ThermostatLangevine &operator=(const ThermostatLangevine &) = delete;
};
#endif
