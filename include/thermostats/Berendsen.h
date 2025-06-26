//===================================================================================
// Реализация термостата Берендсена
//===================================================================================
#ifndef TSBERENDSEN_H
#define TSBERENDSEN_H

#include <sstream>

#include "nlohmann/json.hpp"

#include "classes/Particle.h"
#include "classes/Vector3.h"
#include "core/Settings.h"
#include "core/System.h"

#include "Thermostat.h"

using json = nlohmann::json;


class ThermostatBerendsen : public Thermostat {

private:

  /**
   * @brief Тип термостата.
   */
  const ThermostatType type_{ThermostatType::BERENDSEN};

  /**
   * @brief Флаг включения термостата.
   */
  bool toggle{false};

  /**
   * @brief Предпочтительная температура.
   */
  double pref_temperature{0.0};

  /**
   * @brief Время взаимодействия с резервуаром.
   */
  double tau_t{0.0};

  /**
   * @brief Коэффициент для термостата Берендсена.
   */
  double dt_over_tau{0.0};

public:
  /**
   * @brief Конструктор.
   * @param config - конфигурация.
   * @param settings - настройки.
   */
  ThermostatBerendsen(json &config, Settings &settings) {
    toggle = config.value("toggle", false);
    pref_temperature = config.value("pref_temp", 0.0);
    tau_t = config.value("tau", 0.0);

    dt_over_tau = settings.dt() / tau_t;
  };

  /**
   * @brief Деструктор.
   */
  ~ThermostatBerendsen() = default;

  /**
   * @brief Получение типа термостата.
   * @return Тип термостата.
   */
  inline const ThermostatType getThermostatType() const override {
    return type_;
  }

  /**
   * @brief Применение температурного контроля.
   * @param sys - система частиц.
   */
  void applyTemperatureControl(System &sys) override {
    double lambda =
        sqrt(1 + dt_over_tau * ((pref_temperature / sys.temperature()) - 1));
    for (auto &particle : sys.particles()) {
      Vector3<double> new_velocity = particle.velocity() * lambda;
      particle.setVelocity(new_velocity);
    }
  }

  /**
   * @brief Получение базовой информации о термостате.
   * @return Базовая информация о термостате.
   */
  std::string getData() const override {
    std::ostringstream oss;
    oss.precision(16);
    oss << "Using Berendsen thermostat!"
        << "\n\tToggled: " << toggle
        << "\n\tPrefered temperature: " << pref_temperature
        << "\n\tTau_t: " << tau_t << "\n";
    return oss.str();
  }
  

  // Запрещаем копирование
  ThermostatBerendsen(const ThermostatBerendsen &) = delete;
  ThermostatBerendsen &operator=(const ThermostatBerendsen &) = delete;
};
#endif
