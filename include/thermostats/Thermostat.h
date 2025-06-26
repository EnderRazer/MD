#ifndef THERMOSTAT_H
#define THERMOSTAT_H

#include <string>

#include "core/System.h"

/**
   * @brief Тип термостата.
   */
  enum class ThermostatType { NONE, BERENDSEN, LANGEVIN };
/**
 * @brief Абстрактный класс для термостата.
 * @details Класс предоставляет методы для управления температурой системы.
 */
class Thermostat {
  
private:

  /**
   * @brief Время взаимодействия с резервуаром.
   */
  double tau_t{0.0};

  /**
   * @brief Предпочтительная температура.
   */
  double pref_temperature{0.0};

public:
  

  /**
   * @brief Получение типа термостата.
   * @return Тип термостата.
   */
  inline virtual const ThermostatType getThermostatType() const = 0;

  /**
   * @brief Применение температурного контроля.
   */
  virtual void applyTemperatureControl(System &sys) = 0;

  /**
   * @brief Получение базовой информации о термостате.
   * @return Базовая информация о термостате.
   */
  virtual std::string getData() const = 0;

  /**
   * @brief Деструктор.
   */
  virtual ~Thermostat() = default;
};

#endif
