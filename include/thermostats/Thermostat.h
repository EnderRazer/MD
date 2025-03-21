#ifndef THERMOSTAT_H
#define THERMOSTAT_H

#include <string>

#include "core/System.h"
class Thermostat {
public:
  enum class ThermostatType { BERENDSEN, LANGEVIN };

private:
  double tau_t{0.0};
  double pref_temperature{0.0};

public:
  virtual ~Thermostat() = default;
  virtual void applyTemperatureControl(System &sys) = 0;
  virtual std::string getData() const = 0;

  inline virtual const ThermostatType getThermostatType() const = 0;
};

#endif
