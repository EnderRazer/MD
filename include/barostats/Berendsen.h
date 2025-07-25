#ifndef BSBERENDSEN_H
#define BSBERENDSEN_H

#include "nlohmann/json.hpp"
#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>

#include "classes/Particle.h"
#include "classes/Vector3.h"
#include "core/System.h"

#include "Barostat.h"

using json = nlohmann::json;

class BarostatBerendsen : public Barostat {
private:
  const BarostatType type_{BarostatType::BERENDSEN};
  bool toggle_{false};
  double pref_pressure_{0.0};
  double tau_b_{0.0};

  double dt_over_tau_{0.0};

public:
  BarostatBerendsen(json &config, Settings &settings) {
    toggle_ = config.value("toggle", false);
    pref_pressure_ = config.value("pref_pressure", 0.0);
    tau_b_ = config.value("tau", 0.0);

    dt_over_tau_ = settings.dt() / tau_b_;
  };
  void applyPressureControl(System &sys) override {
    double hi = 1 - dt_over_tau_ * (pref_pressure_ - sys.pressure());
    hi = std::clamp(hi,0.5,1.5);
    double mu = std::cbrt(hi);

    for (Particle &particle : sys.particles()) {
      Vector3<double> new_coord = particle.coord() * mu;
      particle.setCoord(new_coord);
    }
    Vector3<double> new_sizes = sys.dimensions().sizes() * mu;
    sys.dimensions().setSizes(new_sizes);
  }
  std::string getData() const override {
    std::ostringstream oss;
    oss.precision(16);
    oss << "Using Berendsen barostat!"
        << "\n\tToggled: " << toggle_
        << "\n\tPrefered pressure: " << pref_pressure_
        << "\n\tTau_b: " << tau_b_ << "\n";
    return oss.str();
  }

  inline const BarostatType getBarostatType() const override { return type_; }

  // Запрещаем копирование
  BarostatBerendsen(const BarostatBerendsen &) = delete;
  BarostatBerendsen &operator=(const BarostatBerendsen &) = delete;
};
#endif
