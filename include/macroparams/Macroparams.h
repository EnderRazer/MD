#ifndef MACROPARAMS_H
#define MACROPARAMS_H

#include "nlohmann/json.hpp"

#include "core/Settings.h"
#include "core/System.h"

using json = nlohmann::json;

class Macroparams {
private:
  Settings &settings_;
  const double T_CONST;
  bool toggle_{false};

public:
  inline const bool enabled() const { return toggle_; }
  inline double getTemperature(System &sys) const {
    return sys.eTerm() * T_CONST;
  }

  inline double getPressure(System &sys) const {
    double sumMV_xx, sumMV_xy, sumMV_xz;
    double sumMV_yx, sumMV_yy, sumMV_yz;
    double sumMV_zx, sumMV_zy, sumMV_zz;
    
    double sumVirials_xx, sumVirials_xy, sumVirials_xz;
    double sumVirials_yx, sumVirials_yy, sumVirials_yz;
    double sumVirials_zx, sumVirials_zy, sumVirials_zz;

    Particles &particles = sys.particles();

    for (int i = 0; i < particles.size(); i++) {
      sumMV_xx += particles.mass_[i] * (particles.velocity_x_[i] - sys.vcmX()) *
                    (particles.velocity_x_[i] - sys.vcmX()); // xx
      sumMV_xy += particles.mass_[i] * (particles.velocity_x_[i] - sys.vcmX()) *
                    (particles.velocity_y_[i] - sys.vcmY()); // xy
      sumMV_xz += particles.mass_[i] * (particles.velocity_x_[i] - sys.vcmX()) *
                    (particles.velocity_z_[i] - sys.vcmZ()); // xz

      sumMV_yx += particles.mass_[i] * (particles.velocity_y_[i] - sys.vcmY()) *
                    (particles.velocity_x_[i] - sys.vcmX()); // yx
      sumMV_yy += particles.mass_[i] * (particles.velocity_y_[i] - sys.vcmY()) *
                    (particles.velocity_y_[i] - sys.vcmY()); // yy
      sumMV_yz += particles.mass_[i] * (particles.velocity_y_[i] - sys.vcmY()) *
                    (particles.velocity_z_[i] - sys.vcmZ()); // yz

      sumMV_zx += particles.mass_[i] * (particles.velocity_z_[i] - sys.vcmZ()) *
                    (particles.velocity_x_[i] - sys.vcmX()); // zx
      sumMV_zy += particles.mass_[i] * (particles.velocity_z_[i] - sys.vcmZ()) *
                    (particles.velocity_y_[i] - sys.vcmY()); // zy
      sumMV_zz += particles.mass_[i] * (particles.velocity_x_[i] - sys.vcmX()) *
                    (particles.velocity_z_[i] - sys.vcmZ()); // zz

      sumVirials_xx += particles.virial_xx_[i];
      sumVirials_xy += particles.virial_xy_[i];
      sumVirials_xz += particles.virial_xz_[i];

      sumVirials_yx += particles.virial_yx_[i];
      sumVirials_yy += particles.virial_yy_[i];
      sumVirials_yz += particles.virial_yz_[i];

      sumVirials_zx += particles.virial_zx_[i];
      sumVirials_zy += particles.virial_zy_[i];
      sumVirials_zz += particles.virial_zz_[i];
    }
    double volume = sys.dimensions().volume();

    sys.pressureXX() = (sumMV_xx + 0.5 * sumVirials_xx) / volume;
    sys.pressureXY() = (sumMV_xy + 0.5 * sumVirials_xy) / volume;
    sys.pressureXZ() = (sumMV_xz + 0.5 * sumVirials_xz) / volume;

    sys.pressureYX() = (sumMV_yx + 0.5 * sumVirials_yx) / volume;
    sys.pressureYY() = (sumMV_yy + 0.5 * sumVirials_yy) / volume;
    sys.pressureYZ() = (sumMV_yz + 0.5 * sumVirials_yz) / volume;

    sys.pressureZX() = (sumMV_zx + 0.5 * sumVirials_zx) / volume;
    sys.pressureZY() = (sumMV_zy + 0.5 * sumVirials_zy) / volume;
    sys.pressureZZ() = (sumMV_zz + 0.5 * sumVirials_zz) / volume;
    
    // Расчет давления по XX,YY,ZZ компонентам
    return (sys.pressureXX() + sys.pressureYY() + sys.pressureZZ()) / 3;
  }

  Macroparams(json &config, Settings &settings)
      : settings_(settings), toggle_(config.value("toggle", false)),
        T_CONST((2 / (settings_.constants().D *
                      settings_.constants().KBOLTZMAN))) {};

  ~Macroparams() = default;
};
#endif
