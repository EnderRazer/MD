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
      sumMV_xx += particles.mass(i) * (particles.velocityX(i) - sys.vcmX()) *
                    (particles.velocityX(i) - sys.vcmX()); // xx
      sumMV_xy += particles.mass(i) * (particles.velocityX(i) - sys.vcmX()) *
                    (particles.velocityY(i) - sys.vcmY()); // xy
      sumMV_xz += particles.mass(i) * (particles.velocityX(i) - sys.vcmX()) *
                    (particles.velocityZ(i) - sys.vcmZ()); // xz

      sumMV_yx += particles.mass(i) * (particles.velocityY(i) - sys.vcmY()) *
                    (particles.velocityX(i) - sys.vcmX()); // yx
      sumMV_yy += particles.mass(i) * (particles.velocityY(i) - sys.vcmY()) *
                    (particles.velocityY(i) - sys.vcmY()); // yy
      sumMV_yz += particles.mass(i) * (particles.velocityY(i) - sys.vcmY()) *
                    (particles.velocityZ(i) - sys.vcmZ()); // yz

      sumMV_zx += particles.mass(i) * (particles.velocityZ(i) - sys.vcmZ()) *
                    (particles.velocityX(i) - sys.vcmX()); // zx
      sumMV_zy += particles.mass(i) * (particles.velocityZ(i) - sys.vcmZ()) *
                    (particles.velocityY(i) - sys.vcmY()); // zy
      sumMV_zz += particles.mass(i) * (particles.velocityX(i) - sys.vcmX()) *
                    (particles.velocityZ(i) - sys.vcmZ()); // zz

      sumVirials_xx += particles.virialXX(i);
      sumVirials_xy += particles.virialXY(i);
      sumVirials_xz += particles.virialXZ(i);

      sumVirials_yx += particles.virialYX(i);
      sumVirials_yy += particles.virialYY(i);
      sumVirials_yz += particles.virialYZ(i);

      sumVirials_zx += particles.virialZX(i);
      sumVirials_zy += particles.virialZY(i);
      sumVirials_zz += particles.virialZZ(i);
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
