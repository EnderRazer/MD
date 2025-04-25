#ifndef MACROPARAMS_H
#define MACROPARAMS_H

#include "nlohmann/json.hpp"

#include "classes/Matrix3.h"
#include "classes/Vector3.h"
#include "core/Settings.h"
#include "core/System.h"
#include <iostream>

using json = nlohmann::json;

class Macroparams {
private:
  Settings &settings_;
  const double T_CONST;
  bool toggle_{false};

public:
  inline const bool enabled() const { return toggle_; }
  inline double getTemperature(System &sys) const {
    return sys.energies().get(Energy::EnergyType::Thermodynamic) * T_CONST;
  }

  inline double getPressure(System &sys) const {
    Matrix3 sumMV{};
    Matrix3 sumVirials{};
    for (Particle &p : sys.particles()) {
      sumMV.xx() += p.getMass() * (p.velocity().x() - sys.vcm().x()) *
                    (p.velocity().x() - sys.vcm().x()); // xx
      sumMV.xy() += p.getMass() * (p.velocity().x() - sys.vcm().x()) *
                    (p.velocity().y() - sys.vcm().y()); // xy
      sumMV.xz() += p.getMass() * (p.velocity().x() - sys.vcm().x()) *
                    (p.velocity().z() - sys.vcm().z()); // xz

      sumMV.yx() += p.getMass() * (p.velocity().y() - sys.vcm().y()) *
                    (p.velocity().x() - sys.vcm().x()); // yx
      sumMV.yy() += p.getMass() * (p.velocity().y() - sys.vcm().y()) *
                    (p.velocity().y() - sys.vcm().y()); // yy
      sumMV.yz() += p.getMass() * (p.velocity().y() - sys.vcm().y()) *
                    (p.velocity().z() - sys.vcm().z()); // yz

      sumMV.zx() += p.getMass() * (p.velocity().z() - sys.vcm().z()) *
                    (p.velocity().x() - sys.vcm().x()); // zx
      sumMV.zy() += p.getMass() * (p.velocity().z() - sys.vcm().z()) *
                    (p.velocity().y() - sys.vcm().y()); // zy
      sumMV.zz() += p.getMass() * (p.velocity().z() - sys.vcm().z()) *
                    (p.velocity().z() - sys.vcm().z()); // zz

      sumVirials.xx() += p.virials().xx();
      sumVirials.xy() += p.virials().xy();
      sumVirials.xz() += p.virials().xz();

      sumVirials.yx() += p.virials().yx();
      sumVirials.yy() += p.virials().yy();
      sumVirials.yz() += p.virials().yz();

      sumVirials.zx() += p.virials().zx();
      sumVirials.zy() += p.virials().zy();
      sumVirials.zz() += p.virials().zz();
    }
    Matrix3 p_tensors = (sumMV + 0.5 * sumVirials) / sys.dimensions().volume();
    sys.setPressureTensors(p_tensors);
    // Расчет давления по XX,YY,ZZ компонентам
    return (p_tensors.xx() + p_tensors.yy() + p_tensors.zz()) / 3;
  }

  Macroparams(json &config, Settings &settings)
      : settings_(settings), toggle_(config.value("toggle", false)),
        T_CONST((2 / (settings_.constants().D *
                      settings_.constants().KBOLTZMAN))) {};

  ~Macroparams() = default;
};
#endif
