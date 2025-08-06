#ifndef GREEN_KUBO_H
#define GREEN_KUBO_H

#include "classes/Particles.h"
#include "core/System.h"

class GreenKubo {
private:
  /* data */
public:
  inline const double getACFV(System &sys, std::vector<double> &init_vel_x,
                              std::vector<double> &init_vel_y,
                              std::vector<double> &init_vel_z) const {
    double acfv_x = 0, acfv_y = 0, acfv_z = 0;
    int pn = sys.particleNumber();
    Particles &particles = sys.particles();
    for (int i = 0; i < pn; i++) {
      acfv_x += particles.velocity_x_[i] * init_vel_x[i];
      acfv_y += particles.velocity_y_[i] * init_vel_y[i];
      acfv_z += particles.velocity_z_[i] * init_vel_z[i];
    }
    return (acfv_x + acfv_y + acfv_z);
  }
  inline const double getACFP(System &sys, double &init_press_tensor_xy,
                              double &init_press_tensor_yz,
                              double &init_press_tensor_zx) const {
    double acfp_xy = 0, acfp_yz = 0, acfp_zx = 0;
    acfp_xy = init_press_tensor_xy * sys.pressureXY();
    acfp_yz = init_press_tensor_yz * sys.pressureYZ();
    acfp_zx = init_press_tensor_zx * sys.pressureZX();

    return (acfp_xy + acfp_yz + acfp_zx);
  }

  GreenKubo(/* args */) = default;
  ~GreenKubo() = default;
};

#endif
