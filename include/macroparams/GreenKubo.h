#ifndef GREEN_KUBO_H
#define GREEN_KUBO_H

#include "classes/Matrix3.h"
#include "classes/Particle.h"
#include "classes/Vector3.h"
#include "core/System.h"

/**
 * @brief Класс для расчета коэффициентов теплопроводности и вязкости по методу Грина-Кубо.
 * @details Класс предоставляет методы для расчета коэффициентов теплопроводности и вязкости.
 */
class GreenKubo {
private:
  /* data */
public:
  GreenKubo() = default;
  ~GreenKubo() = default;

  /**
   * @brief Расчет автокорреляционной функции скорости.
   * @param sys - система частиц.
   * @param init_vel - начальные скорости частиц.
   * @return Автокорреляционная функция скорости.
   */
  inline const double getACFV(System &sys,
                              std::vector<Vector3<double>> &init_vel) const {
    double acfv_x = 0, acfv_y = 0, acfv_z = 0;
    int pn = sys.particleNumber();
    std::vector<Particle> &particles = sys.particles();
    for (int i = 0; i < pn; i++) {
      acfv_x += particles[i].getVelocityX() * init_vel[i].x();
      acfv_y += particles[i].getVelocityY() * init_vel[i].y();
      acfv_z += particles[i].getVelocityZ() * init_vel[i].z();
    }
    return (acfv_x + acfv_y + acfv_z);
  }

  /**
   * @brief Расчет автокорреляционной функции давления.
   * @param sys - система частиц.
   * @param init_press_tensor - начальный тензор давления.
   * @return Автокорреляционная функция давления.
   */
  inline const double getACFP(System &sys, Matrix3 &init_press_tensor) const {
    double acfp_xy = 0, acfp_yz = 0, acfp_zx = 0;
    const Matrix3 &P = sys.pressureTensors();
    acfp_xy = init_press_tensor.xy() * P.xy();
    acfp_yz = init_press_tensor.yz() * P.yz();
    acfp_zx = init_press_tensor.zx() * P.zx();

    return (acfp_xy + acfp_yz + acfp_zx);
  }

  
};

#endif
