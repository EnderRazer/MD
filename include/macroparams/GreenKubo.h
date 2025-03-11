#ifndef GREEN_KUBO_H
#define GREEN_KUBO_H

class GreenKubo {
private:
  /* data */
public:
  inline const double getACFV(System &sys, Ensemble &ens) const {
    double acfv_x = 0, acfv_y = 0, acfv_z = 0;
    std::vector<Vector3<double>> &V0 = ens.init_velocity_;
    int pn = sys.particleNumber();
    std::vector<Particle> &particles = sys.particles();
    for (int i = 0; i < pn; i++) {
      acfv_x += particles[i].getVelocityX() * V0[i].x();
      acfv_y += particles[i].getVelocityY() * V0[i].y();
      acfv_z += particles[i].getVelocityZ() * V0[i].z();
    }
    return (acfv_x + acfv_y + acfv_z) / (3 * sys.particleNumber());
  }
  inline const double getACFP(System &sys, Ensemble &ens) const {
    double acfp_xy = 0, acfp_yz = 0, acfp_zx = 0;
    const Matrix3 &P0 = ens.init_pressure_tensors_;
    const Matrix3 &P = sys.pressureTensors();
    acfp_xy = P0.xy() * P.xy();
    acfp_yz = P0.yz() * P.yz();
    acfp_zx = P0.zx() * P.zx();

    return (acfp_xy + acfp_yz + acfp_zx) * sys.dimensions().volume() /
           (3 * sys.temperature());
  }
  GreenKubo(/* args */);
  ~GreenKubo();
};

GreenKubo::GreenKubo(/* args */) {}

GreenKubo::~GreenKubo() {}

#endif