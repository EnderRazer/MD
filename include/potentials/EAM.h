#ifndef POTENTIAL_EAM
#define POTENTIAL_EAM

#include <sstream>

#include "nlohmann/json.hpp"

#include "core/Settings.h"
#include "core/System.h"

#include "Potential.h"

using json = nlohmann::json;
// EAM
class EAM : public Potential {
private:
  const PotentialType type_{PotentialType::EAM};

  const double r_e_;
  const double f_e_;
  const double rho_e_;
  const double rho_s_;
  const double rho_n_;
  const double rho_0_;
  const double om_e_;
  const std::vector<double> om_n_;
  const std::vector<double> om_;
  const double alpha_;
  const double beta_;
  const double a_;
  const double b_;
  const double k_;
  const double lambda_;
  const double eta_;
  const double m_;
  const double n_;
  const double energy_unit_;
  const double r_cut_;
  const double r_cut_sqr_;

  inline Mu_RhoF mu_rho_f(double r) const {
    double r_over_R_e = r / r_e_;
    double exp_term1 = exp(-alpha_ * (r_over_R_e - 1));
    double pow_term1 = pow(r_over_R_e - k_, m_);
    double exp_term2 = exp(-beta_ * (r_over_R_e - 1));
    double pow_term2 = pow(r_over_R_e - lambda_, n_);

    return {(f_e_ * exp_term2) / (1 + pow_term2),(a_ * exp_term1 / (1 + pow_term1)) -
           (b_ * exp_term2 / (1 + pow_term2))};
  }
  inline double rho_f(double r) const {
    double r_over_R_e = r / r_e_;
    double exp_term = exp(-beta_ * (r_over_R_e - 1));
    double pow_term = pow(r_over_R_e - lambda_, n_);
    return (f_e_ * exp_term) / (1 + pow_term);
  }
  inline double d_rho_f(double r) const {
    double r_over_Re = r / r_e_;
    double exp_term = exp(-beta_ * (r_over_Re - 1));
    double pow_term = pow(r_over_Re - lambda_, n_);
    double var1 = (-beta_ / r_e_) * exp_term * (1 + pow_term);
    double var2 = (n_ / r_e_) * pow(r_over_Re - lambda_, n_ - 1) * exp_term;
    double var3 = 1 + pow_term;
    return f_e_ * (var1 - var2) / (var3 * var3);
  }

  inline double mu(double r) const {
    double r_over_Re = r / r_e_;
    double exp_term1 = exp(-alpha_ * (r_over_Re - 1));
    double pow_term1 = pow(r_over_Re - k_, m_);
    double exp_term2 = exp(-beta_ * (r_over_Re - 1));
    double pow_term2 = pow(r_over_Re - lambda_, n_);

    return (a_ * exp_term1 / (1 + pow_term1)) -
           (b_ * exp_term2 / (1 + pow_term2));
  }
  inline double d_mu(double r) const {
    double r_over_Re = r / r_e_;
    double exp_term1 = exp(-alpha_ * (r_over_Re - 1));
    double pow_term1 = pow(r_over_Re - k_, m_);
    double var11 = (-a_ * alpha_ / r_e_) * exp_term1 * (1 + pow_term1);
    double var12 = (m_ / r_e_) * pow(r_over_Re - k_, m_ - 1) * a_ * exp_term1;
    double var1 = (var11 - var12) / ((1 + pow_term1) * (1 + pow_term1));

    double exp_term2 = exp(-beta_ * (r_over_Re - 1));
    double pow_term2 = pow(r_over_Re - lambda_, n_);
    double var21 =
        (n_ / r_e_) * pow(r_over_Re - lambda_, n_ - 1) * b_ * exp_term2;
    double var22 = (-b_ * beta_ / r_e_) * exp_term2 * (1 + pow_term2);
    double var2 = (var21 - var22) / ((1 + pow_term2) * (1 + pow_term2));

    return var1 + var2;
  }

  inline double om1(double rho) const {
    double om_sum = 0;
    double rho_over_RHO_n = rho / rho_n_ - 1;
    double rho_pow = 1;

    for (int i = 0; i < 4; i++) {
      om_sum += om_n_[i] * rho_pow;
      rho_pow *= rho_over_RHO_n;
    }

    return om_sum;
  }
  inline double d_om1(double rho) const {
    double d_om_sum = 0;
    double rho_over_RHO_n = rho / rho_n_ - 1;
    double rho_pow = 1;

    for (int i = 1; i < 4; i++) {
      d_om_sum += i * om_n_[i] * rho_pow / rho_n_;
      rho_pow *= rho_over_RHO_n;
    }

    return d_om_sum;
  }

  inline double om2(double rho) const {
    double om_sum = 0;
    double rho_over_RHO_e = rho / rho_e_ - 1;
    double rho_pow = 1;

    for (int i = 0; i < 4; i++) {
      om_sum += om_[i] * rho_pow;
      rho_pow *= rho_over_RHO_e;
    }

    return om_sum;
  }
  inline double d_om2(double rho) const {
    double d_om_sum = 0;
    double rho_over_RHO_e = rho / rho_e_ - 1;
    double rho_pow = 1;

    for (int i = 1; i < 4; i++) {
      d_om_sum += i * om_[i] * rho_pow / rho_e_;
      rho_pow *= rho_over_RHO_e;
    }

    return d_om_sum;
  }

  inline double om3(double rho) const {
    double rho_over_RHO_s = rho / rho_s_;
    double pow_term = pow(rho_over_RHO_s, eta_);
    double log_term = log(pow_term);

    return om_e_ * (1 - log_term) * pow_term;
  }
  inline double d_om3(double rho) const {
    double rho_over_RHO_s = rho / rho_s_;
    double pow_term = pow(rho_over_RHO_s, eta_);
    double log_term = log(pow_term);

    return -om_e_ * eta_ * pow(rho_over_RHO_s, eta_ - 1) * log_term / rho_s_;
  }

  inline double om(double rho) const {
    if (rho < rho_n_) {
      return om1(rho);
    } else if (rho >= rho_0_) {
      return om3(rho);
    } else
      return om2(rho);
  }
  inline double d_om(double rho) const {
    if (rho < rho_n_) {
      return d_om1(rho);
    }
    if (rho >= rho_0_) {
      return d_om3(rho);
    }
    return d_om2(rho);
  }

public:
  ~EAM() = default;

  EAM(json &params_json)
      : r_e_(params_json.value("r_e", 0.0)),
        f_e_(params_json.value("f_e", 0.0)),
        rho_e_(params_json.value("rho_e", 0.0)),
        rho_s_(params_json.value("rho_s", 0.0)), rho_n_(rho_e_ * 0.85),
        rho_0_(rho_e_ * 1.15), om_e_(params_json.value("om_e", 0.0)),
        om_n_(params_json["om_n"].get<std::vector<double>>()),
        om_(params_json["om"].get<std::vector<double>>()),
        alpha_(params_json.value("alpha", 0.0)),
        beta_(params_json.value("beta", 0.0)), a_(params_json.value("a", 0.0)),
        b_(params_json.value("b", 0.0)), k_(params_json.value("k", 0.0)),
        lambda_(params_json.value("lambda", 0.0)),
        eta_(params_json.value("eta", 0.0)), m_(params_json.value("m", 0.0)),
        n_(params_json.value("n", 0.0)),
        energy_unit_(params_json.value("energy_unit", 0.0)),
        r_cut_(params_json.value("r_cut", 0.0)), r_cut_sqr_(r_cut_ * r_cut_) {};

  inline PotentialType getPotentialType() const override {
    return type_;
  } // Getter for type

  std::string getData() const override {
    std::ostringstream oss;
    oss.precision(16);
    oss << "Using EAM potential!"
        << "\n\tR_e = " << r_e_ << "\n\tF_e = " << f_e_
        << "\n\tRho_e = " << rho_e_ << "\n\tRho_s = " << rho_s_
        << "\n\tRho_n = " << rho_n_ << "\n\tRho_0 = " << rho_0_
        << "\n\tOM_e = " << om_e_ << "\n\tOM_n = [" << om_n_[0] << ", "
        << om_n_[1] << ", " << om_n_[2] << ", " << om_n_[3] << "]"
        << "\n\tOM = [" << om_[0] << ", " << om_[1] << ", " << om_[2] << ", "
        << om_[3] << "]"
        << "\n\tALPHA = " << alpha_ << "\n\tBETA = " << beta_
        << "\n\tA = " << a_ << "\n\tB = " << b_ << "\n\tK = " << k_
        << "\n\tLAMBDA = " << lambda_ << "\n\tETA = " << eta_
        << "\n\tM = " << m_ << "\n\tN = " << n_
        << "\n\tENERGY_UNIT = " << energy_unit_ << "\n\tRCUT = " << r_cut_;

    return oss.str();
  }

  inline double getU(double rho_f, double mu) const override {
    return energy_unit_ * (om(rho_f) + 0.5 * mu);
  };

  inline double getFU(double rho_f_i, double rho_f_j,
                      double r_ij) const override {
    return -energy_unit_ *
           ((d_om(rho_f_i) + d_om(rho_f_j)) * d_rho_f(r_ij) + d_mu(r_ij));
  };

  inline double getPairPart(double r) const override { return mu(r); }
  inline double getDensityPart(double r) const override { return rho_f(r); }
  inline Mu_RhoF getPairDesityPart(double r) const override { return mu_rho_f(r); }
  inline double getDerPairPart(double r) const override { return d_mu(r); }
  inline double getDerDensityPart(double r) const override {return d_rho_f(r); }
  inline double getCloud(double rho) const override { return om(rho); }
  inline double getDerCloud(double rho) const override { return d_om(rho); }

  // LJ overrides, not implemented for EAM
  inline double getU(double r) const override {
    throw std::runtime_error("getU(r) not implemented for EAM");
  }; // Потенциальная энергия
  inline double getFU(double r) const override {
    throw std::runtime_error("getFU(r) not implemented for EAM");
  }; // Сила потенциала
  inline const PotentialResult getAll(double r) const override {
    throw std::runtime_error("getAll(r) not implemented for EAM");
  }; // Сила и энергия потенциала

  inline double getRcut() const override { return r_cut_; }
  inline double getSqrRcut() const override { return r_cut_sqr_; }
};
#endif
