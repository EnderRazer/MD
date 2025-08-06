#ifndef POTENTIAL_EAM
#define POTENTIAL_EAM

#include <ostream>
#include <sstream>
#include <vector>

#include "nlohmann/json.hpp"

#include "core/Settings.h"
#include "core/System.h"

#include "Potential.h"

using json = nlohmann::json;
// EAM
class EAM : public Potential {
public:
  const PotentialType type_{PotentialType::EAM};

  const double r_e_;
  const double inv_r_e_;
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
  const int m_;
  const int n_;
  const double energy_unit_;
  const double r_cut_;
  const double r_cut_sqr_;

  const double precision_;

  std::vector<double> t_rho_f_;
  std::vector<double> t_d_rho_f_;
  std::vector<double> t_mu_;
  std::vector<double> t_d_mu_;

  std::vector<double> t_om1_;
  std::vector<double> t_d_om1_;

  std::vector<double> t_om2_;
  std::vector<double> t_d_om2_;

  std::vector<double> t_om3_;
  std::vector<double> t_d_om3_;

  inline double rho_f(double r) const {
    double r_over_R_e = r * inv_r_e_;

    double exp_term = exp(-beta_ * (r_over_R_e - 1));
    double pow_term = pow(r_over_R_e - lambda_, n_);

    return (f_e_ * exp_term) / (1 + pow_term);
  }
  inline double int_rho_f(double r) const {
    double idx = r / precision_;
    int i = static_cast<int>(idx);

    if (i >= t_rho_f_.size() - 1)
      return 0;
    double frac = idx - i;
    return t_rho_f_[i] + frac * (t_rho_f_[i + 1] - t_rho_f_[i]);
  }

  inline double d_rho_f(double r) const {
    double r_over_Re = r * inv_r_e_;

    double exp_term = exp(-beta_ * (r_over_Re - 1));
    double pow_term = pow(r_over_Re - lambda_, n_);

    double var1 = (-beta_ * inv_r_e_) * exp_term * (1 + pow_term);
    double var2 = (n_ * inv_r_e_) * pow(r_over_Re - lambda_, n_ - 1) * exp_term;
    double var3 = 1 + pow_term;

    return f_e_ * (var1 - var2) / (var3 * var3);
  }
  inline double int_d_rho_f(double r) const {
    double idx = r / precision_;
    int i = static_cast<int>(idx);

    if (i >= t_d_rho_f_.size() - 1)
      return 0;
    double frac = idx - i;
    return t_d_rho_f_[i] + frac * (t_d_rho_f_[i + 1] - t_d_rho_f_[i]);
  }

  inline double mu(double r) const {
    double r_over_Re = r * inv_r_e_;

    double exp_term1 = exp(-alpha_ * (r_over_Re - 1));
    double pow_term1 = pow(r_over_Re - k_, m_);

    double exp_term2 = exp(-beta_ * (r_over_Re - 1));
    double pow_term2 = pow(r_over_Re - lambda_, n_);

    return (a_ * exp_term1 / (1 + pow_term1)) -
           (b_ * exp_term2 / (1 + pow_term2));
  }
  inline double int_mu(double r) const {
    double idx = r / precision_;
    int i = static_cast<int>(idx);

    if (i >= t_mu_.size() - 1)
      return 0;
    double frac = idx - i;
    return t_mu_[i] + frac * (t_mu_[i + 1] - t_mu_[i]);
  }

  inline double d_mu(double r) const {
    double r_over_Re = r * inv_r_e_;

    double exp_term1 = exp(-alpha_ * (r_over_Re - 1));
    double pow_term1 = pow(r_over_Re - k_, m_);

    double var11 = (-a_ * alpha_ * inv_r_e_) * exp_term1 * (1 + pow_term1);
    double var12 =
        (m_ * inv_r_e_) * pow(r_over_Re - k_, m_ - 1) * a_ * exp_term1;
    double var1 = (var11 - var12) / ((1 + pow_term1) * (1 + pow_term1));

    double exp_term2 = exp(-beta_ * (r_over_Re - 1));
    double pow_term2 = pow(r_over_Re - lambda_, n_);

    double var21 =
        (n_ * inv_r_e_) * pow(r_over_Re - lambda_, n_ - 1) * b_ * exp_term2;
    double var22 = (-b_ * beta_ * inv_r_e_) * exp_term2 * (1 + pow_term2);
    double var2 = (var21 - var22) / ((1 + pow_term2) * (1 + pow_term2));

    return var1 + var2;
  }
  inline double int_d_mu(double r) const {
    double idx = r / precision_;
    int i = static_cast<int>(idx);

    if (i >= t_d_mu_.size() - 1)
      return 0;
    double frac = idx - i;
    return t_d_mu_[i] + frac * (t_d_mu_[i + 1] - t_d_mu_[i]);
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
      : r_e_(params_json.value("r_e", 0.0)), inv_r_e_(1 / r_e_),
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
        eta_(params_json.value("eta", 0.0)), m_(params_json.value("m", 0)),
        n_(params_json.value("n", 0)),
        energy_unit_(params_json.value("energy_unit", 0.0)),
        r_cut_(params_json.value("r_cut", 0.0)), r_cut_sqr_(r_cut_ * r_cut_),
        precision_(params_json.value("interpolation_dr", 0.0)) {
    int size = r_cut_ / precision_;
    std::cout << "Заполнение табличных значений"
              << "\n\tТочность: " << precision_
              << "\n\tКоличество записей: " << size << std::endl;

    t_rho_f_.reserve(size);
    t_d_rho_f_.reserve(size);

    t_mu_.reserve(size);
    t_d_mu_.reserve(size);

    t_om1_.reserve(size);
    t_d_om1_.reserve(size);

    t_om2_.reserve(size);
    t_d_om2_.reserve(size);

    t_om3_.reserve(size);
    t_d_om3_.reserve(size);

    for (int i = 0; i < size; i++) {
      double r = i * precision_;

      t_rho_f_[i] = rho_f(r);
      t_d_rho_f_[i] = d_rho_f(r);

      t_mu_[i] = mu(r);
      t_d_mu_[i] = d_mu(r);
    }
  };

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

  inline double getFU_fast_branchless(double rho_f_i, double rho_f_j,
                                      double r_ij) const {
    const double r_over_Re = r_ij * inv_r_e_;
    const double t1 = r_over_Re - 1.0;
    const double t_k = r_over_Re - k_;
    const double t_l = r_over_Re - lambda_;

    const double exp_alpha = exp(-alpha_ * t1);
    const double exp_beta = exp(-beta_ * t1);

    const double pow_k_m = pow(t_k, m_);
    const double pow_k_m1 = pow(t_k, m_ - 1);
    const double pow_l_n = pow(t_l, n_);
    const double pow_l_n1 = pow(t_l, n_ - 1);

    // ---------------------------------------------
    // d_om branchless
    // ---------------------------------------------
    auto dom_branchless = [&](double rho) {
      const double v1 = d_om1(rho);
      const double v2 = d_om2(rho);
      const double v3 = d_om3(rho);
      const double m1 = (rho < rho_n_) ? 1.0 : 0.0;
      const double m3 = (rho >= rho_0_) ? 1.0 : 0.0;
      return m1 * v1 + (1.0 - m1 - m3) * v2 + m3 * v3;
    };

    const double dom_sum = dom_branchless(rho_f_i) + dom_branchless(rho_f_j);

    // ---------------------------------------------
    // d_rho_f(r)
    // ---------------------------------------------
    const double one_plus_pow_l_n = 1.0 + pow_l_n;
    const double denom_rho2 = one_plus_pow_l_n * one_plus_pow_l_n;

    const double var1_rho = (-beta_ * inv_r_e_) * exp_beta * one_plus_pow_l_n;
    const double var2_rho = (n_ * inv_r_e_) * pow_l_n1 * exp_beta;

    const double d_rho_f_val = f_e_ * (var1_rho - var2_rho) / denom_rho2;

    // ---------------------------------------------
    // d_mu(r)
    // ---------------------------------------------
    const double one_plus_pow_k_m = 1.0 + pow_k_m;
    const double denom_mu1 = one_plus_pow_k_m * one_plus_pow_k_m;

    const double var11 =
        (-a_ * alpha_ * inv_r_e_) * exp_alpha * one_plus_pow_k_m;
    const double var12 = (m_ * inv_r_e_) * pow_k_m1 * a_ * exp_alpha;
    const double var1_mu = (var11 - var12) / denom_mu1;

    const double var21 = (n_ * inv_r_e_) * pow_l_n1 * b_ * exp_beta;
    const double var22 = (-b_ * beta_ * inv_r_e_) * exp_beta * one_plus_pow_l_n;
    const double var2_mu = (var21 - var22) / denom_rho2;

    const double d_mu_val = var1_mu + var2_mu;

    return -energy_unit_ * (dom_sum * d_rho_f_val + d_mu_val);
  }

  inline double getPairPart(double r) const override { return mu(r); }
  inline double getDensityPart(double r) const override { return rho_f(r); }
  inline double getDerPairPart(double r) const override { return d_mu(r); }
  inline double getDerDensityPart(double r) const override {
    return d_rho_f(r);
  }
  inline double getCloud(double rho) const override { return om(rho); }
  inline double getDerCloud(double rho) const override { return d_om(rho); }

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
