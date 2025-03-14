#ifndef POTENTIAL_LJ
#define POTENTIAL_LJ
using json = nlohmann::json;
// LJ
class LJ : public Potential
{
private:
  const PotentialType type_{PotentialType::LJ};
  const double epsilon_;
  const double sigma_;

  const double sigma_p2_;
  const double sigma_p6_;

  const double r_cut_;
  const double r_cut_sqr_;

  const double epsilon_t4_;
  const double epsilon_t24_;

  const double eps_t4_sigma_p6_;
  const double eps_t4_sigma_p12_;

  double r_cut_pot_;

  const bool precompute_{false};

public:
  ~LJ() = default;

  LJ(json &params_json)
      : epsilon_(params_json.value("epsilon", 0.0)),
        sigma_(params_json.value("sigma", 0.0)), sigma_p2_(sigma_ * sigma_),
        sigma_p6_(sigma_p2_ * sigma_p2_ * sigma_p2_), r_cut_(2.5 * sigma_),
        r_cut_sqr_(r_cut_ * r_cut_), epsilon_t4_(4 * epsilon_),
        epsilon_t24_(24 * epsilon_), eps_t4_sigma_p6_(epsilon_t4_ * sigma_p6_),
        eps_t4_sigma_p12_(epsilon_t4_ * sigma_p6_ * sigma_p6_)
  {
    double r_p6 = r_cut_sqr_ * r_cut_sqr_ * r_cut_sqr_;
    double r_p12 = r_p6 * r_p6;
    r_cut_pot_ = eps_t4_sigma_p12_ / r_p12 - eps_t4_sigma_p6_ / r_p6;
  };

  inline PotentialType getPotentialType() const
  {
    return type_;
  } // Getter for type

  std::string getData() const override
  {
    std::ostringstream oss;
    oss.precision(16);
    oss << "Using LJ potential!\n\tEpsilon: " << epsilon_
        << "\n\tSigma: " << sigma_ << "\n\t" << "Rcut: " << r_cut_
        << "\n\tRcut_pot: " << r_cut_pot_ << "\n";
    return oss.str();
  }

  inline double getU(double r_p2) const override
  {
    double r_p6 = r_p2 * r_p2 * r_p2;
    return eps_t4_sigma_p12_ / (r_p6 * r_p6) - eps_t4_sigma_p6_ / r_p6 -
           r_cut_pot_;
  };
  inline double getFU(double r_p2) const override
  {
    double r_p6 = r_p2 * r_p2 * r_p2;
    return 6 *
           (2 * eps_t4_sigma_p12_ / (r_p6 * r_p6) - eps_t4_sigma_p6_ / r_p6);
  }

  inline const LJResult getAll(double r_p2) const override
  {
    double r_p6 = r_p2 * r_p2 * r_p2;
    double eps_div_r_p12 = eps_t4_sigma_p12_ / (r_p6 * r_p6);
    double eps_div_r_p6 = eps_t4_sigma_p6_ / r_p6;

    double u = eps_div_r_p12 - eps_div_r_p6 - r_cut_pot_;
    double fu = 6 * (2 * eps_div_r_p12 - eps_div_r_p6);
    return {u, fu};
  }

  inline double getU(double rho_f, double mu) const override
  {
    throw std::runtime_error("getU(rho_f, mu) not implemented for LJ");
  }; // Потенциальная энергия
  inline double getFU(double rho_f_i, double rho_f_j, double d_rho_f_ij,
                      double d_mu_ij) const override
  {
    throw std::runtime_error(
        "getFU(rho_f_i, rho_f_j, d_rho_f_ij, d_mu_ij) not implemented for LJ");
  }; // Сила потенциала

  inline double getRcut() const override { return r_cut_; }
  inline double getSqrRcut() const override { return r_cut_sqr_; }
};
#endif