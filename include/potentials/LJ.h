#ifndef POTENTIAL_LJ
#define POTENTIAL_LJ

#include <sstream>

#include "nlohmann/json.hpp"

#include "core/Settings.h"
#include "core/System.h"

#include "Potential.h"

using json = nlohmann::json;

/**
 * @brief Класс для работы с потенциалом Леннард-Джонса.
 * @details Класс предоставляет методы для работы с потенциалом Леннард-Джонса.
 */
class LJ : public Potential {
private:

  /**
   * @brief Тип потенциала.
   */
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
  /**
   * @brief Деструктор.
   */
  ~LJ() = default;

  /**
   * @brief Конструктор.
   */
  LJ(json &params_json)
      : epsilon_(params_json.value("epsilon", 0.0)),
        sigma_(params_json.value("sigma", 0.0)), sigma_p2_(sigma_ * sigma_),
        sigma_p6_(sigma_p2_ * sigma_p2_ * sigma_p2_), r_cut_(2.5 * sigma_),
        r_cut_sqr_(r_cut_ * r_cut_), epsilon_t4_(4 * epsilon_),
        epsilon_t24_(24 * epsilon_), eps_t4_sigma_p6_(epsilon_t4_ * sigma_p6_),
        eps_t4_sigma_p12_(epsilon_t4_ * sigma_p6_ * sigma_p6_) {
    double r_p6 = r_cut_sqr_ * r_cut_sqr_ * r_cut_sqr_;
    double r_p12 = r_p6 * r_p6;
    r_cut_pot_ = eps_t4_sigma_p12_ / r_p12 - eps_t4_sigma_p6_ / r_p6;
  };

  /**
   * @brief Получение типа потенциала.
   * @return Тип потенциала.
   */
  inline PotentialType getPotentialType() const override {
    return type_;
  }

  /**
   * @brief Получение базовой информации о потенциале.
   * @return Базовая информация о потенциале.
   */
  std::string getData() const override {
    std::ostringstream oss;
    oss.precision(16);
    oss << "Using LJ potential!\n\tEpsilon: " << epsilon_
        << "\n\tSigma: " << sigma_ << "\n\t"
        << "Rcut: " << r_cut_ << "\n\tRcut_pot: " << r_cut_pot_ << "\n";
    return oss.str();
  }

  /**
   * @brief Получение потенциальной энергии.
   * @param r_p2 - расстояние между частицами.
   * @return Потенциальная энергия.
   */
  inline double getU(double r_p2) const override {
    double r_p6 = r_p2 * r_p2 * r_p2;
    return eps_t4_sigma_p12_ / (r_p6 * r_p6) - eps_t4_sigma_p6_ / r_p6 -
           r_cut_pot_;
  };

  /**
   * @brief Получение силы потенциала.
   * @param r_p2 - расстояние между частицами в квадрате.
   * @return Сила потенциала.
   */
  inline double getFU(double r_p2) const override {
    double r_p6 = r_p2 * r_p2 * r_p2;
    return 6 *
           (2 * eps_t4_sigma_p12_ / (r_p6 * r_p6) - eps_t4_sigma_p6_ / r_p6);
  }

  /**
   * @brief Получение всех результатов потенциала.
   * @param r_p2 - расстояние между частицами в квадрате.
   * @return Все результаты потенциала.
   */
  inline const PotentialResult getAll(double r_p2) const override {
    double r_p6 = r_p2 * r_p2 * r_p2;
    double eps_div_r_p12 = eps_t4_sigma_p12_ / (r_p6 * r_p6);
    double eps_div_r_p6 = eps_t4_sigma_p6_ / r_p6;

    double u = eps_div_r_p12 - eps_div_r_p6 - r_cut_pot_;
    double fu = 6 * (2 * eps_div_r_p12 - eps_div_r_p6);
    return {u, fu};
  }

  /**
   * @brief Получение радиуса обрезания.
   * @return Радиус обрезания.
   */
  inline double getRcut() const override { return r_cut_; }

  /**
   * @brief Получение квадрата радиуса обрезания.
   * @return Квадрат радиуса обрезания.
   */
  inline double getSqrRcut() const override { return r_cut_sqr_; }

  // EAM overrides, not implemented for LJ
  inline double getU(double rho_f, double mu) const override {
    throw std::runtime_error("getU(rho_f, mu) not implemented for LJ");
  };
  inline double getFU(double rho_f_i, double rho_f_j, double r) const override {
    throw std::runtime_error("getFU(rho_f_i, rho_f_j, r) not implemented for LJ");
  };
  inline double getPairPart(double r) const override {
    throw std::runtime_error("getPairPart(r) not implemented for LJ");
  };
  inline double getDensityPart(double r) const override {
    throw std::runtime_error("getDensityPart(r) not implemented for LJ");
  };
  inline Mu_RhoF getPairDesityPart(double r) const override {
    throw std::runtime_error("getPairDesityPart(r) not implemented for LJ");
  };
  inline double getDerPairPart(double r) const override {
    throw std::runtime_error("getDerPairPart(r) not implemented for LJ");
  };
  inline double getDerDensityPart(double r) const override {
    throw std::runtime_error("getDerDensityPart(r) not implemented for LJ");
  };
  inline double getCloud(double rho) const override {
    throw std::runtime_error("getCloud(rho) not implemented for LJ");
  }
  inline double getDerCloud(double rho) const override {
    throw std::runtime_error("getDerCloud(rho) not implemented for LJ");
  }
};
#endif
