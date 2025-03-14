// Абстрактный класс потенциалов
#ifndef POTENTIAL_PARAMS_H
#define POTENTIAL_PARAMS_H

struct LJResult
{
  double u{0.0};
  double fu{0.0};
};
class Potential
{
public:
  enum class PotentialType
  {
    LJ,
    EAM,
    UNDEFINED
  }; // Strongly-typed enum
private:
  PotentialType type_{PotentialType::UNDEFINED};
  const double r_cut_{0.0};
  const double r_cut_sqr_{0.0};

public:
  virtual ~Potential() = 0;

  inline virtual double getU(double r) const = 0;  // Потенциальная энергия
  inline virtual double getFU(double r) const = 0; // Сила потенциала
  inline virtual const LJResult getAll(double r) const = 0;

  inline virtual double getU(double rho_f,
                             double mu) const = 0; // Потенциальная энергия
  inline virtual double getFU(double rho_f_i, double rho_f_j, double d_rho_f_ij,
                              double d_mu_ij) const = 0; // Сила потенциала
  inline virtual std::string getData() const = 0;

  inline virtual double getRcut() const = 0;
  inline virtual double getSqrRcut() const = 0;

  // inline virtual void precompute() = 0;

  inline virtual PotentialType getPotentialType() const = 0;
};

inline Potential::~Potential() {}

#endif // POTENTIAL_PARAMS_H