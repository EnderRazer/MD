// Абстрактный класс потенциалов
#ifndef POTENTIAL_PARAMS_H
#define POTENTIAL_PARAMS_H

#include <string>

struct PotentialResult {
  double u{0.0};
  double fu{0.0};
};
class Potential {
public:
  enum class PotentialType { LJ, EAM, UNDEFINED }; // Strongly-typed enum
private:
  PotentialType type_{PotentialType::UNDEFINED};
  const double r_cut_{0.0};
  const double r_cut_sqr_{0.0};

public:
  virtual ~Potential() = 0;

  inline virtual double getU(double r) const = 0;  // Потенциальная энергия
  inline virtual double getFU(double r) const = 0; // Сила потенциала
  inline virtual const PotentialResult getAll(double r) const = 0;

  inline virtual double getU(double rho_f,
                             double mu) const = 0; // Потенциальная энергия
  inline virtual double getFU(double rho_f_i, double rho_f_j,
                              double r) const = 0; // Сила потенциала
  inline virtual double getDensityPart(double r) const = 0;
  inline virtual double getPairPart(double r) const = 0;
  inline virtual std::string getData() const = 0;

  inline virtual double getRcut() const = 0;
  inline virtual double getSqrRcut() const = 0;

  // inline virtual void precompute() = 0;

  inline virtual PotentialType getPotentialType() const = 0;
};

inline Potential::~Potential() {}

#endif // POTENTIAL_PARAMS_H
