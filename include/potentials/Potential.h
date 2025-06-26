// Абстрактный класс потенциалов
#ifndef POTENTIAL_PARAMS_H
#define POTENTIAL_PARAMS_H

#include <string>

/**
 * @brief Тип потенциала.
 */
enum class PotentialType { NONE, LJ, EAM };

/**
 * @brief Структура для хранения результатов потенциала.
 */
struct PotentialResult {
  double u{0.0};
  double fu{0.0};
};

/**
 * @brief Структура для хранения результатов потенциала.
 */
struct Mu_RhoF {
  double mu;
  double rho_f;
};

/**
 * @brief Абстрактный класс потенциалов.
 * @details Класс предоставляет методы для работы с потенциалами.
 */
class Potential {
private:
  /**
   * @brief Тип потенциала.
   */
  PotentialType type_{PotentialType::NONE};

  /**
   * @brief Радиус обрезания.
   */
  const double r_cut_{0.0};

  /**
   * @brief Квадрат радиуса обрезания.
   */
  const double r_cut_sqr_{0.0};

public:
  /**
   * @brief Деструктор.
   */
  virtual ~Potential() = 0;

  /**
   * @brief Получение потенциальной энергии.
   */
  inline virtual double getU(double r) const = 0;

  /**
   * @brief Получение силы потенциала.
   */
  inline virtual double getFU(double r) const = 0;

  /**
   * @brief Получение всех результатов потенциала.
   */
  inline virtual const PotentialResult getAll(double r) const = 0;

  /**
   * @brief Получение потенциальной энергии.
   */
  inline virtual double getU(double rho_f, double mu) const = 0;

  /**
   * @brief Получение силы потенциала.
   */
  inline virtual double getFU(double rho_f_i, double rho_f_j, double r) const = 0;

  /**
   * @brief Получение электронной плотности.
   */
  inline virtual double getDensityPart(double r) const = 0;

  /**
   * @brief Получение парного потенциала.
   */
  inline virtual double getPairPart(double r) const = 0;

  /**
   * @brief Получение электронной плотности и парного потенциала.
   */
  inline virtual Mu_RhoF getPairDesityPart(double r) const = 0;

  inline virtual double getDerPairPart(double r) const = 0;
  inline virtual double getDerDensityPart(double r) const = 0;

  inline virtual double getCloud(double rho) const = 0;
  inline virtual double getDerCloud(double rho) const = 0;

/**
   * @brief Получение базовой информации о потенциале.
   * @return Базовая информация о потенциале.
   */
  inline virtual std::string getData() const = 0;

  /**
   * @brief Получение радиуса обрезания.
   */
  inline virtual double getRcut() const = 0;

  /**
   * @brief Получение квадрата радиуса обрезания.
   */
  inline virtual double getSqrRcut() const = 0;

  /**
   * @brief Получение типа потенциала.
   */
  inline virtual PotentialType getPotentialType() const = 0;
};

inline Potential::~Potential() {}

#endif // POTENTIAL_PARAMS_H
