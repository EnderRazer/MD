#ifndef BAROSTAT_H
#define BAROSTAT_H

#include <string>

#include "core/System.h"

/**
 * @brief Абстрактный класс для баростатов.
 *
 * Абстрактный класс для баростатов, который определяет методы для применения баростатов к системе.
 */
class Barostat {
public:
  /**
   * @brief Тип баростата.
   *
   * Тип баростата.
   */
  enum class BarostatType { NONE, BERENDSEN };
  /**
   * @brief Деструктор.
   */
  virtual ~Barostat() = default;
  /**
   * @brief Применение давления.
   */
  virtual void applyPressureControl(System &sys) = 0;
  /**
   * @brief Получение базовой информации о баростате.
   */
  virtual std::string getData() const = 0;
  /**
   * @brief Получение типа баростата.
   */
  inline virtual const BarostatType getBarostatType() const = 0;
};

#endif // BAROSTAT_H
