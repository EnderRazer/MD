#ifndef GENERATOR_H
#define GENERATOR_H

#include "CoordGenerator.h"
#include "SpeedGenerator.h"

using json = nlohmann::json;

/**
 * @brief Класс для генерации координат и скоростей частиц.
 * @details Класс предоставляет методы для генерации координат и скоростей частиц.
 */
class Generator {
private:
  /**
   * @brief Генератор скоростей.
   */
  SpeedGenerator speedGenerator_;

  /**
   * @brief Генератор координат.
   */
  CoordinatesGenerator coordGenerator_;

public:
  Generator() = default;
  ~Generator() = default;

  // Запрещаем копирование
  Generator(const Generator &) = delete;
  Generator &operator=(const Generator &) = delete;

/**
   * @brief Конструктор.
   * @param settings - настройки.
   */
  Generator(Settings &settings)
      : speedGenerator_(SpeedGenerator(settings)),
        coordGenerator_(CoordinatesGenerator()) {}
        
  /**
   * @brief Генерация координат.
   * @param sys - система частиц.
   * @param settings - настройки.
   */
  void generateCoords(System &sys, Settings &settings) {
    coordGenerator_.startingCoords(sys, settings);
  }

  /**
   * @brief Генерация скоростей.
   * @param sys - система частиц.
   * @param settings - настройки.
   */
  void generateSpeeds(System &sys, Settings &settings) {
    speedGenerator_.startingSpeeds(sys, settings);
  }
  
};
#endif