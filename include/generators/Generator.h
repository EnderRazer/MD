#ifndef GENERATOR_H
#define GENERATOR_H

#include "CoordGenerator.h"
#include "SpeedGenerator.h"

using json = nlohmann::json;

class Generator {
private:
  SpeedGenerator speedGenerator_;
  CoordinatesGenerator coordGenerator_;

public:
  void generateCoords(System &sys, Settings &settings) {
    coordGenerator_.startingCoords(sys, settings);
  }
  void generateSpeeds(System &sys, Settings &settings) {
    speedGenerator_.startingSpeeds(sys, settings);
  }
  Generator(Settings &settings)
      : speedGenerator_(SpeedGenerator(settings)),
        coordGenerator_(CoordinatesGenerator()) {}
  ~Generator() = default;

  // Запрещаем копирование
  Generator(const Generator &) = delete;
  Generator &operator=(const Generator &) = delete;
};
#endif