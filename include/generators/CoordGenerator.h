#ifndef COORD_GENERATOR_H
#define COORD_GENERATOR_H

#include "classes/Vector3.h"
#include "core/Settings.h"
#include "core/System.h"

/**
 * @brief Структура для хранения пользовательской структуры материала.
 * @details Структура для хранения пользовательской структуры материала.
 */
struct CustomLayout {
  Vector3<double> coords;
  Vector3<double> velocities;
};

/**
 * @brief Класс для генерации координат частиц.
 * @details Класс предоставляет методы для генерации координат частиц.
 */
class CoordinatesGenerator {
private:
  /**
   * @brief Генерация решетки.
   * @param count - счетчик частиц.
   * @param sys - система частиц.
   * @param settings - настройки.
   * @param limit - ограничения.
   * @param offset - смещение.
   */
  void generateLattice(int &count, System &sys, Settings &settings,
                       const Vector3<int> &limit,
                       const Vector3<double> &offset = {0, 0, 0}) {
    for (int x = 0; x < limit.x(); x++) {
      for (int y = 0; y < limit.y(); y++) {
        for (int z = 0; z < limit.z(); z++) {
          Vector3<double> coord(x + offset.x(), y + offset.y(), z + offset.z());
          sys.particles()[count].setCoord(coord *
                                          sys.dimensions().cristLength());
          count++;
        }
      }
    }
  }

  /**
   * @brief Генерация куба.
   * @param sys - система частиц.
   * @param settings - настройки.
   */
  void generatePrimitiveCube(System &sys, Settings &settings) {
    Vector3<int> limit = settings.hasPbc() ? sys.dimensions().numCrists()
                                           : sys.dimensions().numCrists() + 1;
    int count = 0;
    generateLattice(count, sys, settings, limit);
    assert(count == sys.particleNumber() && "Particle number mismatch");
  }

  /**
   * @brief Генерация FCC структуры.
   * @param sys - система частиц.
   * @param settings - настройки.
   */
  void generateFCC(System &sys, Settings &settings) {
    Vector3<int> offset = sys.dimensions().numVoid()/2;
    Vector3<int> limit = settings.hasPbc() ? sys.dimensions().numCrists()
                                           : sys.dimensions().numCrists() + 1;

    int count = 0;

    // FCC: Vertices
    generateLattice(count, sys, settings, limit, offset);
    std::cout << "||Vertices|| Generated: " << count << std::endl;
    // FCC: Face centers
    generateLattice(
        count, sys, settings,
        Vector3<int>{limit.x() - 1, limit.y() - 1, limit.z()},
        Vector3<double>{offset.x() + 0.5, offset.y() + 0.5, offset.z() + 0.0});
    std::cout << "||XY Face centers|| Generated: " << count << std::endl;
    generateLattice(
        count, sys, settings,
        Vector3<int>{limit.x() - 1, limit.y(), limit.z() - 1},
        Vector3<double>{offset.x() + 0.5, offset.y() + 0.0, offset.z() + 0.5});
    std::cout << "||XZ Face centers|| Generated: " << count << std::endl;
    generateLattice(
        count, sys, settings,
        Vector3<int>{limit.x(), limit.y() - 1, limit.z() - 1},
        Vector3<double>{offset.x() + 0.0, offset.y() + 0.5, offset.z() + 0.5});
    std::cout << "||YZ Face centers|| Generated: " << count << std::endl;

    assert(count == sys.particleNumber() && "Particle number mismatch");
  }

  /**
   * @brief Генерация пользовательской структуры.
   * @param sys - система частиц.
   * @param settings - настройки.
   */
  void generateCustom(System &sys, Settings &settings) {
    json custom = settings.customLayout();
    std::cout << custom << std::endl;
    if (custom.empty())
      throw std::invalid_argument("Custom layout is empty");
    std::vector<CustomLayout> customLayouts;
    for (const auto &item : custom["particles"]) {
      std::cout << item << std::endl;
      CustomLayout layout;
      layout.coords = item["position"].get<std::vector<double>>();
      layout.velocities = item["velocity"].get<std::vector<double>>();
      customLayouts.push_back(layout);
    }
    sys.particles().resize(customLayouts.size());
    for (int i = 0; i < customLayouts.size(); i++) {
      sys.particles()[i].setCoord(customLayouts[i].coords);
      sys.particles()[i].setVelocity(customLayouts[i].velocities);
    }
    std::cout << sys.getParticlesCoordsInfo() << std::endl;
    std::cout << sys.getParticlesVelocityInfo() << std::endl;
  }

public:
  CoordinatesGenerator() = default;
  ~CoordinatesGenerator() = default;

  // Запрещаем копирование
  CoordinatesGenerator(const CoordinatesGenerator &) = delete;
  CoordinatesGenerator &operator=(const CoordinatesGenerator &) = delete;

  /**
   * @brief Генерация начальных координат.
   * @param sys - система частиц.
   * @param settings - настройки.
   */
  void startingCoords(System &sys, Settings &settings) {
    if (settings.structType() == std::string("CP")) {
      std::cout << "Generating CP structure" << std::endl;
      generatePrimitiveCube(sys, settings);
    } else if (settings.structType() == std::string("FCC")) {
      std::cout << "Generating FCC structure" << std::endl;
      generateFCC(sys, settings);
    } else if (settings.structType() == std::string("Custom")) {
      std::cout << "Generating custom structure" << std::endl;
      generateCustom(sys, settings);
    } else {
      throw std::runtime_error("Unknown structure type: " +
                               settings.structType());
    }
  }
};

#endif
