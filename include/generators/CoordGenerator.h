#ifndef COORD_GENERATOR_H
#define COORD_GENERATOR_H

#include "classes/Particles.h"
#include "core/Settings.h"
#include "core/System.h"
#include <vector>

struct CustomLayout {
  double coord_x{}, coord_y{}, coord_z{};
  double velocity_x{}, velocity_y{}, velocity_z{};
};

class CoordinatesGenerator {
private:
  void generateLattice(int &count, System &sys, Settings &settings,
                       const int &limit_x, const int &limit_y, const int &limit_z,
                       const double &offset_x = 0.0, const double &offset_y = 0, const double &offset_z = 0) {
    double coord_x = 0, coord_y = 0, coord_z = 0;
    Particles &particles = sys.particles();
    double &crist_length = sys.dimensions().cristLength();
    for (int x = 0; x < limit_x; x++) {
      for (int y = 0; y < limit_y; y++) {
        for (int z = 0; z < limit_z; z++) {
          coord_x = x + offset_x;
          coord_y = y + offset_y;
          coord_z = z + offset_z;
          
          particles.coordX(count) = coord_x*crist_length;
          particles.coordY(count) = coord_y*crist_length;
          particles.coordZ(count) = coord_z*crist_length;
          count++;
        }
      }
    }
  }

  void generatePrimitiveCube(System &sys, Settings &settings) {
    int limit_x = settings.hasPbc() ? sys.dimensions().numCristX()
                                           : sys.dimensions().numCristX() + 1;
    int limit_y = settings.hasPbc() ? sys.dimensions().numCristY()
                                           : sys.dimensions().numCristY() + 1;
    int limit_z = settings.hasPbc() ? sys.dimensions().numCristZ()
                                           : sys.dimensions().numCristZ() + 1;

    int count = 0;
    generateLattice(count, sys, settings, limit_x, limit_y,limit_z);
    assert(count == sys.particleNumber() && "Particle number mismatch");
  }
  void generateFCC(System &sys, Settings &settings) {
    int offset_x = 0.5 * sys.dimensions().numVoidX();
    int offset_y = 0.5 * sys.dimensions().numVoidX();
    int offset_z = 0.5 * sys.dimensions().numVoidX();
    
    int limit_x = settings.hasPbc() ? sys.dimensions().numCristX()
                                           : sys.dimensions().numCristX() + 1;
    int limit_y = settings.hasPbc() ? sys.dimensions().numCristY()
                                           : sys.dimensions().numCristY() + 1;
    int limit_z = settings.hasPbc() ? sys.dimensions().numCristZ()
                                           : sys.dimensions().numCristZ() + 1;

    int count = 0;

    // FCC: Vertices
    generateLattice(count, sys, settings, limit_x, limit_y, limit_z, offset_x, offset_y, offset_z);
    std::cout << "||Vertices|| Generated: " << count << std::endl;
    // FCC: Face centers
    generateLattice(
        count, sys, settings,
        limit_x - 1, limit_y - 1, limit_z,
        offset_x + 0.5, offset_y + 0.5, offset_z + 0.0);
    std::cout << "||XY Face centers|| Generated: " << count << std::endl;
    generateLattice(
        count, sys, settings,
        limit_x - 1, limit_y, limit_z - 1,
        offset_x + 0.5, offset_y + 0.0, offset_z + 0.5);
    std::cout << "||XZ Face centers|| Generated: " << count << std::endl;
    generateLattice(
        count, sys, settings,
        limit_x, limit_y - 1, limit_z - 1,
        offset_x + 0.0, offset_y + 0.5, offset_z + 0.5);
    std::cout << "||YZ Face centers|| Generated: " << count << std::endl;

    assert(count == sys.particleNumber() && "Particle number mismatch");

    //std::cout << sys.getParticlesCoordsInfo() << std::endl;
  }

  void generateCustom(System &sys, Settings &settings) {
    json custom = settings.customLayout();
    std::cout << custom << std::endl;
    if (custom.empty())
      throw std::invalid_argument("Custom layout is empty");
    std::vector<CustomLayout> layouts;
    CustomLayout layout;
    for (const auto &item : custom["particles"]) {
      std::cout << item << std::endl;
      std::vector<double> positions = item["position"].get<std::vector<double>>();
      std::vector<double> velocities = item["velocity"].get<std::vector<double>>();

      layout.coord_x = positions[0];
      layout.coord_x = positions[1];
      layout.coord_x = positions[2];

      layout.velocity_x = velocities[0];
      layout.velocity_x = velocities[1];
      layout.velocity_x = velocities[2];
      
      layouts.push_back(layout);
    }
    int size = layouts.size();
    sys.particles().resize(size);
    Particles &particles = sys.particles();
    for (int i = 0; i < size; i++) {
      particles.coordX(i) = layouts[i].coord_x;
      particles.coordY(i) = layouts[i].coord_y;
      particles.coordZ(i) = layouts[i].coord_z;

      particles.velocityX(i) = layouts[i].velocity_x;
      particles.velocityY(i) = layouts[i].velocity_y;
      particles.velocityZ(i) = layouts[i].velocity_z;
    }
    
    std::cout << sys.getParticlesCoordsInfo() << std::endl;
    std::cout << sys.getParticlesVelocityInfo() << std::endl;
  }

public:
  CoordinatesGenerator(/* args */) = default;
  ~CoordinatesGenerator() = default;

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
  // Запрещаем копирование
  CoordinatesGenerator(const CoordinatesGenerator &) = delete;
  CoordinatesGenerator &operator=(const CoordinatesGenerator &) = delete;
};

#endif
