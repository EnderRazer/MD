#ifndef COORD_GENERATOR_H
#define COORD_GENERATOR_H

struct CustomLayout {
  Vector3<double> coords;
  Vector3<double> velocities;
};

class CoordinatesGenerator {
private:
  void generateLattice(int &count, System &sys, Settings &settings,
                       const Vector3<int> &limit,
                       const Vector3<double> &offset = {0, 0, 0}) {
    for (int x = 0; x < limit.x(); x++) {
      for (int y = 0; y < limit.y(); y++) {
        for (int z = 0; z < limit.z(); z++) {
          Vector3<double> coord(x + offset.x(), y + offset.y(), z + offset.z());
          sys.particles().at(count).setCoord(coord *
                                             sys.dimensions().cristLength());
          count++;
        }
      }
    }
  }

  void generatePrimitiveCube(System &sys, Settings &settings) {
    Vector3<int> limit = settings.hasPbc() ? sys.dimensions().numCrists()
                                           : sys.dimensions().numCrists() + 1;
    int count = 0;
    generateLattice(count, sys, settings, limit);
    assert(count == sys.particleNumber() && "Particle number mismatch");
  }
  void generateFCC(System &sys, Settings &settings) {
    Vector3<int> limit = settings.hasPbc() ? sys.dimensions().numCrists()
                                           : sys.dimensions().numCrists() + 1;

    int count = 0;

    // FCC: Vertices
    generateLattice(count, sys, settings, limit);

    // FCC: Face centers
    generateLattice(count, sys, settings, limit,
                    Vector3<double>{0.5, 0.5, 0.0});
    generateLattice(count, sys, settings, limit,
                    Vector3<double>{0.5, 0.0, 0.5});
    generateLattice(count, sys, settings, limit,
                    Vector3<double>{0.0, 0.5, 0.5});

    assert(count == sys.particleNumber() && "Particle number mismatch");
  }

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
  CoordinatesGenerator(/* args */);
  ~CoordinatesGenerator();

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

CoordinatesGenerator::CoordinatesGenerator(/* args */) {}

CoordinatesGenerator::~CoordinatesGenerator() {}
#endif