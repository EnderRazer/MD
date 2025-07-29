#ifndef COORDINATE_ALGORITHM_H
#define COORDINATE_ALGORITHM_H

#include "classes/Dimensions.h"
#include "classes/Particles.h"
#include "classes/Timer.h"
#include "core/Settings.h"

class CoordinateAlgorithm {
private:
  Settings &settings_;
  double dt_{0.0};


public:
  inline void compute(Particles &p) {
    Timer timer{1};
    timer.start();
    p.computeCoordinatesSIMD(dt_);
    timer.stop();
    std::cout << "Coord time elapsed = " << timer.elapsed() << "ms"<<std::endl;
  }
  
  inline void applyPBC(Dimensions &dim, Particles &p) {
    const double lx = dim.sizeX();
    const double ly = dim.sizeY();
    const double lz = dim.sizeZ();
    for(int i=0;i<p.size();i++){
      p.coordX(i) -= lx * std::floor(p.coordX(i) / lx);
      p.coordY(i) -= ly * std::floor(p.coordY(i) / ly);
      p.coordZ(i) -= lz * std::floor(p.coordZ(i) / lz);
    }
  }
  CoordinateAlgorithm() = delete;
  CoordinateAlgorithm(Settings &settings) : settings_(settings), dt_(settings_.dt()) {}
  ~CoordinateAlgorithm() = default;

  // Запрещаем копирование
  CoordinateAlgorithm(const CoordinateAlgorithm &) = delete;
  CoordinateAlgorithm &operator=(const CoordinateAlgorithm &) = delete;
};
#endif // COORDINATE_ALGORITHM_H
