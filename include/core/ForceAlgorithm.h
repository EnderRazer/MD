#ifndef FORCE_ALGORITHM_H
#define FORCE_ALGORITHM_H

class ForceAlgorithm {
private:
  Settings &settings_;
  Dimensions dim_;
  std::unique_ptr<Potential> potential_;

  // Предрасчеты для EAM
  std::vector<double> rho_;
  std::vector<double> mu_;

  inline Vector3<double> mirror_vector(Vector3<double> &vec, double lx,
                                       double ly, double lz) {
    double halfLX = lx / 2;
    double halfLY = ly / 2;
    double halfLZ = lz / 2;

    if (vec.x() > halfLX)
      vec.x() -= lx;
    if (vec.x() <= -halfLX)
      vec.x() += lx;

    if (vec.y() > halfLY)
      vec.y() -= ly;
    if (vec.y() <= -halfLY)
      vec.y() += ly;

    if (vec.z() > halfLZ)
      vec.z() -= lz;
    if (vec.z() <= -halfLZ)
      vec.z() += lz;
    return vec;
  }

public:
  inline ForceCalcValues compute(Particle &p1, Particle &p2) {
    ForceCalcValues output;
    Vector3<double> rVec = p1.coord() - p2.coord();
    output.rVec = rVec;
    // Отражение частицы
    if (settings_.hasPbc()) {
      double lx = dim_.lx();
      double ly = dim_.ly();
      double lz = dim_.lz();
      rVec = mirror_vector(rVec, lx, ly, lz);
    }
    output.rVec_mirrored = rVec;
    double lengthSqr = rVec.lengthSquared();
    // Проверка на соседство
    if (lengthSqr > potential_->getSqrRcut())
      return output;
    // Расчет силы
    double U = 0.0;
    double FU = 0.0;
    switch (potential_->getPotentialType()) {
    case Potential::PotentialType::LJ:
      U = potential_->getU(lengthSqr);
      FU = potential_->getFU(lengthSqr);
      break;
    case Potential::PotentialType::EAM:
      U = potential_->getU(lengthSqr);
      FU = potential_->getFU(lengthSqr);
    default:
      std::cout << "Unknown type of potential" << std::endl;
      break;
    }
    output.e_pot = U / 2;
    output.force = rVec * FU / lengthSqr;
    output.virials = Matrix3().outerProduct(rVec, output.force);
    return output;
  }
  ForceAlgorithm() = delete;
  ForceAlgorithm(Settings &settings, std::unique_ptr<Potential> &&potential,
                 Dimensions &dim)
      : settings_(settings), potential_(std::move(potential)), dim_(dim) {}
  ~ForceAlgorithm() = default;

  // Запрещаем копирование
  ForceAlgorithm(const ForceAlgorithm &) = delete;
  ForceAlgorithm &operator=(const ForceAlgorithm &) = delete;
};
#endif // FORCE_ALGORITHM_H