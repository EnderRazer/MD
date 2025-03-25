#ifndef CELL_LIST_H
#define CELL_LIST_H

#include <cassert>
#include <iostream>
#include <sstream>

#include "classes/Particle.h"
#include "classes/Vector3.h"

class CellList {
public:
  CellList(double cutoff, const Vector3<double> &boxSize)
      : cutoff_(cutoff), boxSize_(boxSize) {
    numCells_.x() =
        std::max(1, static_cast<int>(std::floor(boxSize_.x() / cutoff_)));
    numCells_.y() =
        std::max(1, static_cast<int>(std::floor(boxSize_.y() / cutoff_)));
    numCells_.z() =
        std::max(1, static_cast<int>(std::floor(boxSize_.z() / cutoff_)));
    cellSize_.x() = boxSize_.x() / numCells_.x();
    cellSize_.y() = boxSize_.y() / numCells_.y();
    cellSize_.z() = boxSize_.z() / numCells_.z();
    totalCells_ = numCells_.x() * numCells_.y() * numCells_.z();
    cells_.resize(totalCells_);
  }

  void build(const std::vector<Particle> &particles) {
    clearCells();
    for (int i = 0; i < particles.size(); ++i) {
      int cellIdx = getCellIndex(particles[i].coord());
      if (!(cellIdx >= 0 && cellIdx < totalCells_)) {
        std::cout << "Cell index out of bounds: " << cellIdx
                  << " For particle: " << i << " " << particles[i].coord()
                  << std::endl;
        throw std::length_error("Cell index out of bounds");
      }
      cells_[cellIdx].push_back(i);
    }
  }

  std::vector<int> getNeighbors(const std::vector<Particle> &particles,
                                const int index) const {
    std::vector<int> neighbors;
    int cellIdx = getCellIndex(particles[index].coord());
    int cz = cellIdx / (numCells_.x() * numCells_.y());
    int cy = (cellIdx % (numCells_.x() * numCells_.y())) / numCells_.x();
    int cx = (cellIdx % (numCells_.x() * numCells_.y())) % numCells_.x();

    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dz = -1; dz <= 1; ++dz) {
          int nx = (cx + dx + numCells_.x()) % numCells_.x();
          int ny = (cy + dy + numCells_.y()) % numCells_.y();
          int nz = (cz + dz + numCells_.z()) % numCells_.z();
          int neighborIdx =
              nx + ny * numCells_.x() + nz * numCells_.x() * numCells_.y();
          assert(neighborIdx >= 0 && neighborIdx < totalCells_);
          /*
          for (int i = 0; i < cells_[neighborIdx].size(); ++i) {
            if (index != cells_[neighborIdx][i]) {
              Vector3<double> rVec = particles[cells_[neighborIdx][i]].coord() -
                                     particles[index].coord();
              if (rVec.length() < cutoff_)
                neighbors.push_back(cells_[neighborIdx][i]);
            }
          }
          */
          neighbors.insert(neighbors.end(), cells_[neighborIdx].begin(),
                           cells_[neighborIdx].end());
        }
      }
    }
    return neighbors;
  }

  inline const std::string getData() const {
    std::ostringstream oss;
    oss << "CellList Data:";
    oss << "\n\tCutoff: " << cutoff_;
    oss << "\n\tBox Size: " << boxSize_;
    oss << "\n\tCell Size: " << cellSize_;
    oss << "\n\tTotal Cells: " << totalCells_;
    oss << "\n\tNumber of Cells: " << numCells_;

    return oss.str();
  }

  inline const int getNumCells() const { return totalCells_; };

private:
  double cutoff_;
  Vector3<double> boxSize_;
  Vector3<double> cellSize_;
  Vector3<int> numCells_;
  int totalCells_;
  std::vector<std::vector<int>> cells_;

  void clearCells() {
    for (auto &cell : cells_) {
      cell.clear();
    }
  }

  int getCellIndex(const Vector3<double> &pos) const {
    int x = std::min(static_cast<int>(std::floor(pos.x() / cellSize_.x())),
                     numCells_.x() - 1);
    int y = std::min(static_cast<int>(std::floor(pos.y() / cellSize_.y())),
                     numCells_.y() - 1);
    int z = std::min(static_cast<int>(std::floor(pos.z() / cellSize_.z())),
                     numCells_.z() - 1);
    return x + y * numCells_.x() + z * numCells_.x() * numCells_.y();
  }
};

#endif
