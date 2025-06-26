#ifndef CELL_LIST_H
#define CELL_LIST_H

#include <cassert>
#include <iostream>
#include <sstream>

#include "classes/Particle.h"
#include "classes/Vector3.h"

/**
 * @brief Класс для списка ячеек.
 *
 * Класс для списка ячеек, который определяет методы для построения списка ячеек и получения соседей.
 */
class CellList {
public:
  /**
   * @brief Конструктор.
   *
   * Конструктор для списка ячеек.
   * @param cutoff - радиус обрезания.
   * @param boxSize - размер ячейки.
   */
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

  /**
   * @brief Построение списка ячеек.
   *
   * Построение списка ячеек.
   * @param particles - вектор частиц.
   */
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

  /**
   * @brief Получение списка соседей.
   *
   * Получение списка соседей для частицы.
   * @param particles - вектор частиц.
   * @param index - индекс частицы.
   * @return вектор соседей.
   */
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

          // NEW with distance check aka Verlet List
          for (int i = 0; i < cells_[neighborIdx].size(); ++i) {
            if (index != cells_[neighborIdx][i]) {
              Vector3<double> rVec = particles[cells_[neighborIdx][i]].coord() -
                                     particles[index].coord();
              
              // Mirrorig vector (PBC)
              if (rVec.x() > boxSize_.x() / 2)
                rVec.x() -= boxSize_.x();
              if (rVec.x() <= -boxSize_.x() / 2)
                rVec.x() += boxSize_.x();

              if (rVec.y() > boxSize_.y() / 2)
                rVec.y() -= boxSize_.y();
              if (rVec.y() <= -boxSize_.y() / 2)
                rVec.y() += boxSize_.y();

              if (rVec.z() > boxSize_.z() / 2)
                rVec.z() -= boxSize_.z();
              if (rVec.z() <= -boxSize_.z() / 2)
                rVec.z() += boxSize_.z();

              if (rVec.length() < cutoff_)
                neighbors.push_back(cells_[neighborIdx][i]);
            }
          }
          // OLD without distance check
          // neighbors.insert(neighbors.end(), cells_[neighborIdx].begin(),
          //                  cells_[neighborIdx].end());
        }
      }
    }
    return neighbors;
  }

  /**
   * @brief Получение базовой информации о списке ячеек.
   *
   * Получение базовой информации о списке ячеек.
   * @return строка с информацией о списке ячеек.
   */
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

  /**
   * @brief Получение количества ячеек.
   *
   * Получение количества ячеек.
   * @return количество ячеек.
   */
  inline const int getNumCells() const { return totalCells_; };

private:
  /**
   * @brief Радиус обрезания.
   */
  double cutoff_;

  /**
   * @brief Размер бокса.
   */
  Vector3<double> boxSize_;

  /**
   * @brief Размер ячейки.
   */
  Vector3<double> cellSize_;

  /**
   * @brief Количество ячеек.
   */
  Vector3<int> numCells_;

  /**
   * @brief Общее количество ячеек.
   */
  int totalCells_;
  
  /**
   * @brief Список ячеек.
   */
  std::vector<std::vector<int>> cells_;

  /**
   * @brief Очистка списка ячеек.
   */
  void clearCells() {
    for (auto &cell : cells_) {
      cell.clear();
    }
  }

  /**
   * @brief Получение индекса ячейки.
   *
   * Получение индекса ячейки для частицы.
   * @param pos - координаты частицы.
   * @return индекс ячейки.
   */
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
