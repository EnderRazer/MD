#ifndef CELL_LIST_H
#define CELL_LIST_H

#include <_abort.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <ostream>
#include <sstream>
#include <vector>

#include "classes/Particles.h"

class CellList {
public:
  CellList() = default;

  void rebuild(const double &cutoff, const double &boxSize_x,
               const double &boxSize_y, const double &boxSize_z) {
    cutoff_ = cutoff;
    cutoff_sqr_ = cutoff * cutoff;

    box_size_x_ = boxSize_x, box_size_y_ = boxSize_y, box_size_z_ = boxSize_z;
    half_box_size_x_ = 0.5 * boxSize_x, half_box_size_y_ = 0.5 * boxSize_y,
    half_box_size_z_ = 0.5 * boxSize_z;

    inv_box_size_x_ = 1 / boxSize_x, inv_box_size_y_ = 1 / boxSize_y,
    inv_box_size_z_ = 1 / boxSize_z;

    num_cells_x_ = std::max(1, static_cast<int>(box_size_x_ / cutoff_));
    num_cells_y_ = std::max(1, static_cast<int>(box_size_y_ / cutoff_));
    num_cells_z_ = std::max(1, static_cast<int>(box_size_z_ / cutoff_));
    num_cells_xy_ = num_cells_x_ * num_cells_y_;

    cell_size_x_ = box_size_x_ / num_cells_x_,
    cell_size_y_ = box_size_y_ / num_cells_y_,
    cell_size_z_ = box_size_z_ / num_cells_z_;
    inv_cell_size_x_ = 1 / cell_size_x_, inv_cell_size_y_ = 1 / cell_size_y_,
    inv_cell_size_z_ = 1 / cell_size_z_;

    totalCells_ = num_cells_x_ * num_cells_y_ * num_cells_z_;
  }
  void distribute_particles(const Particles &particles) {
    // 1. Очистка
    flat_cells_.clear();
    cell_offsets_.clear();

    flat_cells_.resize(particles.size());
    cell_offsets_.resize(totalCells_ + 1); // +1 для удобства границы

    cell_counts_.assign(totalCells_, 0);

    // 2. Подсчет сколько частиц в каждой ячейке
    for (int i = 0; i < particles.size(); ++i) {
      int cellIdx = getCellIndex(particles.coord_x_[i], particles.coord_y_[i],
                                 particles.coord_z_[i]);
      ++cell_counts_[cellIdx];
    }
    // 3. Построение offset-таблицыf
    cell_offsets_[0] = 0;
    for (int i = 0; i < totalCells_; ++i) {
      cell_offsets_[i + 1] = cell_offsets_[i] + cell_counts_[i];
    }

    // 4. Заполнение частиц
    std::vector<size_t> current_offsets =
        cell_offsets_; // текущая позиция вставки
    for (int i = 0; i < particles.size(); ++i) {
      int cellIdx = getCellIndex(particles.coord_x_[i], particles.coord_y_[i],
                                 particles.coord_z_[i]);
      flat_cells_[current_offsets[cellIdx]++] = i;
    }
  }

  inline int maxParticlesInCell() {
    return *std::max_element(cell_counts_.begin(), cell_counts_.end());
  }
  inline int minParticlesInCell() {
    return *std::min_element(cell_counts_.begin(), cell_counts_.end());
  }
  inline int avgParticlesInCell() {
    double avg = 0;
    for (int i = 0; i < cell_counts_.size(); i++) {
      avg += cell_counts_[i];
    }
    return avg / cell_counts_.size();
  }

  int getNeighbors(std::vector<int> &neighbors, Particles &particles,
                   const int index, bool pbc) const {
    neighbors.clear();
    int count = 0;
    double xi = particles.coord_x_[index];
    double yi = particles.coord_y_[index];
    double zi = particles.coord_z_[index];

    double *__restrict__ xn = particles.coord_x_.data();
    double *__restrict__ yn = particles.coord_y_.data();
    double *__restrict__ zn = particles.coord_z_.data();

    int cellIdx = getCellIndex(xi, yi, zi);
    int cz = cellIdx / num_cells_xy_;
    int cy = (cellIdx - cz * num_cells_xy_) / num_cells_x_;
    int cx = cellIdx - cy * num_cells_x_ - cz * num_cells_xy_;

    // Предварительно вычислим 27 соседних ячеек
    int neighbor_cells[27];
    int n_cells = 0;
    for (int dx = -1; dx <= 1; ++dx)
      for (int dy = -1; dy <= 1; ++dy)
        for (int dz = -1; dz <= 1; ++dz) {
          int nx = cx + dx, ny = cy + dy, nz = cz + dz;
          if (pbc) {
            nx = (nx + num_cells_x_) % num_cells_x_;
            ny = (ny + num_cells_y_) % num_cells_y_;
            nz = (nz + num_cells_z_) % num_cells_z_;
          } else {
            if (nx >= num_cells_x_ || nx < 0 || ny >= num_cells_y_ || ny < 0 ||
                nz >= num_cells_z_ || nz < 0)
              continue;
          }
          neighbor_cells[n_cells] = nx + ny * num_cells_x_ + nz * num_cells_xy_;
          n_cells++;
        }
    for (int ni = 0; ni < n_cells; ++ni) {
      int nId = neighbor_cells[ni];
      int start = cell_offsets_[nId], end = cell_offsets_[nId + 1];

      for (int k = start; k < end; ++k) {
        int j = flat_cells_[k];
        if (j == index)
          continue;

        double dx = xn[j] - xi, dy = yn[j] - yi, dz = zn[j] - zi;

        if (pbc) {
          dx -= std::round(dx / box_size_x_) * box_size_x_;
          dy -= std::round(dy / box_size_y_) * box_size_y_;
          dz -= std::round(dz / box_size_z_) * box_size_z_;
        }

        double r2 = dx * dx + dy * dy + dz * dz;
        if (r2 < cutoff_sqr_) {
          neighbors.push_back(j);
          count++;
        }
      }
    }

    return count;
  }

  int getNeighborsUniq(std::vector<int> &neighbors, Particles &particles,
                       const int index, bool pbc) const {
    neighbors.clear();
    int count = 0;
    double xi = particles.coord_x_[index];
    double yi = particles.coord_y_[index];
    double zi = particles.coord_z_[index];

    double *__restrict__ xn = particles.coord_x_.data();
    double *__restrict__ yn = particles.coord_y_.data();
    double *__restrict__ zn = particles.coord_z_.data();

    int cellIdx = getCellIndex(xi, yi, zi);
    int cz = cellIdx / num_cells_xy_;
    int cy = (cellIdx - cz * num_cells_xy_) / num_cells_x_;
    int cx = cellIdx - cy * num_cells_x_ - cz * num_cells_xy_;

    // Предварительно вычислим 27 соседних ячеек
    int neighbor_cells[27];
    int n_cells = 0;
    for (int dx = -1; dx <= 1; ++dx)
      for (int dy = -1; dy <= 1; ++dy)
        for (int dz = -1; dz <= 1; ++dz) {
          int nx = cx + dx, ny = cy + dy, nz = cz + dz;
          if (pbc) {
            nx = (nx + num_cells_x_) % num_cells_x_;
            ny = (ny + num_cells_y_) % num_cells_y_;
            nz = (nz + num_cells_z_) % num_cells_z_;
          } else {
            if (nx >= num_cells_x_ || nx < 0 || ny >= num_cells_y_ || ny < 0 ||
                nz >= num_cells_z_ || nz < 0)
              continue;
          }
          neighbor_cells[n_cells] = nx + ny * num_cells_x_ + nz * num_cells_xy_;
          n_cells++;
        }
    for (int ni = 0; ni < n_cells; ++ni) {
      int nId = neighbor_cells[ni];
      int start = cell_offsets_[nId], end = cell_offsets_[nId + 1];

      for (int k = start; k < end; ++k) {
        int j = flat_cells_[k];
        if (j <= index)
          continue;

        double dx = xn[j] - xi, dy = yn[j] - yi, dz = zn[j] - zi;

        if (pbc) {
          dx -= std::round(dx / box_size_x_) * box_size_x_;
          dy -= std::round(dy / box_size_y_) * box_size_y_;
          dz -= std::round(dz / box_size_z_) * box_size_z_;
        }

        double r2 = dx * dx + dy * dy + dz * dz;
        if (r2 < cutoff_sqr_) {
          neighbors.push_back(j);
          count++;
        }
      }
    }

    return count;
  }

  inline const std::string getData() const {
    std::ostringstream oss;
    oss << "CellList Data:";
    oss << "\n\tCutoff: " << cutoff_;
    oss << "\n\tBox Size: (" << box_size_x_ << ", " << box_size_y_ << ", "
        << box_size_z_ << ")";
    oss << "\n\tCell Size: (" << cell_size_x_ << ", " << cell_size_y_ << ", "
        << cell_size_z_ << ")";
    oss << "\n\tTotal Cells: " << totalCells_;
    oss << "\n\tNumber of Cells: (" << num_cells_x_ << ", " << num_cells_y_
        << ", " << num_cells_z_ << ")";

    return oss.str();
  }

  inline int getNumCells() const { return totalCells_; };

  inline int getCellIndex(const double &pos_x, const double &pos_y,
                          const double &pos_z) const {
    int ix =
        std::min(static_cast<int>(pos_x * inv_cell_size_x_), num_cells_x_ - 1);
    int iy =
        std::min(static_cast<int>(pos_y * inv_cell_size_y_), num_cells_y_ - 1);
    int iz =
        std::min(static_cast<int>(pos_z * inv_cell_size_z_), num_cells_z_ - 1);
    return ix + iy * num_cells_x_ + iz * num_cells_xy_;
  }

private:
  double cutoff_{0.0};
  double cutoff_sqr_{0.0};

  double box_size_x_{0.0}, box_size_y_{0.0}, box_size_z_{0.0};
  double half_box_size_x_{0.0}, half_box_size_y_{0.0}, half_box_size_z_{0.0};
  double inv_box_size_x_{0.0}, inv_box_size_y_{0.0}, inv_box_size_z_{0.0};

  double cell_size_x_{0.0}, cell_size_y_{0.0}, cell_size_z_{0.0};
  double inv_cell_size_x_{0.0}, inv_cell_size_y_{0.0}, inv_cell_size_z_{0.0};

  int num_cells_x_{0}, num_cells_y_{0}, num_cells_z_{0}, num_cells_xy_{0};

  int totalCells_{0};

  std::vector<size_t> flat_cells_;
  std::vector<size_t> cell_offsets_; // размер или offset начала
  std::vector<size_t>
      cell_counts_; // сколько в каждой ячейке (временно во время распределения)

  CellList(const CellList &) = delete;
  CellList &operator=(const CellList &) = delete;
};

#endif
