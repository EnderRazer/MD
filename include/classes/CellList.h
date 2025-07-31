#ifndef CELL_LIST_H
#define CELL_LIST_H

#include <algorithm>
#include <cassert>
#include <iostream>
#include <ostream>
#include <sstream>

#include "classes/Particles.h"

class CellList {
public:
  CellList() = default;
  
  void rebuild(const double &cutoff, const double &boxSize_x,
             const double &boxSize_y, const double &boxSize_z) {
  cutoff_ = cutoff;
  cutoff_sqr_ = cutoff*cutoff;

  box_size_x_ = boxSize_x, box_size_y_ = boxSize_y, box_size_z_ = boxSize_z;

  inv_box_size_x_ = 1/boxSize_x, inv_box_size_y_ = 1/boxSize_y, inv_box_size_z_ = 1/boxSize_z;

  num_cells_x_ = std::max(1, static_cast<int>(box_size_x_ / cutoff_));
  num_cells_y_ = std::max(1, static_cast<int>(box_size_y_ / cutoff_));
  num_cells_z_ = std::max(1, static_cast<int>(box_size_z_ / cutoff_));
  num_cells_xy_ = num_cells_x_*num_cells_y_;

  cell_size_x_ = box_size_x_ / num_cells_x_, cell_size_y_ = box_size_y_ / num_cells_y_, cell_size_z_ = box_size_z_ / num_cells_z_;
  inv_cell_size_x_ = 1/cell_size_x_, inv_cell_size_y_ = 1/cell_size_y_, inv_cell_size_z_ = 1/cell_size_z_;

  totalCells_ = num_cells_x_ * num_cells_y_ * num_cells_z_;
}
  void distribute_particles(const Particles &particles) {
    // 1. Очистка
    flat_cells_.clear();
    cell_counts_.assign(totalCells_, 0);  // временно

    // 2. Подсчет сколько частиц в каждой ячейке
    for (int i = 0; i < particles.size(); ++i) {
      int cellIdx = getCellIndex(particles.coordX(i), particles.coordY(i), particles.coordZ(i));
      ++cell_counts_[cellIdx];
    }
    // 3. Построение offset-таблицы
    cell_offsets_.resize(totalCells_ + 1);  // +1 для удобства границы
    cell_offsets_[0] = 0;
    for (int i = 0; i < totalCells_; ++i) {
      cell_offsets_[i + 1] = cell_offsets_[i] + cell_counts_[i];
    }

    // 4. Заполнение частиц
    flat_cells_.resize(particles.size());
    std::vector<size_t> current_offsets = cell_offsets_;  // текущая позиция вставки
    for (int i = 0; i < particles.size(); ++i) {
      int cellIdx = getCellIndex(particles.coordX(i), particles.coordY(i), particles.coordZ(i));
      flat_cells_[current_offsets[cellIdx]++] = i;
    }
  }

  int maxNeighborsCount() {
    return *std::max_element(cell_counts_.begin(),cell_counts_.end())*27;
  }

  int getNeighbors(std::vector<int>&neighbors, Particles &particles, const int index, bool pbc) const {
    int count=0;
    double xi = particles.coord_x_[index];
    double yi = particles.coord_y_[index];
    double zi = particles.coord_z_[index];

    double* __restrict__ xn = particles.coord_x_.data();
    double* __restrict__ yn = particles.coord_y_.data();
    double* __restrict__ zn = particles.coord_z_.data();
    
    int cellIdx = getCellIndex(xi, yi, zi);
    int cz = cellIdx / num_cells_xy_;
    int cy = (cellIdx - cz * num_cells_xy_) / num_cells_x_;
    int cx = cellIdx - cy * num_cells_x_ - cz * num_cells_xy_;

    int nx,ny,nz,nId;
    int start,end,j;

    double rVec_x,rVec_y,rVec_z;
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dz = -1; dz <= 1; ++dz) {
          //Если не ПГУ то не смотрим отраженные ячейки
          int inside =
            ((unsigned)(cx + dx) < (unsigned)num_cells_x_) &
            ((unsigned)(cy + dy) < (unsigned)num_cells_y_) &
            ((unsigned)(cz + dz) < (unsigned)num_cells_z_);
          if (!pbc && !inside) continue;

          nx = (cx + dx + num_cells_x_) % num_cells_x_;
          ny = (cy + dy + num_cells_y_) % num_cells_y_;
          nz = (cz + dz + num_cells_z_) % num_cells_z_;
          nId = nx + ny * num_cells_x_ + nz * num_cells_xy_;

          // NEW with distance check aka Verlet List
          start = cell_offsets_[nId], end = cell_offsets_[nId + 1];

          #pragma omp parallel for
          for (size_t k = start; k < end; ++k) {
            j = flat_cells_[k];
            
            rVec_x = xn[j] - xi;
            rVec_y = yn[j] - yi;
            rVec_z = zn[j] - zi;
            // Mirrorig vector (PBC)
            rVec_x -= pbc * box_size_x_ * (long long)(rVec_x * inv_box_size_x_ + (rVec_x >= 0 ? 0.5 : -0.5));
            rVec_y -= pbc * box_size_y_ * (long long)(rVec_y * inv_box_size_y_ + (rVec_y >= 0 ? 0.5 : -0.5));
            rVec_z -= pbc * box_size_z_ * (long long)(rVec_z * inv_box_size_z_ + (rVec_z >= 0 ? 0.5 : -0.5));
              
            double length_sqr = rVec_x*rVec_x + rVec_y*rVec_y + rVec_z*rVec_z;
            int mask = (length_sqr < cutoff_sqr_) & (j != index);
            neighbors[count]=j;
            count+=mask;
          }
        }
      }
    }
    return count;
  }

  inline const std::string getData() const {
    std::ostringstream oss;
    oss << "CellList Data:";
    oss << "\n\tCutoff: " << cutoff_;
    oss << "\n\tBox Size: (" << box_size_x_ << ", "<<box_size_y_ << ", "<<box_size_z_<<")";
    oss << "\n\tCell Size: (" << cell_size_x_ << ", "<<cell_size_y_ << ", "<<cell_size_z_<<")";
    oss << "\n\tTotal Cells: " << totalCells_;
    oss << "\n\tNumber of Cells: (" << num_cells_x_ << ", "<<num_cells_y_ << ", "<<num_cells_z_<<")";

    return oss.str();
  }

  inline int getNumCells() const { return totalCells_; };

private:
  double cutoff_{0.0};
  double cutoff_sqr_{0.0};

  double box_size_x_{0.0},box_size_y_{0.0},box_size_z_{0.0};
  double inv_box_size_x_{0.0},inv_box_size_y_{0.0},inv_box_size_z_{0.0};

  double cell_size_x_{0.0},cell_size_y_{0.0},cell_size_z_{0.0};
  double inv_cell_size_x_{0.0},inv_cell_size_y_{0.0},inv_cell_size_z_{0.0};

  int num_cells_x_{0},num_cells_y_{0},num_cells_z_{0}, num_cells_xy_{0};

  int totalCells_{0};

  std::vector<int> flat_cells_;
  std::vector<size_t> cell_offsets_;  // размер или offset начала
  std::vector<size_t> cell_counts_;   // сколько в каждой ячейке (временно во время распределения)

  int getCellIndex(const double &pos_x,const double &pos_y,const double &pos_z) const {
    int ix = std::min(int(pos_x * inv_cell_size_x_), num_cells_x_ - 1);
    int iy = std::min(int(pos_y * inv_cell_size_y_), num_cells_y_ - 1);
    int iz = std::min(int(pos_z * inv_cell_size_z_), num_cells_z_ - 1);

    return ix + iy * num_cells_x_ + iz * num_cells_x_ * num_cells_y_;
  }

  CellList(const CellList &) = delete;
  CellList &operator=(const CellList &) = delete;
};

#endif
