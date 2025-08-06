#ifndef CELL_LIST_H
#define CELL_LIST_H

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

  inline int maxNeighborsCount() {
    return *std::max_element(cell_counts_.begin(), cell_counts_.end());
  }
  inline int minNeighborsCount() {
    return *std::min_element(cell_counts_.begin(), cell_counts_.end());
  }
  inline int avgNeighborsCount() {
    double avg = 0;
    for (int i = 0; i < cell_counts_.size(); i++) {
      avg += cell_counts_[i];
    }
    return avg / cell_counts_.size();
  }

  int getNeighbors(std::vector<int> &neighbors, Particles &particles,
                   const int index, bool pbc) const {
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
          if (r2 > cutoff_sqr_ * 2) {
            printf("Suspicious neighbor i=%d j=%d r2=%.3f cutoff^2=%.3f\n",
                   index, j, r2, cutoff_sqr_);
          }
          neighbors[count] = j;
          count++;
        }
      }
    }

    return count;
    /*
    double rVec_x,rVec_y,rVec_z;
    if(pbc)
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dz = -1; dz <= 1; ++dz) {
          nx = (cx + dx + num_cells_x_) % num_cells_x_;
          ny = (cy + dy + num_cells_y_) % num_cells_y_;
          nz = (cz + dz + num_cells_z_) % num_cells_z_;
          nId = nx + ny * num_cells_x_ + nz * num_cells_xy_;

          // NEW with distance check aka Verlet List
          start = cell_offsets_[nId], end = cell_offsets_[nId + 1];

          //#pragma omp parallel for
          for (size_t k = start; k < end; ++k) {
            j = flat_cells_[k];

            rVec_x = xn[j] - xi;
            rVec_y = yn[j] - yi;
            rVec_z = zn[j] - zi;
            // Отражение частицы
            rVec_x -= (rVec_x > half_box_size_x_) * box_size_x_;
            rVec_x += (rVec_x < -half_box_size_x_) * box_size_x_;

            rVec_y -= (rVec_y > half_box_size_y_) * box_size_y_;
            rVec_y += (rVec_y < -half_box_size_y_) * box_size_y_;

            rVec_z -= (rVec_z > half_box_size_z_) * box_size_z_;
            rVec_z += (rVec_z < -half_box_size_z_) * box_size_z_;

            double length_sqr = rVec_x*rVec_x + rVec_y*rVec_y + rVec_z*rVec_z;
            if ((length_sqr < cutoff_sqr_) && (j != index)) {
              neighbors[count++] = j;
            }
          }
        }
      }
    }
    else
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dz = -1; dz <= 1; ++dz) {
          //Если не ПГУ то не смотрим отраженные ячейки
          int inside =
            ((unsigned)(cx + dx) < (unsigned)num_cells_x_) &
            ((unsigned)(cy + dy) < (unsigned)num_cells_y_) &
            ((unsigned)(cz + dz) < (unsigned)num_cells_z_);
          if (!inside) continue;

          nx = (cx + dx + num_cells_x_) % num_cells_x_;
          ny = (cy + dy + num_cells_y_) % num_cells_y_;
          nz = (cz + dz + num_cells_z_) % num_cells_z_;
          nId = nx + ny * num_cells_x_ + nz * num_cells_xy_;

          // NEW with distance check aka Verlet List
          start = cell_offsets_[nId], end = cell_offsets_[nId + 1];

          //#pragma omp parallel for
          for (size_t k = start; k < end; ++k) {
            j = flat_cells_[k];

            rVec_x = xn[j] - xi;
            rVec_y = yn[j] - yi;
            rVec_z = zn[j] - zi;

            double length_sqr = rVec_x*rVec_x + rVec_y*rVec_y + rVec_z*rVec_z;

            if ((length_sqr < cutoff_sqr_) && (j != index)) {
              neighbors[count] = j;
              count++;
            }
          }
        }
      }
    };
    return count;
    */
  }

  void buildAllNeighbors(std::vector<std::vector<int>> &neighbors,
                         const Particles &particles, bool pbc) {
    neighbors.assign(particles.size(), {}); // очистка

    const double *__restrict__ xn = particles.coord_x_.data();
    const double *__restrict__ yn = particles.coord_y_.data();
    const double *__restrict__ zn = particles.coord_z_.data();

#pragma omp parallel for schedule(dynamic)
    for (int cell = 0; cell < totalCells_; ++cell) {
      int cz = cell / num_cells_xy_;
      int cy = (cell - cz * num_cells_xy_) / num_cells_x_;
      int cx = cell - cy * num_cells_x_ - cz * num_cells_xy_;

      // Собираем список соседних ячеек (включая себя)
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
              if (nx >= num_cells_x_ || nx < 0 || ny >= num_cells_y_ ||
                  ny < 0 || nz >= num_cells_z_ || nz < 0)
                continue;
            }
            neighbor_cells[n_cells++] =
                nx + ny * num_cells_x_ + nz * num_cells_xy_;
          }

      // Все частицы текущей ячейки
      int start_i = cell_offsets_[cell];
      int end_i = cell_offsets_[cell + 1];

      for (int ni = 0; ni < n_cells; ++ni) {
        int nId = neighbor_cells[ni];
        int start_j = cell_offsets_[nId];
        int end_j = cell_offsets_[nId + 1];

        for (int k = start_i; k < end_i; ++k) {
          int i = flat_cells_[k];
          double xi = xn[i], yi = yn[i], zi = zn[i];

          for (int l = start_j; l < end_j; ++l) {
            int j = flat_cells_[l];
            if (j == i)
              continue;

            double dx = xn[j] - xi;
            double dy = yn[j] - yi;
            double dz = zn[j] - zi;

            if (pbc) {
              dx -= std::round(dx / box_size_x_) * box_size_x_;
              dy -= std::round(dy / box_size_y_) * box_size_y_;
              dz -= std::round(dz / box_size_z_) * box_size_z_;
            }

            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 < cutoff_sqr_) {
              neighbors[i].push_back(j);
            }
          }
        }
      }
    }
  }

  void buildAllNeighborsFlat(std::vector<int> &neighbors_flat_index,
                             std::vector<double> &neighbors_flat_distance,
                             std::vector<int> &neighbors_offset,
                             const Particles &particles, bool pbc) {
    const int N = particles.size();
    neighbors_offset.assign(N + 1, 0);

    const double *__restrict__ xn = particles.coord_x_.data();
    const double *__restrict__ yn = particles.coord_y_.data();
    const double *__restrict__ zn = particles.coord_z_.data();

    std::vector<int> local_counts(N, 0);

// -------- 1-й проход: считаем количество соседей --------
#pragma omp parallel
    {
      std::vector<int> private_counts(N, 0);

      int cz, cy, cx, nx, ny, nz;
      int neighbor_cells[27], n_cells;
      int start_i, end_i, start_j, end_j;
      int nId, i, j;
      double xi, yi, zi, dx, dy, dz, r2;
#pragma omp for schedule(dynamic)
      for (int cell = 0; cell < totalCells_; ++cell) {
        cz = cell / num_cells_xy_;
        cy = (cell - cz * num_cells_xy_) / num_cells_x_;
        cx = cell - cy * num_cells_x_ - cz * num_cells_xy_;

        n_cells = 0;
        for (int dx = -1; dx <= 1; ++dx)
          for (int dy = -1; dy <= 1; ++dy)
            for (int dz = -1; dz <= 1; ++dz) {
              nx = cx + dx, ny = cy + dy, nz = cz + dz;
              if (pbc) {
                nx = (nx + num_cells_x_) % num_cells_x_;
                ny = (ny + num_cells_y_) % num_cells_y_;
                nz = (nz + num_cells_z_) % num_cells_z_;
              } else {
                if (nx >= num_cells_x_ || nx < 0 || ny >= num_cells_y_ ||
                    ny < 0 || nz >= num_cells_z_ || nz < 0)
                  continue;
              }
              neighbor_cells[n_cells++] =
                  nx + ny * num_cells_x_ + nz * num_cells_xy_;
            }

        start_i = cell_offsets_[cell];
        end_i = cell_offsets_[cell + 1];

        for (int ni = 0; ni < n_cells; ++ni) {
          nId = neighbor_cells[ni];
          start_j = cell_offsets_[nId];
          end_j = cell_offsets_[nId + 1];

          for (int k = start_i; k < end_i; ++k) {
            i = flat_cells_[k];
            xi = xn[i], yi = yn[i], zi = zn[i];

            for (int l = start_j; l < end_j; ++l) {
              j = flat_cells_[l];
              if (j == i)
                continue;

              dx = xn[j] - xi;
              dy = yn[j] - yi;
              dz = zn[j] - zi;

              if (pbc) {
                dx -= std::round(dx / box_size_x_) * box_size_x_;
                dy -= std::round(dy / box_size_y_) * box_size_y_;
                dz -= std::round(dz / box_size_z_) * box_size_z_;
              }

              r2 = dx * dx + dy * dy + dz * dz;
              if (r2 < cutoff_sqr_) {
                private_counts[i]++;
              }
            }
          }
        }
      }

// Сливаем локальные счетчики
#pragma omp critical
      {
        for (int i = 0; i < N; ++i)
          local_counts[i] += private_counts[i];
      }
    }

    // Префиксная сумма
    for (int i = 0; i < N; ++i)
      neighbors_offset[i + 1] = neighbors_offset[i] + local_counts[i];

    neighbors_flat_index.resize(neighbors_offset[N]);
    neighbors_flat_distance.resize(neighbors_offset[N]);

// -------- 2-й проход: заполняем массив соседей --------
#pragma omp parallel
    {
      std::vector<int> write_pos = neighbors_offset; // локальная копия

      int cz, cy, cx, nx, ny, nz;
      int neighbor_cells[27], n_cells;
      int start_i, end_i, start_j, end_j, pos;
      int nId, i, j;
      double xi, yi, zi, dx, dy, dz, r2, r;
#pragma omp for schedule(dynamic)
      for (int cell = 0; cell < totalCells_; ++cell) {
        cz = cell / num_cells_xy_;
        cy = (cell - cz * num_cells_xy_) / num_cells_x_;
        cx = cell - cy * num_cells_x_ - cz * num_cells_xy_;

        n_cells = 0;
        for (int dx = -1; dx <= 1; ++dx)
          for (int dy = -1; dy <= 1; ++dy)
            for (int dz = -1; dz <= 1; ++dz) {
              nx = cx + dx, ny = cy + dy, nz = cz + dz;
              if (pbc) {
                nx = (nx + num_cells_x_) % num_cells_x_;
                ny = (ny + num_cells_y_) % num_cells_y_;
                nz = (nz + num_cells_z_) % num_cells_z_;
              } else {
                if (nx >= num_cells_x_ || nx < 0 || ny >= num_cells_y_ ||
                    ny < 0 || nz >= num_cells_z_ || nz < 0)
                  continue;
              }
              neighbor_cells[n_cells++] =
                  nx + ny * num_cells_x_ + nz * num_cells_xy_;
            }

        start_i = cell_offsets_[cell];
        end_i = cell_offsets_[cell + 1];

        for (int ni = 0; ni < n_cells; ++ni) {
          nId = neighbor_cells[ni];
          start_j = cell_offsets_[nId];
          end_j = cell_offsets_[nId + 1];

          for (int k = start_i; k < end_i; ++k) {
            i = flat_cells_[k];
            xi = xn[i], yi = yn[i], zi = zn[i];

            for (int l = start_j; l < end_j; ++l) {
              j = flat_cells_[l];
              if (j == i)
                continue;

              dx = xn[j] - xi;
              dy = yn[j] - yi;
              dz = zn[j] - zi;

              if (pbc) {
                dx -= std::round(dx / box_size_x_) * box_size_x_;
                dy -= std::round(dy / box_size_y_) * box_size_y_;
                dz -= std::round(dz / box_size_z_) * box_size_z_;
              }

              r2 = dx * dx + dy * dy + dz * dz;
              if (r2 < cutoff_sqr_) {
                pos = write_pos[i]++;
                r = std::sqrt(r2);
                neighbors_flat_index[pos] = j;
                neighbors_flat_distance[pos] = r;
              }
            }
          }
        }
      }
    }
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

  int getCellIndex(const double &pos_x, const double &pos_y,
                   const double &pos_z) const {
    int ix = std::min(int(pos_x * inv_cell_size_x_), num_cells_x_ - 1);
    int iy = std::min(int(pos_y * inv_cell_size_y_), num_cells_y_ - 1);
    int iz = std::min(int(pos_z * inv_cell_size_z_), num_cells_z_ - 1);

    return ix + iy * num_cells_x_ + iz * num_cells_xy_;
  }

  CellList(const CellList &) = delete;
  CellList &operator=(const CellList &) = delete;
};

#endif
