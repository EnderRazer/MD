#ifndef MDALGORITHMS_H
#define MDALGORITHMS_H

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <vector>
// #include <sleef.h>

#include "backups/BackupManager.h"
#include "barostats/Barostat.h"
#include "classes/CellList.h"
#include "classes/Dimensions.h"
#include "classes/EnsembleManager.h"
#include "classes/Particles.h"
#include "classes/ThreadPool.h"
#include "classes/Timer.h"
#include "macroparams/Macroparams.h"
#include "output/OutputManager.h"
#include "potentials/EAM.h"
#include "potentials/LJ.h"
#include "potentials/Potential.h"
#include "thermostats/Thermostat.h"

using json = nlohmann::json;

class MDAlgorithms {
private:
  Timer timer_{1};
  Settings &settings_;
  ThreadPool &threadPool_;

  BackupManager &backupManager_;
  OutputManager &outputManager_;
  CellList &cell_list_;

  std::unique_ptr<EnsembleManager> ensemble_manager_;
  std::unique_ptr<Thermostat> thermostat_;
  std::unique_ptr<Barostat> barostat_;
  std::unique_ptr<Potential> potential_;

  Macroparams macroparams_;

public:
  void CalculateCoordinates(System &sys_) {
    Particles &particles = sys_.particles();
    Dimensions &dim = sys_.dimensions();
    int pn = particles.size();
    double dt = settings_.dt();
    int pbc = settings_.hasPbc();
    int tn = settings_.threads();

    const double lx = dim.sizeX();
    const double ly = dim.sizeY();
    const double lz = dim.sizeZ();

    double *__restrict__ x_ = particles.coord_x_.data();
    double *__restrict__ y_ = particles.coord_y_.data();
    double *__restrict__ z_ = particles.coord_z_.data();

    double *__restrict__ vx_ = particles.velocity_x_.data();
    double *__restrict__ vy_ = particles.velocity_y_.data();
    double *__restrict__ vz_ = particles.velocity_z_.data();

#pragma omp parallel for simd schedule(static)                                 \
    aligned(x_, y_, z_, vx_, vy_, vz_ : 64)
    for (int i = 0; i < pn; ++i) {
      x_[i] += vx_[i] * dt;
      y_[i] += vy_[i] * dt;
      z_[i] += vz_[i] * dt;

      x_[i] -= pbc * lx * std::floor(x_[i] / lx);
      y_[i] -= pbc * ly * std::floor(y_[i] / ly);
      z_[i] -= pbc * lz * std::floor(z_[i] / lz);
    }
  }

  void CalculateVelocities(System &sys_) {
    Particles &particles = sys_.particles();
    int pn = particles.size();

    double mt = settings_.dt() / (2 * settings_.mass());

    double *__restrict__ vx_ = particles.velocity_x_.data();
    double *__restrict__ vy_ = particles.velocity_y_.data();
    double *__restrict__ vz_ = particles.velocity_z_.data();

    double *__restrict__ fx_ = particles.force_x_.data();
    double *__restrict__ fy_ = particles.force_y_.data();
    double *__restrict__ fz_ = particles.force_z_.data();

#pragma omp parallel for simd schedule(static)
    for (int i = 0; i < pn; ++i) {
      vx_[i] += fx_[i] * mt;
      vy_[i] += fy_[i] * mt;
      vz_[i] += fz_[i] * mt;
    }
  }

  void CalculateEAM(System &sys_, EAM &pot_) {
    Dimensions &dim = sys_.dimensions();
    Particles &particles = sys_.particles();
    const int pn = particles.size();

    int tn = settings_.threads();
    int pbc = settings_.hasPbc();

    double lx = dim.sizeX(), ly = dim.sizeY(), lz = dim.sizeZ();
    double half_lx = dim.halfSizeX(), half_ly = dim.halfSizeY(),
           half_lz = dim.halfSizeZ();

    double *__restrict__ x_ = particles.coord_x_.data();
    double *__restrict__ y_ = particles.coord_y_.data();
    double *__restrict__ z_ = particles.coord_z_.data();

    double *__restrict__ fx_ = particles.force_x_.data();
    double *__restrict__ fy_ = particles.force_y_.data();
    double *__restrict__ fz_ = particles.force_z_.data();

    double *__restrict__ vxx_ = particles.virial_xx_.data();
    double *__restrict__ vxy_ = particles.virial_xy_.data();
    double *__restrict__ vxz_ = particles.virial_xz_.data();

    double *__restrict__ vyx_ = particles.virial_yx_.data();
    double *__restrict__ vyy_ = particles.virial_yy_.data();
    double *__restrict__ vyz_ = particles.virial_yz_.data();

    double *__restrict__ vzx_ = particles.virial_zx_.data();
    double *__restrict__ vzy_ = particles.virial_zy_.data();
    double *__restrict__ vzz_ = particles.virial_zz_.data();

    double *__restrict__ ed_ = particles.electron_density_.data();
    double *__restrict__ pp_ = particles.pair_potential_.data();
    double *__restrict__ ce_ = particles.ce_.data();
    double *__restrict__ d_ce_ = particles.d_ce_.data();

    double *__restrict__ epot_ = particles.e_pot_.data();

    Timer timer{1};

    // Создаем список ячеек
    timer.start();
    std::vector<int> n_flat_index;
    std::vector<double> n_flat_distance;
    std::vector<int> n_offsets;
    cell_list_.rebuild(pot_.getRcut(), dim.sizeX(), dim.sizeY(), dim.sizeZ());
    cell_list_.distribute_particles(particles);
    cell_list_.buildAllNeighborsFlat(n_flat_index, n_flat_distance, n_offsets,
                                     particles, pbc);
    timer.stop();
    std::cout << "Cell list build elapsed time: " << timer.elapsed() << " ms"
              << std::endl;

    timer.start();
// Предрасчеты для потенциала EAM
#pragma omp parallel num_threads(tn)
    {
      int ni, ns;
      int start, end;

      double dx, dy, dz;
      double r2, r;

      double r_over_Re, t1, t2, t3, exp1, exp2, pow1, pow2;
      double mu, rho_f;

#pragma omp for schedule(dynamic)
      for (int i = 0; i < pn; i++) {
        mu = 0.0, rho_f = 0.0;

        start = n_offsets[i];
        end = n_offsets[i + 1];
        ns = end - start;

        // #pragma omp simd aligned(x_, y_, z_ : 64) reduction(+ : mu, rho_f)
        for (int j = start; j < end; j++) {
          r = n_flat_distance[j];
          // Расчеты ЕАМ
          r_over_Re = r * pot_.inv_r_e_;
          t1 = r_over_Re - 1;
          t2 = r_over_Re - pot_.k_;
          t3 = r_over_Re - pot_.lambda_;

          exp1 = exp(-pot_.alpha_ * t1);
          exp2 = exp(-pot_.beta_ * t1);

          pow1 = std::pow(t2, pot_.m_);
          pow2 = std::pow(t3, pot_.n_);

          rho_f += (pot_.f_e_ * exp2) / (1 + pow2);
          mu += (pot_.a_ * exp1 / (1 + pow1)) - (pot_.b_ * exp2 / (1 + pow2));
        }
        ed_[i] = rho_f;
        pp_[i] = mu;
      }
    }
    timer.stop();
    std::cout << "Precalc EAM time elapsed = " << timer.elapsed() << " ms"
              << std::endl;

    timer.start();
    // #pragma omp parallel num_threads(tn)
    {
      double om_val, dom_val;

      double rho_over_RHO_n, rho_pow1, rho_pow2, rho_pow3;
      double rho_over_RHO_s, pow_term, log_term, pow_term_eta_minus1;
      double rho_over_RHO_e;
      // #pragma omp for schedule(dynamic)
      for (int i = 0; i < pn; i++) {
        om_val = 0.0;
        dom_val = 0.0;

        if (ed_[i] < pot_.rho_n_) {
          // ---- om1 и d_om1 ----
          rho_over_RHO_n = ed_[i] / pot_.rho_n_ - 1.0;

          // Horner для многочлена 3-й степени
          rho_pow1 = rho_over_RHO_n;
          rho_pow2 = rho_pow1 * rho_over_RHO_n;
          rho_pow3 = rho_pow2 * rho_over_RHO_n;

          // om1
          om_val = pot_.om_n_[0] + pot_.om_n_[1] * rho_pow1 +
                   pot_.om_n_[2] * rho_pow2 + pot_.om_n_[3] * rho_pow3;

          // d_om1
          dom_val =
              (1.0 * pot_.om_n_[1] * 1.0 + 2.0 * pot_.om_n_[2] * rho_pow1 +
               3.0 * pot_.om_n_[3] * rho_pow2) /
              pot_.rho_n_;

        } else if (ed_[i] >= pot_.rho_0_) {
          // ---- om3 и d_om3 ----
          rho_over_RHO_s = ed_[i] / pot_.rho_s_;

          // pow и log
          pow_term = std::pow(rho_over_RHO_s, pot_.eta_);
          log_term = std::log(pow_term); // = eta_ * log(rho/rho_s)

          // om3
          om_val = pot_.om_e_ * (1.0 - log_term) * pow_term;

          // d_om3
          // my_pow_n(x, eta_-1) = pow_term / x
          pow_term_eta_minus1 = pow_term / rho_over_RHO_s;
          dom_val = -pot_.om_e_ * pot_.eta_ * pow_term_eta_minus1 * log_term /
                    pot_.rho_s_;

        } else {
          // ---- om2 и d_om2 ----
          rho_over_RHO_e = ed_[i] / pot_.rho_e_ - 1.0;

          rho_pow1 = rho_over_RHO_e;
          rho_pow2 = rho_pow1 * rho_over_RHO_e;
          rho_pow3 = rho_pow2 * rho_over_RHO_e;

          // om2
          om_val = pot_.om_[0] + pot_.om_[1] * rho_pow1 +
                   pot_.om_[2] * rho_pow2 + pot_.om_[3] * rho_pow3;

          // d_om2
          dom_val = (1.0 * pot_.om_[1] * 1.0 + 2.0 * pot_.om_[2] * rho_pow1 +
                     3.0 * pot_.om_[3] * rho_pow2) /
                    pot_.rho_e_;
        }

        epot_[i] = pot_.energy_unit_ * (om_val + 0.5 * pp_[i]);
        ce_[i] = om_val;
        d_ce_[i] = dom_val;
      }
    }
    timer.stop();
    std::cout << "Middle EAM time elapsed = " << timer.elapsed() << " ms"
              << std::endl;

    // Основной блок расчетов
    timer.start();
#pragma omp parallel num_threads(tn)
    {
      int ns, ni, start, end;

      double ix, iy, iz;
      double idom, jdom;
      double dx, dy, dz;
      double r2;
      double r, inv_r;

      double fx, fy, fz;
      double vxx, vxy, vxz;
      double vyx, vyy, vyz;
      double vzx, vzy, vzz;

      double FU, Fx, Fy, Fz;

      double r_over_Re, t1, t2, t3;
      double exp1, exp2, pow1, pow1m1, pow2, pow2m1;
      double plus_pow1, inv_plus_pow1, inv_plus_pow1_sqr, plus_pow2,
          inv_plus_pow2, inv_plus_pow2_sqr;
      double var11, var12, var1, var21, var22, var2;
      double var1_rho, var2_rho;
      double d_mu, d_rho_f;

#pragma omp for schedule(dynamic)
      for (int i = 0; i < pn; i++) {
        // Обнуление сил и вириалов
        fx = 0.0, fy = 0.0, fz = 0.0;
        vxx = 0.0, vxy = 0.0, vxz = 0.0;
        vyx = 0.0, vyy = 0.0, vyz = 0.0;
        vzx = 0.0, vzy = 0.0, vzz = 0.0;

        ix = x_[i], iy = y_[i], iz = z_[i], idom = d_ce_[i];

        start = n_offsets[i];
        end = n_offsets[i + 1];
        ns = end - start;

        //#pragma omp simd reduction(+ : fx, fy, fz, vxx, vxy, vxz, vyx, vyy, vyz, vzx,  \
                               vzy, vzz)
        for (int j = start; j < end; j++) {
          ni = n_flat_index[j];
          jdom = d_ce_[ni];
          // Векторное расстояние
          dx = ix - x_[ni];
          dy = iy - y_[ni];
          dz = iz - z_[ni];
          // Отражение частицы
          dx -= pbc * (dx > half_lx) * lx;
          dx += pbc * (dx < -half_lx) * lx;

          dy -= pbc * (dy > half_ly) * ly;
          dy += pbc * (dy < -half_ly) * ly;

          dz -= pbc * (dz > half_lz) * lz;
          dz += pbc * (dz < -half_lz) * lz;

          // Расстояние между частицами
          // r2 = dx*dx + dy*dy + dz*dz;
          // r  = sqrt(r2);
          r = n_flat_distance[j];
          inv_r = 1 / r;
          // Сила потенциала

          r_over_Re = r * pot_.inv_r_e_;
          t1 = r_over_Re - 1.0;
          t2 = r_over_Re - pot_.k_;
          t3 = r_over_Re - pot_.lambda_;

          // Общие экспоненты
          exp1 = exp(-pot_.alpha_ * t1);
          exp2 = exp(-pot_.beta_ * t1);

          // Общие степени
          pow1 = std::pow(t2, pot_.m_); // (r/Re - k)^m
          plus_pow1 = (1.0 + pow1);
          inv_plus_pow1 = 1.0 / plus_pow1;
          inv_plus_pow1_sqr = inv_plus_pow1 * inv_plus_pow1;
          pow1m1 = std::pow(t2, pot_.m_ - 1); // (r/Re - k)^(m-1)
          pow2 = std::pow(t3, pot_.n_);       // (r/Re - λ)^n
          plus_pow2 = (1.0 + pow2);
          inv_plus_pow2 = 1.0 / plus_pow2;
          inv_plus_pow2_sqr = inv_plus_pow2 * inv_plus_pow2;
          pow2m1 = std::pow(t3, pot_.n_ - 1); // (r/Re - λ)^(n-1)

          // ------------------- d_mu -------------------
          var11 = (-pot_.a_ * pot_.alpha_ * pot_.inv_r_e_) * exp1 * plus_pow1;
          var12 = (pot_.m_ * pot_.inv_r_e_) * pow1m1 * pot_.a_ * exp1;
          var1 = (var11 - var12) * inv_plus_pow1_sqr;

          var21 = (pot_.n_ * pot_.inv_r_e_) * pow2m1 * pot_.b_ * exp2;
          var22 = (-pot_.b_ * pot_.beta_ * pot_.inv_r_e_) * exp2 * plus_pow2;
          var2 = (var21 - var22) * inv_plus_pow2_sqr;

          // ------------------- d_rho_f -------------------
          var1_rho = (-pot_.beta_ * pot_.inv_r_e_) * exp2 * plus_pow2;
          var2_rho = (pot_.n_ * pot_.inv_r_e_) * pow2m1 * exp2;
          // var3_rho = 1.0 + pow2;

          d_mu = var1 + var2;
          d_rho_f = pot_.f_e_ * (var1_rho - var2_rho) * inv_plus_pow2_sqr;

          FU = -pot_.energy_unit_ * (((idom + jdom) * d_rho_f) + d_mu);

          Fx = dx * FU * inv_r;
          Fy = dy * FU * inv_r;
          Fz = dz * FU * inv_r;

          fx += Fx;
          fy += Fy;
          fz += Fz;
          vxx += dx * Fx;
          vxy += dx * Fy;
          vxz += dx * Fz;
          vyx += dy * Fx;
          vyy += dy * Fy;
          vyz += dy * Fz;
          vzx += dz * Fx;
          vzy += dz * Fy;
          vzz += dz * Fz;
        }
        fx_[i] = fx;
        fy_[i] = fy;
        fz_[i] = fz;
        vxx_[i] = vxx;
        vxy_[i] = vxy;
        vxz_[i] = vxz;
        vyx_[i] = vyx;
        vyy_[i] = vyy;
        vyz_[i] = vyz;
        vzx_[i] = vzx;
        vzy_[i] = vzy;
        vzz_[i] = vzz;
      }
    }
    timer.stop();
    std::cout << "Main EAM time elapsed " << timer.elapsed() << " ms"
              << std::endl;

    // Основной блок расчетов
    /*
    if (pot_.getPotentialType() == Potential::PotentialType::LJ) {
      #pragma omp parallel num_threads(tn)
      {
        std::vector<int> neighbors;
        neighbors.resize(cell_list_.maxNeighborsCount());

        double e_pot;
        double fx,fy,fz;
        double vxx, vxy, vxz;
        double vyx, vyy, vyz;
        double vzx, vzy, vzz;
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < pn; i++) {
          //Обнуление сил и вириалов
          e_pot = 0.0;
          fx = fy = fz = 0.0;
          vxx = vxy = vxz = vyx = vyy = vyz = vzx = vzy = vzz = 0.0;

          int ns = cell_list_.getNeighbors(neighbors, particles, i,
    settings_.hasPbc());
          //std::cout << "Neighbors count for particle "<<i<<" :
    "<<neighbors.size() << std::endl;
          //#pragma omp simd reduction(+:e_pot, fx, fy, fz, vxx, vxy, vxz, vyx,
    vyy, vyz, vzx, vzy, vzz) for (int j = 0; j < ns;j++) { int ni =
    neighbors[j];

            double dx = x_[i] - x_[ni];
            double dy = y_[i] - y_[ni];
            double dz = z_[i] - z_[ni];

            // Отражение частицы
            dx -= pbc * lx * std::nearbyint(dx * inv_lx);
            dy -= pbc * ly * std::nearbyint(dy * inv_ly);
            dz -= pbc * lz * std::nearbyint(dz * inv_lz);

            double r2 = dx*dx + dy*dy + dz*dz;
            double r = std::sqrt(r2);
            PotentialResult res = pot_.getAll(r);
            // U = potential_->getU(lengthSqr);
            // FU = potential_->getFU(lengthSqr);
            double U = res.u / 2;
            double FU = res.fu;

            double Fx = dx * FU / r;
            double Fy = dy * FU / r;
            double Fz = dz * FU / r;

            e_pot += U;
            fx += Fx;
            fy += Fy;
            fz += Fz;
            vxx += dx * Fx; vxy += dx * Fy; vxz += dx * Fz;
            vyx += dy * Fx; vyy += dy * Fy; vyz += dy * Fz;
            vzx += dz * Fx; vzy += dz * Fy; vzz += dz * Fz;
          }

          epot_[i] = e_pot;
          fx_[i] = fx; fy_[i] = fy; fz_[i] = fz;
          vxx_[i] = vxx; vxy_[i] = vxy; vxz_[i] = vxz;
          vyx_[i] = vyx; vyy_[i] = vyy; vyz_[i] = vyz;
          vzx_[i] = vzx; vzy_[i] = vzy; vzz_[i] = vzz;
        }
      }
    }
    */
  }

  void initialStep(System &sys_) {
    backupManager_.createBackup(sys_);
    // Расчет сил
    timer_.start();
    switch (potential_->getPotentialType()) {
    case Potential::PotentialType::EAM:
      CalculateEAM(sys_, *static_cast<EAM *>(potential_.get()));
      break;
    case Potential::PotentialType::LJ:
      // CalculateForces(sys_, *static_cast<LJ*>(potential_.get()));
      break;
    default:
      throw std::invalid_argument("Unknown potential");
      break;
    }
    timer_.stop();
    std::cout << "CalculateForce elapsed time: " << timer_.elapsed() << " ms"
              << std::endl;
    // CalculateForces(sys_);
    //  Обновление центра масс
    timer_.start();
    sys_.updateVCM();
    timer_.stop();
    std::cout << "updateVCM elapsed time: " << timer_.elapsed() << " ms"
              << std::endl;
    // Обновление энергий
    timer_.start();
    sys_.updateEnergy();
    timer_.stop();
    std::cout << "updateEnergy elapsed time: " << timer_.elapsed() << " ms"
              << std::endl;
    // Обновление импульса
    timer_.start();
    sys_.updatePulse();
    timer_.stop();
    std::cout << "updatePulse elapsed time: " << timer_.elapsed() << " ms"
              << std::endl;
    // Расчет макропараметров (Температура, давление)
    timer_.start();
    sys_.setTemperature(macroparams_.getTemperature(sys_));
    sys_.setPressure(macroparams_.getPressure(sys_));
    if (ensemble_manager_ && ensemble_manager_->enabled()) {
      std::cout << "Расчет транспортных коэффициентов" << std::endl;
      ensemble_manager_->accumulateGreenCubo(sys_);
    }
    timer_.stop();
    std::cout << "calculateMacroparams elapsed time: " << timer_.elapsed()
              << " ms" << std::endl;
    // Вывод в файл
    outputManager_.writeSystemProperties(sys_);
    outputManager_.writeStepData(sys_);
  }

  bool advanceStep(System &sys_) {
    sys_.advanceStep();
    // Термостат Ланжевена
    if (thermostat_ && thermostat_->getThermostatType() ==
                           Thermostat::ThermostatType::LANGEVIN) {
      timer_.start();
      std::cout << "Applying Langevin thermostat" << std::endl;
      thermostat_->applyTemperatureControl(sys_);
      timer_.stop();
      std::cout << "applyTemperatureControl elapsed time: " << timer_.elapsed()
                << " ms" << std::endl;
    }

    timer_.start();
    // Расчет первой половины скоростей
    CalculateVelocities(sys_);
    timer_.stop();
    std::cout << "CalculateVelocities elapsed time: " << timer_.elapsed()
              << " ms" << std::endl;

    timer_.start();
    // Расчет новых координат
    CalculateCoordinates(sys_);
    timer_.stop();
    std::cout << "CalculateCoordinates elapsed time: " << timer_.elapsed()
              << " ms" << std::endl;

    // Баростат Берендсена
    if (barostat_ &&
        barostat_->getBarostatType() == Barostat::BarostatType::BERENDSEN) {
      timer_.start();
      std::cout << "Applying Berendsen barostat" << std::endl;
      barostat_->applyPressureControl(sys_);
      timer_.stop();
      std::cout << "applyPressureControl elapsed time: " << timer_.elapsed()
                << " ms" << std::endl;
    }

    // Расчет сил
    timer_.start();
    switch (potential_->getPotentialType()) {
    case Potential::PotentialType::EAM:
      CalculateEAM(sys_, *static_cast<EAM *>(potential_.get()));
      break;
    case Potential::PotentialType::LJ:
      // CalculateForces(sys_, *static_cast<LJ*>(potential_.get()));
      break;
    default:
      throw std::invalid_argument("Unknown potential");
      break;
    }
    timer_.stop();
    std::cout << "CalculateForces elapsed time: " << timer_.elapsed() << " ms"
              << std::endl;

    // Расчет второй половины скоростей
    timer_.start();
    CalculateVelocities(sys_);
    timer_.stop();
    std::cout << "CalculateVelocities elapsed time: " << timer_.elapsed()
              << " ms" << std::endl;

    // Термостат Берендсена
    if (thermostat_ && thermostat_->getThermostatType() ==
                           Thermostat::ThermostatType::BERENDSEN) {
      timer_.start();
      std::cout << "Applying Berendsen thermostat" << std::endl;
      thermostat_->applyTemperatureControl(sys_);
      timer_.stop();
      std::cout << "applyTemperatureControl elapsed time: " << timer_.elapsed()
                << " ms" << std::endl;
    }
    // Делаем бекап после расчета основных параметров
    if (sys_.currentStep() % backupManager_.frequency() == 0) {
      backupManager_.createBackup(sys_);
    }
    // Обновление центра масс
    timer_.start();
    sys_.updateVCM();
    timer_.stop();
    std::cout << "updateVCM elapsed time: " << timer_.elapsed() << " ms"
              << std::endl;
    // Обновление энергий
    timer_.start();
    sys_.updateEnergy();
    timer_.stop();
    std::cout << "updateEnergy elapsed time: " << timer_.elapsed() << " ms"
              << std::endl;
    // Обновление импульса
    timer_.start();
    sys_.updatePulse();
    timer_.stop();
    std::cout << "updatePulse elapsed time: " << timer_.elapsed() << " ms"
              << std::endl;
    // Расчет макропараметров (Температура, давление)
    timer_.start();
    sys_.setTemperature(macroparams_.getTemperature(sys_));
    sys_.setPressure(macroparams_.getPressure(sys_));
    if (ensemble_manager_ && ensemble_manager_->enabled()) {
      std::cout << "Расчет транспортных коэффициентов" << std::endl;
      ensemble_manager_->accumulateGreenCubo(sys_);
    }
    timer_.stop();
    std::cout << "calculateMacroparams elapsed time: " << timer_.elapsed()
              << " ms" << std::endl;

    if (sys_.currentStep() % outputManager_.frequency() == 0) {
      // Вывод в файл
      outputManager_.writeSystemProperties(sys_);
      outputManager_.writeAvgSystemProperties(sys_);
      outputManager_.writeStepData(sys_);
    }
    return false;
  }

  MDAlgorithms() = delete;
  MDAlgorithms(Settings &settings, System &sys, BackupManager &backupManager,
               OutputManager &outputManager, ThreadPool &threadPool,
               CellList &cellList,
               std::unique_ptr<EnsembleManager> &&ensemble_manager,
               std::unique_ptr<Thermostat> &&thermostat,
               std::unique_ptr<Barostat> &&barostat,
               std::unique_ptr<Potential> &&potential, json macroparams_config)
      : settings_(settings), threadPool_(threadPool),
        ensemble_manager_(std::move(ensemble_manager)),
        backupManager_(backupManager), outputManager_(outputManager),
        cell_list_(cellList), thermostat_(std::move(thermostat)),
        barostat_(std::move(barostat)), potential_(std::move(potential)),
        macroparams_(macroparams_config, settings) {};

  ~MDAlgorithms() = default;

  // Запрещаем копирование
  MDAlgorithms(const MDAlgorithms &) = delete;
  MDAlgorithms &operator=(const MDAlgorithms &) = delete;
};
#endif
