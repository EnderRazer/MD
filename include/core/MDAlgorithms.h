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
  CellList      &cell_list_;

  std::unique_ptr<EnsembleManager> ensemble_manager_;
  std::unique_ptr<Thermostat> thermostat_;
  std::unique_ptr<Barostat> barostat_;
  std::unique_ptr<Potential> potential_;

  //ForceAlgorithm force_;
  //CoordinateAlgorithm coords_;
  //VelocityAlgorithm vels_;
  Macroparams macroparams_;

public:
  void CalculateCoordinates(System &sys_) {
    Timer timer{1};
    timer.start();
    Particles &particles = sys_.particles();
    Dimensions &dim = sys_.dimensions();
    double dt = settings_.dt();
    int pbc = settings_.hasPbc();
    
    const double lx = dim.sizeX();
    const double ly = dim.sizeY();
    const double lz = dim.sizeZ();

    #pragma clang loop vectorize(enable)
    for (size_t i = 0; i < particles.size(); ++i) {
      particles.coord_x_[i] += particles.velocity_x_[i] * dt;
      particles.coord_y_[i] += particles.velocity_y_[i] * dt;
      particles.coord_z_[i] += particles.velocity_z_[i] * dt;

      particles.coord_x_[i] -= pbc * lx * std::floor(particles.coord_x_[i] / lx);
      particles.coord_y_[i] -= pbc * ly * std::floor(particles.coord_y_[i] / ly);
      particles.coord_z_[i] -= pbc * lz * std::floor(particles.coord_z_[i] / lz);
    }

    //coords_.compute(particles);
    /*
    for(int i=0;i<pn;i++){
      coords_.compute(particles, i);
    }
    */
    /*
    const int blockSize =
        pn / settings_.threads(); // или любая удобная величина
    std::vector<std::future<void>> futures;
    for (int start = 0; start < pn; start += blockSize) {
      int end = std::min(start + blockSize, pn);
      // enqueue задачу на вычисление для блока [start, end)
      futures.push_back(
          threadPool_.enqueue([this, &particles, &dim, start, end]() {
            for (int i = start; i < end; i++) {
              coords_.compute(particles, i);
            }
          }));
    }
    // Ждём завершения всех
    for (auto &f : futures) {
      f.get();
    }
    */
    timer.stop();
    std::cout<<"Coord elapsed time: "<<timer.elapsed()<<" ms"<<std::endl;
  }

  void CalculateVelocities(System &sys_) {
    Timer timer{1};
    timer.start();
    Particles &particles = sys_.particles();
    double mt = settings_.dt() / (2 * settings_.mass());
    //int pn = particles.size();
    //vels_.compute(particles);
    #pragma clang loop vectorize(enable)
    for (size_t i = 0; i < particles.size(); ++i) {
      particles.velocity_x_[i] += particles.force_x_[i]*mt;
      particles.velocity_y_[i] += particles.force_y_[i]*mt;
      particles.velocity_z_[i] += particles.force_z_[i]*mt;
    }
    
    /*
    const int blockSize =
        pn / settings_.threads(); // или любая удобная величина
    std::vector<std::future<void>> futures;

    for (int start = 0; start < pn; start += blockSize) {
      int end = std::min(start + blockSize, pn);
      // enqueue задачу на вычисление для блока [start, end)
      futures.push_back(threadPool_.enqueue([this, &particles, start, end]() {
        std::cout << "Thread: " << start << " - " << end << std::endl;
        for (int i = start; i < end; i++) {
          vels_.compute(particles, i);
        }
      }));
    }
    // Ждём завершения всех
    for (auto &f : futures) {
      f.get();
    }
    */
    timer.stop();
    std::cout<<"Velocity elapsed time: "<<timer.elapsed()<<" ms"<<std::endl;
  }
  template<typename Potential>
  void CalculateForces(System &sys_, Potential &pot_) {
    Dimensions &dim = sys_.dimensions();
    Particles &particles = sys_.particles();
    int pn = particles.size();

    double* __restrict__ x_ = particles.coord_x_.data();
    double* __restrict__ y_ = particles.coord_y_.data();
    double* __restrict__ z_ = particles.coord_z_.data();

    double* __restrict__ fx_ = particles.force_x_.data();
    double* __restrict__ fy_ = particles.force_y_.data();
    double* __restrict__ fz_ = particles.force_z_.data();

    double* __restrict__ vxx_ = particles.virial_xx_.data();
    double* __restrict__ vxy_ = particles.virial_xy_.data();
    double* __restrict__ vxz_ = particles.virial_xz_.data();

    double* __restrict__ vyx_ = particles.virial_yx_.data();
    double* __restrict__ vyy_ = particles.virial_yy_.data();
    double* __restrict__ vyz_ = particles.virial_yz_.data();

    double* __restrict__ vzx_ = particles.virial_zx_.data();
    double* __restrict__ vzy_ = particles.virial_zy_.data();
    double* __restrict__ vzz_ = particles.virial_zz_.data();

    double* __restrict__ ed_ = particles.electron_density_.data();
    double* __restrict__ pp_ = particles.pair_potential_.data();

    double* __restrict__ epot_ = particles.e_pot_.data();

    // Создаем список ячеек
    Timer timer{1};
    timer.start();
    cell_list_.rebuild(pot_.getRcut(), dim.sizeX(), dim.sizeY(), dim.sizeZ());
    //CellList cellList(force_.getCutOff(), dim.sizeX(), dim.sizeY(), dim.sizeZ());
    //std::cout << cell_list_.getData() << std::endl;
    cell_list_.distribute_particles(particles);
    timer.stop();
    std::cout<<"Cell list build elapsed time: "<<timer.elapsed()<<" ms"<<std::endl;
    
    
    std::vector<int> neighbors;
    neighbors.resize(pn);
    int pbc = settings_.hasPbc();
    int ns = 0;

    double dx = 0.0, dy = 0.0 , dz = 0.0;
    double fx = 0.0, fy = 0.0, fz = 0.0;
    double vxx = 0.0, vxy = 0.0,vxz = 0.0;
    double vyx = 0.0, vyy = 0.0,vyz = 0.0;
    double vzx = 0.0, vzy = 0.0,vzz = 0.0;

    double lx = dim.sizeX(), ly = dim.sizeY(), lz = dim.sizeZ();
    double inv_lx = 1/dim.sizeX(), inv_ly = 1/dim.sizeY(), inv_lz = 1/dim.sizeZ();
    double length,lengthSqr;
    double U = 0,FU = 0;
    double mu,rho_f;

    if (pot_.getPotentialType() == Potential::PotentialType::EAM) {
      timer.start();
      // Предрасчеты для потенциала EAM
      for (int i = 0; i < pn; i++) {
        mu = 0.0, rho_f = 0.0;
        ns = cell_list_.getNeighbors(neighbors, particles, i, pbc);
        // Векторизация по соседям
        #pragma omp simd reduction(+:mu, rho_f)
        for (int j = 0;j < ns;j++) {
          int ni = neighbors[j];
          //Векторное расстояние
          double dx = x_[i] - x_[ni];
          double dy = y_[i] - y_[ni];
          double dz = z_[i] - z_[ni];
          // Отражение частицы
          dx -= pbc * lx * (long long)(dx * inv_lx + (dx >= 0 ? 0.5 : -0.5));
          dy -= pbc * ly * (long long)(dy * inv_ly + (dy >= 0 ? 0.5 : -0.5));
          dz -= pbc * lz * (long long)(dz * inv_lz + (dz >= 0 ? 0.5 : -0.5));
          //Расстояние между частицами
          double r2 = dx*dx + dy*dy + dz*dz;
          double r  = sqrt(r2);
          //Расчеты ЕАМ
          mu += pot_.getPairPart(length);
          rho_f += pot_.getDensityPart(length);
        }
        ed_[i] = rho_f;
        pp_[i] = mu;
      }
      timer.stop();
      std::cout<<"Precalc EAM time elapsed " << timer.elapsed() << " ms" << std::endl;
      // Основной блок расчетов
      timer.start();
      for (int i = 0; i < pn; i++) {
        //Обнуление сил и вириалов
        fx = 0.0,fy = 0.0,fz = 0.0;
        vxx = 0.0, vxy = 0.0, vxz = 0.0;
        vyx = 0.0, vyy = 0.0, vyz = 0.0;
        vzx = 0.0, vzy = 0.0, vzz = 0.0;

        ns = cell_list_.getNeighbors(neighbors, particles,i, pbc);

        #pragma omp simd reduction(+:fx, fy, fz, vxx, vxy, vxz, vyx, vyy, vyz, vzx, vzy, vzz)
        for (int j = 0; j < ns;j++) {
          //Сосед
          int ni = neighbors[j];
          //Векторное расстояние
          double dx = x_[i] - x_[ni];
          double dy = y_[i] - y_[ni];
          double dz = z_[i] - z_[ni];
          // Отражение частицы
          dx -= pbc * lx * (long long)(dx * inv_lx + (dx >= 0 ? 0.5 : -0.5));
          dy -= pbc * ly * (long long)(dy * inv_ly + (dy >= 0 ? 0.5 : -0.5));
          dz -= pbc * lz * (long long)(dz * inv_lz + (dz >= 0 ? 0.5 : -0.5));
          //Расстояние между частицами
          double r2 = dx*dx + dy*dy + dz*dz;
          double r  = sqrt(r2);
          //Сила потенциала
          double FU = pot_.getFU(ed_[i], ed_[ni], r);    
           
          double Fx = dx * FU / r;
          double Fy = dy * FU / r;
          double Fz = dz * FU / r;

          fx += Fx; fy += Fy; fz += Fz;
          vxx += dx*Fx; vxy += dx*Fy; vxz += dx*Fz;
          vyx += dy*Fx; vyy += dy*Fy; vyz += dy*Fz;
          vzx += dz*Fx; vzy += dz*Fy; vzz += dz*Fz;
        }
        epot_[i] = pot_.getU(ed_[i], pp_[i]);
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
      timer.stop();
      std::cout<<"Main EAM time elapsed " << timer.elapsed() << " ms" << std::endl;
    }
    // Основной блок расчетов
    if (pot_.getPotentialType() == Potential::PotentialType::LJ) {
      for (int i = 0; i < pn; i++) {
        //Обнуление сил и вириалов
        particles.e_pot_[i] = 0.0;
        particles.force_x_[i] = 0.0;
        particles.force_y_[i] = 0.0;
        particles.force_z_[i] = 0.0;
        particles.virial_xx_[i] = 0.0;
        particles.virial_xy_[i] = 0.0;
        particles.virial_xz_[i] = 0.0;
        particles.virial_yx_[i] = 0.0;
        particles.virial_yy_[i] = 0.0;
        particles.virial_yz_[i] = 0.0;
        particles.virial_zx_[i] = 0.0;
        particles.virial_zy_[i] = 0.0;
        particles.virial_zz_[i] = 0.0;

        ns = cell_list_.getNeighbors(neighbors, particles, i, settings_.hasPbc());
        //std::cout << "Neighbors count for particle "<<i<<" : "<<neighbors.size() << std::endl;
        for (int j = 0; j < ns;j++) {
          int ni = neighbors[j];
          
          dx = particles.coord_x_[i] - particles.coord_x_[ni];
          dy = particles.coord_y_[i] - particles.coord_y_[ni];
          dz = particles.coord_z_[i] - particles.coord_z_[ni];

          // Отражение частицы
          dx -= pbc * lx * std::nearbyint(dx * inv_lx);
          dy -= pbc * ly * std::nearbyint(dy * inv_ly);
          dz -= pbc * lz * std::nearbyint(dz * inv_lz);

          lengthSqr = dx*dx + dy*dy + dz*dz;
          if (lengthSqr > pot_.getSqrRcut())
            continue;
          length = std::sqrt(lengthSqr);
          PotentialResult res = pot_.getAll(lengthSqr);
          // U = potential_->getU(lengthSqr);
          // FU = potential_->getFU(lengthSqr);
          U = res.u / 2;
          FU = res.fu;
          
          particles.e_pot_[i] += U;
          particles.force_x_[i] += dx * FU / length;
          particles.force_y_[i] += dy * FU / length;
          particles.force_z_[i] += dz * FU / length;
          particles.virial_xx_[i] += dx * particles.force_x_[i];
          particles.virial_xy_[i] += dx * particles.force_y_[i];
          particles.virial_xz_[i] += dx * particles.force_z_[i];
          particles.virial_yx_[i] += dy * particles.force_x_[i];
          particles.virial_yy_[i] += dy * particles.force_y_[i];
          particles.virial_yz_[i] += dy * particles.force_z_[i];
          particles.virial_zx_[i] += dz * particles.force_x_[i];
          particles.virial_zy_[i] += dz * particles.force_y_[i];
          particles.virial_zz_[i] += dz * particles.force_z_[i];
        }
      }
    }
  }

  void initialStep(System &sys_) {
    backupManager_.createBackup(sys_);
    // Расчет сил
    switch (potential_->getPotentialType()) {
    case Potential::PotentialType::EAM:
      CalculateForces(sys_, *static_cast<EAM*>(potential_.get()));
      break;
    case Potential::PotentialType::LJ:
      CalculateForces(sys_, *static_cast<LJ*>(potential_.get()));
      break;
    default:
      throw std::invalid_argument("Unknown potential");
      break;
    }
    //CalculateForces(sys_);
    // Обновление центра масс
    sys_.updateVCM();
    // Обновление энергий
    sys_.updateEnergy();
    // Обновление импульса
    sys_.updatePulse();
    // Расчет макропараметров (Температура, давление)
    sys_.setTemperature(macroparams_.getTemperature(sys_));
    sys_.setPressure(macroparams_.getPressure(sys_));

    if (ensemble_manager_ && ensemble_manager_->enabled()) {
      std::cout << "Расчет транспортных коэффициентов" << std::endl;
      ensemble_manager_->accumulateGreenCubo(sys_);
    }
    // Вывод в файл
    outputManager_.writeSystemProperties(sys_);
    outputManager_.writeStepData(sys_);
  }
  bool advanceStep(System &sys_) {
    sys_.advanceStep();
    // Термостат Ланжевена
    if (thermostat_ && thermostat_->getThermostatType() ==
                           Thermostat::ThermostatType::LANGEVIN) {
      //std::cout << "Applying Langevin thermostat" << std::endl;
      thermostat_->applyTemperatureControl(sys_);
    }

    // Расчет первой половины скоростей
    CalculateVelocities(sys_);

    // Расчет новых координат
    CalculateCoordinates(sys_);

    // Баростат Берендсена
    if (barostat_ &&
        barostat_->getBarostatType() == Barostat::BarostatType::BERENDSEN) {
      //std::cout << "Applying Berendsen barostat" << std::endl;
      barostat_->applyPressureControl(sys_);
    }

    // Расчет сил
    switch (potential_->getPotentialType()) {
    case Potential::PotentialType::EAM:
      CalculateForces(sys_, *static_cast<EAM*>(potential_.get()));
      break;
    case Potential::PotentialType::LJ:
      CalculateForces(sys_, *static_cast<LJ*>(potential_.get()));
      break;
    default:
      throw std::invalid_argument("Unknown potential");
      break;
    }
    // Расчет сил
    //CalculateForces(sys_);

    // Расчет второй половины скоростей
    CalculateVelocities(sys_);

    // Термостат Берендсена
    if (thermostat_ && thermostat_->getThermostatType() ==
                           Thermostat::ThermostatType::BERENDSEN) {
      std::cout << "Applying Berendsen thermostat" << std::endl;
      thermostat_->applyTemperatureControl(sys_);
    }
    // Делаем бекап после расчета основных параметров
    if (sys_.currentStep() % backupManager_.frequency() == 0) {
      backupManager_.createBackup(sys_);
    }
    // Обновление центра масс
    sys_.updateVCM();
    // Обновление энергий
    sys_.updateEnergy();
    // Обновление импульса
    sys_.updatePulse();
    // Расчет макропараметров (Температура, давление)
    sys_.setTemperature(macroparams_.getTemperature(sys_));
    sys_.setPressure(macroparams_.getPressure(sys_));
    // Транспортные коэффициенты
    if (macroparams_.enabled()) {
      if (ensemble_manager_ && ensemble_manager_->enabled() &&
          !ensemble_manager_->completed()) {
        ensemble_manager_->accumulateGreenCubo(sys_);
      }
      if (ensemble_manager_->completed()) {
        return true;
      }
    }

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
               OutputManager &outputManager, ThreadPool &threadPool, CellList &cellList,
               std::unique_ptr<EnsembleManager> &&ensemble_manager,
               std::unique_ptr<Thermostat> &&thermostat,
               std::unique_ptr<Barostat> &&barostat, 
               std::unique_ptr<Potential>&&potential,
               json macroparams_config)
      : settings_(settings), threadPool_(threadPool),
        ensemble_manager_(std::move(ensemble_manager)),
        backupManager_(backupManager), outputManager_(outputManager),
        cell_list_(cellList),
        thermostat_(std::move(thermostat)), barostat_(std::move(barostat)),
        potential_(std::move(potential)),
        macroparams_(macroparams_config, settings) {};

  ~MDAlgorithms() = default;

  // Запрещаем копирование
  MDAlgorithms(const MDAlgorithms &) = delete;
  MDAlgorithms &operator=(const MDAlgorithms &) = delete;
};
#endif
