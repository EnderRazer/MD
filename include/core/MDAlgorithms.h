#ifndef MDALGORITHMS_H
#define MDALGORITHMS_H

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <vector>

#include "CoordinateAlgorithm.h"
#include "ForceAlgorithm.h"
#include "VelocityAlgorithm.h"
#include "backups/BackupManager.h"
#include "barostats/Barostat.h"
#include "classes/CellList.h"
#include "classes/EnsembleManager.h"
#include "classes/Particle.h"
#include "classes/ThreadPool.h"
#include "classes/Timer.h"
#include "macroparams/Macroparams.h"
#include "output/OutputManager.h"
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

  std::unique_ptr<EnsembleManager> ensemble_manager_;
  std::unique_ptr<Thermostat> thermostat_;
  std::unique_ptr<Barostat> barostat_;

  ForceAlgorithm force_;
  CoordinateAlgorithm coords_;
  VelocityAlgorithm vels_;
  Macroparams macroparams_;

public:
  void CalculateCoordinates(System &sys_) {
    int pn = sys_.particleNumber();
    std::vector<Particle> &particles = sys_.particles();
    Dimensions dim = sys_.dimensions();

    const int blockSize =
        pn / settings_.threads(); // или любая удобная величина
    std::vector<std::future<void>> futures;
    for (int start = 0; start < pn; start += blockSize) {
      int end = std::min(start + blockSize, pn);
      // enqueue задачу на вычисление для блока [start, end)
      futures.push_back(
          threadPool_.enqueue([this, &particles, &dim, start, end]() {
            for (int i = start; i < end; i++) {
              coords_.compute(particles[i]);
              if (settings_.hasPbc())
                coords_.applyPBC(particles[i], dim);
            }
          }));
    }
    // Ждём завершения всех
    for (auto &f : futures) {
      f.get();
    }
  }

  void CalculateVelocities(System &sys_) {
    int pn = sys_.particleNumber();
    std::vector<Particle> &particles = sys_.particles();

    const int blockSize =
        pn / settings_.threads(); // или любая удобная величина
    std::vector<std::future<void>> futures;

    for (int start = 0; start < pn; start += blockSize) {
      int end = std::min(start + blockSize, pn);
      // enqueue задачу на вычисление для блока [start, end)
      futures.push_back(threadPool_.enqueue([this, &particles, start, end]() {
        for (int i = start; i < end; i++) {
          vels_.compute(particles[i]);
        }
      }));
    }
    // Ждём завершения всех
    for (auto &f : futures) {
      f.get();
    }
  }

  void CalculateForces(System &sys_) {
    int pn = sys_.particleNumber();
    std::vector<Particle> &particles = sys_.particles();

    // Создаем список ячеек
    CellList cellList(force_.getCutOff(), sys_.dimensions().sizes());
   //std::cout << cellList.getData() << std::endl;
    cellList.build(particles);

    const int blockSize = cellList.getNumCells() /
                          settings_.threads(); // или любая удобная величина
    std::vector<std::future<void>> futures;

    // Предрасчеты для потенциала EAM
    if (force_.getPotentialType() == Potential::PotentialType::EAM) {
      for (int start = 0; start < pn; start += blockSize) {
        int end = std::min(start + blockSize, pn);
        // enqueue задачу на вычисление для блока [start, end)
        futures.push_back(threadPool_.enqueue([this, &pn, &particles, &cellList,
                                               start, end]() {
          for (int i = start; i < end; i++) {
            particles[i].setElectronDensity(0.0);
            particles[i].setPairPotential(0.0);
            std::vector<int> neighbors = cellList.getNeighbors(particles, i, settings_.hasPbc());
            for (int j : neighbors) {
              if (i == j)
                continue;
              force_.preCompute(particles[i], particles[j]);
            }
          }
        }));
      }
      // Ждём завершения всех
      for (auto &f : futures) {
        f.get();
      }
    }
    // Основной блок расчетов
    futures.clear();
    for (int start = 0; start < pn; start += blockSize) {
      int end = std::min(start + blockSize, pn);
      // enqueue задачу на вычисление для блока [start, end)
      futures.push_back(
          threadPool_.enqueue([this, &pn, &particles, &cellList, start, end]() {
            for (int i = start; i < end; i++) {
              ForceCalcValues result_interaction_i;
              ForceCalcValues interaction_ij;
              std::vector<int> neighbors = cellList.getNeighbors(particles, i, settings_.hasPbc());
              for (int j : neighbors) {
                if (i == j)
                  continue;
                interaction_ij = force_.compute(particles[i], particles[j]);
                result_interaction_i.interaction_count +=
                    interaction_ij.interaction_count;
                // При ЕАМ возвращается e_pot = 0, так как вычисляется позже
                // сразу на всю частицу без суммы частей
                result_interaction_i.e_pot += interaction_ij.e_pot;
                result_interaction_i.force += interaction_ij.force;
                result_interaction_i.virials += interaction_ij.virials;
              }
              if (force_.getPotentialType() == Potential::PotentialType::EAM) {
                // Потенциальная энергия в ЕАМ считается 1 раз сразу для всей
                // частицы
                result_interaction_i.e_pot =
                    force_.PotentialEnergy_EAM(particles[i]);
              }
              particles[i].applyForceInteraction(result_interaction_i);
            }
          }));
    }
    // Ждём завершения всех
    for (auto &f : futures) {
      f.get();
    }
    /*
        double max_Fx = -INFINITY, max_Fy = -INFINITY, max_Fz = -INFINITY;
        double min_Fx = INFINITY, min_Fy = INFINITY, min_Fz = INFINITY;
        double sum_Fx = 0, sum_Fy = 0, sum_Fz = 0;

        double max_VirXX = -INFINITY, max_VirYY = -INFINITY, max_VirZZ =
       -INFINITY; double min_VirXX = INFINITY, min_VirYY = INFINITY, min_VirZZ =
       INFINITY; double sum_VirXX = 0, sum_VirYY = 0, sum_VirZZ = 0; for (int i
       = 0; i < pn; i++) { if (particles[i].getForceX() > max_Fx) max_Fx =
       particles[i].getForceX(); if (particles[i].getForceX() < min_Fx) min_Fx =
       particles[i].getForceX();

          if (particles[i].getForceY() > max_Fy)
            max_Fy = particles[i].getForceY();
          if (particles[i].getForceY() < min_Fy)
            min_Fy = particles[i].getForceY();

          if (particles[i].getForceZ() > max_Fz)
            max_Fz = particles[i].getForceZ();
          if (particles[i].getForceZ() < min_Fz)
            min_Fz = particles[i].getForceZ();

          sum_Fx += particles[i].getForceX();
          sum_Fy += particles[i].getForceY();
          sum_Fz += particles[i].getForceZ();

          if (particles[i].virials().xx() > max_VirXX)
            max_VirXX = particles[i].virials().xx();
          if (particles[i].virials().xx() < min_VirXX)
            min_VirXX = particles[i].virials().xx();

          if (particles[i].virials().yy() > max_VirYY)
            max_VirYY = particles[i].virials().yy();
          if (particles[i].virials().yy() < min_VirYY)
            min_VirYY = particles[i].virials().yy();

          if (particles[i].virials().zz() > max_VirZZ)
            max_VirZZ = particles[i].virials().zz();
          if (particles[i].virials().zz() < min_VirZZ)
            min_VirZZ = particles[i].virials().zz();

          sum_VirXX += particles[i].virials().xx();
          sum_VirYY += particles[i].virials().yy();
          sum_VirZZ += particles[i].virials().zz();
        }

        std::cout << "Force max: \n\tFx_max: " << max_Fx << "\n\tFy_max: " <<
       max_Fy
                  << "\n\tFz_max: " << max_Fz
                  << "\nForce min: \n\tFx_min: " << min_Fx
                  << "\n\tFy_min: " << min_Fy << "\n\tFz_min: " << min_Fz
                  << "\nForce sum: \n\tFx_sum: " << sum_Fx
                  << "\n\tFy_sum: " << sum_Fy << "\n\tFz_sum: " << sum_Fz <<
       "\n";

        std::cout << "Virial max: \n\tVirXX_max: " << max_VirXX
                  << "\n\tVirYY_max: " << max_VirYY
                  << "\n\tVirZZ_max: " << max_VirZZ
                  << "\nVirial min: \n\tVirXX_min: " << min_VirXX
                  << "\n\tVirYY_min: " << min_VirYY
                  << "\n\tVirZZ_min: " << min_VirZZ
                  << "\nVirial sum: \n\tVirXX_sum: " << sum_VirXX
                  << "\n\tVirYY_sum: " << sum_VirYY
                  << "\n\tVirZZ_sum: " << sum_VirZZ << "\n";
                  */
  }

  void initialStep(System &sys_) {
    backupManager_.createBackup(sys_);
    // Расчет сил
    CalculateForces(sys_);
    // Обновление центра масс
    sys_.updateVCM();
    // Обновление энергий
    sys_.updateEnergy();
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
    CalculateForces(sys_);

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
    sys_.updateEnergyAvg();
    // Обновление импульса
    sys_.updatePulse();
    sys_.updatePulseAvg();
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
               OutputManager &outputManager, ThreadPool &threadPool,
               std::unique_ptr<EnsembleManager> &&ensemble_manager,
               std::unique_ptr<Potential> &&potential,
               std::unique_ptr<Thermostat> &&thermostat,
               std::unique_ptr<Barostat> &&barostat, json macroparams_config)
      : settings_(settings), threadPool_(threadPool),
        ensemble_manager_(std::move(ensemble_manager)),
        backupManager_(backupManager), outputManager_(outputManager),
        thermostat_(std::move(thermostat)), barostat_(std::move(barostat)),
        force_(settings_, std::move(potential), sys), coords_(settings_),
        vels_(settings_), macroparams_(macroparams_config, settings) {};

  ~MDAlgorithms() = default;

  // Запрещаем копирование
  MDAlgorithms(const MDAlgorithms &) = delete;
  MDAlgorithms &operator=(const MDAlgorithms &) = delete;
};
#endif
