#ifndef MDALGORITHMS_H
#define MDALGORITHMS_H

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <memory>
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
#include "helpers/Timer.h"
#include "macroparams/Macroparams.h"
#include "output/OutputManager.h"
#include "potentials/Potential.h"
#include "thermostats/Thermostat.h"

using json = nlohmann::json;

/**
 * @brief Класс для работы с алгоритмами молекулярной динамики.
 * @details Класс предоставляет методы для работы с алгоритмами молекулярной
 * динамики.
 */
class MDAlgorithms {
private:
  /**
   * @brief Таймер для измерения времени выполнения.
   */
  Timer timer_{1};

  /**
   * @brief Настройки.
   */
  Settings &settings_;

  /**
   * @brief Пул потоков.
   */
  ThreadPool &threadPool_;

  /**
   * @brief Менеджер бекапов.
   */
  BackupManager &backupManager_;

  /**
   * @brief Менеджер вывода.
   */
  OutputManager &outputManager_;

  /**
   * @brief Менеджер ансамблей.
   */
  std::unique_ptr<EnsembleManager> ensemble_manager_;

  /**
   * @brief Термостат.
   */
  std::unique_ptr<Thermostat> thermostat_;

  /**
   * @brief Баростат.
   */
  std::unique_ptr<Barostat> barostat_;

  /**
   * @brief Потенциал
   */
  std::unique_ptr<Potential> potential_;
  /**
   * @brief Алгоритмы расчета сил.
   */
  ForceAlgorithm force_;

  /**
   * @brief Алгоритмы расчета координат.
   */
  CoordinateAlgorithm coords_;

  /**
   * @brief Алгоритмы расчета скоростей.
   */
  VelocityAlgorithm vels_;

  /**
   * @brief Макропараметры.
   */
  Macroparams macroparams_;

public:
  /**
   * @brief Расчет координат.
   * @param sys_ - система частиц.
   * @details Расчет координат частиц.
   */
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

  /**
   * @brief Расчет скоростей.
   * @param sys_ - система частиц.
   * @details Расчет скоростей частиц.
   */
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

  /**
   * @brief Расчет сил.
   * @param sys_ - система частиц.
   * @details Расчет сил частиц.
   */
  void CalculateForces(System &sys_) {
    int pn = sys_.particleNumber();
    std::vector<Particle> &particles = sys_.particles();

    // Создаем список ячеек
    CellList cellList(potential_->getRcut(), sys_.dimensions().sizes());
    std::cout << cellList.getData() << std::endl;
    cellList.build(particles);

    const int blockSize = cellList.getNumCells() /
                          settings_.threads(); // или любая удобная величина
    std::vector<std::future<void>> futures;

    // Предрасчеты для потенциала EAM
    if (potential_->getPotentialType() == PotentialType::EAM) {
      EAM* eam_potential = dynamic_cast<EAM*>(potential_.get());
      for (int start = 0; start < pn; start += blockSize) {
        int end = std::min(start + blockSize, pn);
        // enqueue задачу на вычисление для блока [start, end)
        futures.push_back(threadPool_.enqueue([this,&eam_potential, &pn, &particles, &cellList,
                                               start, end]() {
          for (int i = start; i < end; i++) {
            particles[i].setElectronDensity(0.0);
            particles[i].setPairPotential(0.0);
            std::vector<int> neighbors = cellList.getNeighbors(particles, i);
            for (int j : neighbors) {
              if (i == j)
                continue;
              force_.preComputeEAM(*eam_potential,particles[i], particles[j],settings_.hasPbc());
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
    if(potential_->getPotentialType() == PotentialType::EAM){
      EAM* eam_potential = dynamic_cast<EAM*>(potential_.get());
      for (int start = 0; start < pn; start += blockSize) {
        int end = std::min(start + blockSize, pn);
        // enqueue задачу на вычисление для блока [start, end)
        futures.push_back(
            threadPool_.enqueue([this,&eam_potential, &pn, &particles, &cellList, start, end]() {
              for (int i = start; i < end; i++) {
                ForceCalcValues result_interaction_i;
                ForceCalcValues interaction_ij;
                std::vector<int> neighbors = cellList.getNeighbors(particles, i);
                for (int j : neighbors) {
                  if (i == j)
                    continue;
                  interaction_ij = force_.compute(*eam_potential,particles[i], particles[j],settings_.hasPbc());
                  result_interaction_i.interaction_count +=
                      interaction_ij.interaction_count;
                  // При ЕАМ возвращается e_pot = 0, так как вычисляется позже
                  // сразу на всю частицу без суммы частей
                  result_interaction_i.e_pot += interaction_ij.e_pot;
                  result_interaction_i.force += interaction_ij.force;
                  result_interaction_i.virials += interaction_ij.virials;
                }
                  // Потенциальная энергия в ЕАМ считается 1 раз сразу для всей
                  // частицы
                result_interaction_i.e_pot =
                      eam_potential->getU(particles[i].electron_density(),particles[i].pairPotential());
                particles[i].applyForceInteraction(result_interaction_i);
              }
            }));
        }
    }
    // Ждём завершения всех
    for (auto &f : futures) {
      f.get();
    }

    futures.clear();
    if(potential_->getPotentialType() == PotentialType::LJ){
      LJ* lj_potential = dynamic_cast<LJ*>(potential_.get());
      for (int start = 0; start < pn; start += blockSize) {
        int end = std::min(start + blockSize, pn);
        // enqueue задачу на вычисление для блока [start, end)
        futures.push_back(
            threadPool_.enqueue([this,&lj_potential, &pn, &particles, &cellList, start, end]() {
              for (int i = start; i < end; i++) {
                ForceCalcValues result_interaction_i;
                ForceCalcValues interaction_ij;
                std::vector<int> neighbors = cellList.getNeighbors(particles, i);
                for (int j : neighbors) {
                  if (i == j)
                    continue;
                  interaction_ij = force_.compute(*lj_potential,particles[i], particles[j],settings_.hasPbc());
                  result_interaction_i.interaction_count +=
                      interaction_ij.interaction_count;
                  result_interaction_i.e_pot += interaction_ij.e_pot;
                  result_interaction_i.force += interaction_ij.force;
                  result_interaction_i.virials += interaction_ij.virials;
                }
                particles[i].applyForceInteraction(result_interaction_i);
              }
            }));
        }
    }
    // Ждём завершения всех
    for (auto &f : futures) {
      f.get();
    }
  }

  /**
   * @brief Нулевой шаг.
   * @param sys_ - система частиц.
   */
  void initialStep(System &sys_) {
    std::cout << "Создание бекапа" << std::endl;
    backupManager_.createBackup(sys_);
    // Расчет сил
    std::cout << "Расчет сил" << std::endl;
    CalculateForces(sys_);
    // Обновление центра масс
    std::cout << "Расчет центра масс" << std::endl;
    sys_.updateVCM();
    // Обновление энергий
    std::cout << "Обновление энергий" << std::endl;
    sys_.updateEnergy();
    // Расчет макропараметров (Температура, давление)
    std::cout << "Расчет температуры" << std::endl;
    sys_.setTemperature(macroparams_.getTemperature(sys_));
    std::cout << "Расчет давления" << std::endl;
    sys_.setPressure(macroparams_.getPressure(sys_));

    if (ensemble_manager_ && ensemble_manager_->enabled()) {
      std::cout << "Расчет транспортных коэффициентов" << std::endl;
      ensemble_manager_->accumulateGreenCubo(sys_);
    }
    // Вывод в файл
    std::cout << "Вывод в файл" << std::endl;
    outputManager_.writeSystemProperties(sys_);
    outputManager_.writeStepData(sys_);
  }

  /**
   * @brief Расчет шага.
   * @param sys_ - система частиц.
   * @details Расчет шага.
   */
  bool advanceStep(System &sys_) {
    sys_.advanceStep();
    // Термостат Ланжевена
    if (thermostat_ && thermostat_->getThermostatType() ==
                           ThermostatType::LANGEVIN) {
      std::cout << "Applying Langevin thermostat" << std::endl;
      thermostat_->applyTemperatureControl(sys_);
    }

    // Расчет первой половины скоростей
    CalculateVelocities(sys_);

    // Расчет новых координат
    CalculateCoordinates(sys_);

    // Баростат Берендсена
    if (barostat_ &&
        barostat_->getBarostatType() == Barostat::BarostatType::BERENDSEN) {
      std::cout << "Applying Berendsen thermostat" << std::endl;
      barostat_->applyPressureControl(sys_);
    }

    // Расчет сил
    CalculateForces(sys_);

    // Расчет второй половины скоростей
    CalculateVelocities(sys_);

    // Термостат Берендсена
    if (thermostat_ && thermostat_->getThermostatType() ==
                           ThermostatType::BERENDSEN) {
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

  /**
   * @brief Конструктор.
   * @param settings - настройки.
   * @param sys - система частиц.
   * @param backupManager - менеджер бекапов.
   * @param outputManager - менеджер вывода.
   * @param threadPool - пул потоков.
   * @param ensemble_manager - менеджер ансамблей.
   * @param potential - потенциал.
   * @param thermostat - термостат.
   * @param barostat - баростат.
   * @param macroparams - макропараметры.
   */
  MDAlgorithms(Settings &settings, System &sys, BackupManager &backupManager,
               OutputManager &outputManager, ThreadPool &threadPool,
               std::unique_ptr<EnsembleManager> &&ensemble_manager,
               std::unique_ptr<Potential> &&potential,
               std::unique_ptr<Thermostat> &&thermostat,
               std::unique_ptr<Barostat> &&barostat, Macroparams &macroparams)
      : settings_(settings), threadPool_(threadPool),
        ensemble_manager_(std::move(ensemble_manager)),
        backupManager_(backupManager), outputManager_(outputManager),
        thermostat_(std::move(thermostat)), barostat_(std::move(barostat)),
        potential_(std::move(potential)),
        force_(settings_, sys), coords_(settings_),
        vels_(settings_), macroparams_(macroparams) {};

  MDAlgorithms() = delete;
  ~MDAlgorithms() = default;

  // Запрещаем копирование
  MDAlgorithms(const MDAlgorithms &) = delete;
  MDAlgorithms &operator=(const MDAlgorithms &) = delete;
};
#endif
