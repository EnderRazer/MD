#ifndef MDALGORITHMS_H
#define MDALGORITHMS_H

#include "CoordinateAlgorithm.h"
#include "ForceAlgorithm.h"
#include "VelocityAlgorithm.h"

class MDAlgorithms
{
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
  void CalculateCoordinates(System &sys_)
  {
    int pn = sys_.particleNumber();
    std::vector<Particle> &particles = sys_.particles();
    Dimensions dim = sys_.dimensions();

    const int blockSize =
        pn / settings_.threads(); // или любая удобная величина
    std::vector<std::future<void>> futures;
    for (int start = 0; start < pn; start += blockSize)
    {
      int end = std::min(start + blockSize, pn);
      // enqueue задачу на вычисление для блока [start, end)
      futures.push_back(
          threadPool_.enqueue([this, &particles, &dim, start, end]()
                              {
            for (int i = start; i < end; i++) {
              coords_.compute(particles[i]);
              if (settings_.hasPbc())
                coords_.applyPBC(particles[i], dim);
            } }));
    }
    // Ждём завершения всех
    for (auto &f : futures)
    {
      f.get();
    }
  }
  void CalculateVelocities(System &sys_)
  {
    int pn = sys_.particleNumber();
    std::vector<Particle> &particles = sys_.particles();

    const int blockSize =
        pn / settings_.threads(); // или любая удобная величина
    std::vector<std::future<void>> futures;

    for (int start = 0; start < pn; start += blockSize)
    {
      int end = std::min(start + blockSize, pn);
      // enqueue задачу на вычисление для блока [start, end)
      futures.push_back(threadPool_.enqueue([this, &particles, start, end]()
                                            {
        for (int i = start; i < end; i++) {
          vels_.compute(particles[i]);
        } }));
    }
    // Ждём завершения всех
    for (auto &f : futures)
    {
      f.get();
    }
  }
  void CalculateForces(System &sys_)
  {
    int pn = sys_.particleNumber();
    std::vector<Particle> &particles = sys_.particles();

    const int blockSize =
        pn / settings_.threads(); // или любая удобная величина
    std::vector<std::future<void>> futures;

    for (int start = 0; start < pn; start += blockSize)
    {
      int end = std::min(start + blockSize, pn);
      // enqueue задачу на вычисление для блока [start, end)
      futures.push_back(
          threadPool_.enqueue([this, &pn, &particles, start, end]()
                              {
            for (int i = start; i < end; i++) {
              ForceCalcValues result_interaction_i;
              ForceCalcValues interaction_ij;
              for (int j = 0; j < pn; j++) {
                if (i == j)
                  continue;
                interaction_ij = force_.compute(particles[i], particles[j]);
                result_interaction_i.e_pot += interaction_ij.e_pot;
                result_interaction_i.force += interaction_ij.force;
                result_interaction_i.virials += interaction_ij.virials;
              }
              particles[i].applyForceInteraction(result_interaction_i);
            } }));
    }
    // Ждём завершения всех
    for (auto &f : futures)
    {
      f.get();
    }
  }

  void initialStep(System &sys_)
  {
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

    if (ensemble_manager_ && ensemble_manager_->enabled())
    {
      std::cout << "Расчет транспортных коэффициентов" << std::endl;
      ensemble_manager_->accumulateGreenCubo(sys_);
      outputManager_.writeEnsemblesDetailed(ensemble_manager_->getEnsembles());
    }
    // Вывод в файл
    outputManager_.writeSystemProperties(sys_);
    outputManager_.writeStepData(sys_);
  }
  bool advanceStep(System &sys_)
  {
    sys_.advanceStep();
    // Термостат Ланжевена
    if (thermostat_ && thermostat_->getThermostatType() ==
                           Thermostat::ThermostatType::LANGEVIN)
    {
      std::cout << "Applying Langevin thermostat" << std::endl;
      thermostat_->applyTemperatureControl(sys_);
    }

    // Расчет первой половины скоростей
    CalculateVelocities(sys_);

    // Расчет новых координат
    CalculateCoordinates(sys_);

    // Баростат Берендсена
    if (barostat_ &&
        barostat_->getBarostatType() == Barostat::BarostatType::BERENDSEN)
    {
      std::cout << "Applying Berendsen thermostat" << std::endl;
      barostat_->applyPressureControl(sys_);
    }

    timer_.start();
    // Расчет сил
    CalculateForces(sys_);
    timer_.stop();
    std::cout << "Force calc time: " << timer_.elapsed() << std::endl;
    force_avg_time += timer_.elapsed();

    // Расчет второй половины скоростей
    CalculateVelocities(sys_);

    // Термостат Берендсена
    if (thermostat_ && thermostat_->getThermostatType() ==
                           Thermostat::ThermostatType::BERENDSEN)
    {
      std::cout << "Applying Berendsen thermostat" << std::endl;
      thermostat_->applyTemperatureControl(sys_);
    }
    // Делаем бекап после расчета основных параметров
    if (sys_.currentStep() % backupManager_.frequency() == 0)
    {
      backupManager_.createBackup(sys_);
    }
    // Обновление центра масс
    sys_.updateVCM();
    // Обновление энергий
    sys_.updateEnergy();
    sys_.updateEnergyAvg();
    // Обновление импульса
    sys_.updatePulse();
    // Расчет макропараметров (Температура, давление)
    sys_.setTemperature(macroparams_.getTemperature(sys_));
    sys_.updateTemperatureAvg();
    sys_.setPressure(macroparams_.getPressure(sys_));
    sys_.updatePressureAvg();
    // Транспортные коэффициенты
    if (macroparams_.enabled())
    {
      if (ensemble_manager_ && ensemble_manager_->enabled() && !ensemble_manager_->completed())
      {
        ensemble_manager_->accumulateGreenCubo(sys_);
      }
    }

    if (sys_.currentStep() % outputManager_.frequency() == 0)
    {

      // Вывод в файл
      outputManager_.writeSystemProperties(sys_);
      outputManager_.writeStepData(sys_);
      if (macroparams_.enabled() && ensemble_manager_->enabled())
      {
        // Вывод в файлы ансамблей
        outputManager_.writeEnsemblesDetailed(ensemble_manager_->getEnsembles());

        // Вывод усредненных значений ансамблей
        if (ensemble_manager_->completed())
        {
          outputManager_.writeAvgEnsembleDetailed(ensemble_manager_.get());
          return true;
        }
      }
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
        ensemble_manager_(std::move(ensemble_manager)), backupManager_(backupManager),
        outputManager_(outputManager), thermostat_(std::move(thermostat)),
        barostat_(std::move(barostat)),
        force_(settings_, std::move(potential), sys.dimensions()),
        coords_(settings_), vels_(settings_),
        macroparams_(macroparams_config, settings) {};

  ~MDAlgorithms() = default;

  // Запрещаем копирование
  MDAlgorithms(const MDAlgorithms &) = delete;
  MDAlgorithms &operator=(const MDAlgorithms &) = delete;
};
#endif