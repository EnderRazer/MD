#ifndef ENSEMBLE_MANAGER_H
#define ENSEMBLE_MANAGER_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "classes/Particles.h"
#include "nlohmann/json.hpp"

#include "classes/IncrementalIntegrator.h" //Интегральная функция
#include "classes/ThreadPool.h"
#include "core/Settings.h"
#include "core/System.h"
//#include "macroparams/Einstein.h"
#include "macroparams/GreenKubo.h"
#include "output/OutputManager.h"

struct Ensemble {
  int id_{};  // Идентификатор ансамбля
  int cur_{}; // Внутренний счетчик ансамбля

  std::string filename_{};      // Имя файла для записи данных ансамбля
  std::ofstream output_file_{}; // Файл для записи данных ансамбля

  std::vector<double> init_velocity_x_{};  // Начальные компоненты скорости
  std::vector<double> init_velocity_y_{};  // (для расчетов коэффициента диффузии)
  std::vector<double> init_velocity_z_{};

  std::vector<double> accum_acfv_{}; // Автокорреляционные функции скорости
  std::vector<double> coef_diff_integral_{}; // Интеграл от АКФ

  // Начальные тензоры давления системы
  double init_pressure_xx_{}, init_pressure_xy_{}, init_pressure_xz_{};
  double init_pressure_yx_{}, init_pressure_yy_{}, init_pressure_yz_{};
  double init_pressure_zx_{}, init_pressure_zy_{}, init_pressure_zz_{}; 

  std::vector<double> accum_acfp_{}; // Автокорреляционные функции тензоров давления
  std::vector<double> coef_visc_integral_{}; // Интеграл от АКФ
};

using json = nlohmann::json;

class EnsembleManager {
private:
  ThreadPool &threadPool_;       // Пул потоков для параллельных вычислений
  OutputManager &outputManager_; // Менеджер вывода данных

  IncrementalIntegrator integral_diff_; // Интеграл от АКФ
  IncrementalIntegrator integral_visc_; // Интеграл от АКФ
  GreenKubo gk; // Методы Грин-Кубо для транспортных коэффициентов

  bool toggle_{false};
  bool finished_{false}; // Флаг завершения работы
  int nsteps_{0};        // Количество шагов

  int ensemble_threads_num_; // Количество параллельных "потоков" ансамблей
  int ensemble_size_;        // Размер ансамбля (кол-во шагов)
  int ensemble_offset_;      // Отступ между ансамблями
  int ensemble_count_;       // Сколько всего полных ансамблей собирать

  int current_ens_count{0};
  std::vector<Ensemble> ensembles_; // Вектор ансамблей
  Ensemble avg_ensemble_;           // Сборник ансамблей

  inline void flushEnsemble(Ensemble &ens) {
    ens.cur_ = 0;
    for (int i = 0; i < ensemble_size_; i++) {
      // Записываем
      avg_ensemble_.accum_acfv_[i] += ens.accum_acfv_[i];
      avg_ensemble_.accum_acfp_[i] += ens.accum_acfp_[i];
      avg_ensemble_.coef_diff_integral_[i] += ens.coef_diff_integral_[i];
      avg_ensemble_.coef_visc_integral_[i] += ens.coef_visc_integral_[i];
      // Очищаем
      ens.accum_acfv_[i] = 0.0;
      ens.accum_acfp_[i] = 0.0;
      ens.coef_diff_integral_[i] = 0.0;
      ens.coef_visc_integral_[i] = 0.0;
    }
    ens.filename_ = "";
    ens.output_file_.close();
  };

  inline void initEnsemble(int thread_finished_count, Ensemble &ens,
                           System &sys) {
    if (current_ens_count > 0)
      ens.id_ += ensembles_.size();
    if (ens.id_ >= ensemble_count_)
      return;
    std::cout << "Initializing ensemble with ID: " << ens.id_ << std::endl;

    // Если не успеваем рассчитать полный ансамбль, то и не считаем, чтобы не
    // забирать ресурсы
    // Или уже рассчитано нужное количество ансамблей
    if (sys.currentStep() + ensemble_size_ > nsteps_) {
      ens.cur_ = -INT32_MAX; // Просто очень большое отрицательное число до
                             // которого не дойдем
      return;
    }
    const int pn = sys.particleNumber();
    const Particles &particles = sys.particles();
    for (int j = 0; j < pn; j++) {
      ens.init_velocity_x_[j] = particles.velocity_x_[j];
      ens.init_velocity_y_[j] = particles.velocity_y_[j];
      ens.init_velocity_z_[j] = particles.velocity_z_[j];
    }
    ens.init_pressure_xx_ = sys.pressureXX();
    ens.init_pressure_xy_ = sys.pressureXY();
    ens.init_pressure_xz_ = sys.pressureXZ();

    ens.init_pressure_yx_ = sys.pressureYX();
    ens.init_pressure_yy_ = sys.pressureYY();
    ens.init_pressure_yz_ = sys.pressureYZ();

    ens.init_pressure_zx_ = sys.pressureZX();
    ens.init_pressure_zy_ = sys.pressureZY();
    ens.init_pressure_zz_ = sys.pressureZZ();

    ens.filename_ = "ensemble_" + std::to_string(ens.id_) + ".csv";
    ens.output_file_ = std::ofstream(
        outputManager_.getEnsembleDir() + "/" + ens.filename_, std::ios::app);
  }

  void writeEnsembleDetailed(Ensemble &ens) const {
    std::vector<std::string> headers = {"Step", "ACFV", "CDiff", "ACFP",
                                        "CVisc"};

    std::vector<double> row = {
        static_cast<double>(ens.cur_), ens.accum_acfv_[ens.cur_],
        ens.coef_diff_integral_[ens.cur_], ens.accum_acfp_[ens.cur_],
        ens.coef_visc_integral_[ens.cur_]};

    bool writeHeaders = (ens.cur_ == 0);

    outputManager_.writeRowToStream(ens.output_file_, headers, row,
                                    writeHeaders);
  }

  void writeAvgEnsembleDetailed(Ensemble &avg_ensemble) const {
    std::vector<std::string> headers = {"Step", "ACFV", "CDiff", "ACFP",
                                        "CVisc"};
    std::vector<std::vector<double>> data;

    // Prepare all data rows
    for (int i = 0; i < ensemble_size_; i++) {
      std::vector<double> row = {
          static_cast<double>(i), avg_ensemble.accum_acfv_[i],
          avg_ensemble.coef_diff_integral_[i], avg_ensemble.accum_acfp_[i],
          avg_ensemble.coef_visc_integral_[i]};
      data.push_back(row);
    }

    // Use output manager to write the data to the output_file_ stream
    outputManager_.writeDataToStream(avg_ensemble.output_file_, headers, data,
                                     true);
  }

public:
  EnsembleManager(json &config, Settings &settings, System &sys,
                  ThreadPool &threadPool, OutputManager &outputManager)
      : integral_diff_(settings.dt()), integral_visc_(settings.dt()),
        threadPool_(threadPool), outputManager_(outputManager) {
    toggle_ = config.value("toggle", false);
    if (toggle_) {
      nsteps_ = settings.steps();
      ensemble_threads_num_ = config.value("threads", 1);
      ensemble_size_ = config.value("size", 0);
      ensemble_offset_ = config.value("offset", 0);
      ensemble_count_ = config.value("quantity", 0);
      ensembles_.resize(ensemble_threads_num_);

      avg_ensemble_.filename_ = "ensembles_average.csv";
      avg_ensemble_.output_file_ = std::ofstream(
          outputManager_.getEnsembleDir() + "/" + avg_ensemble_.filename_);
      avg_ensemble_.accum_acfv_.resize(ensemble_size_, 0.0);
      avg_ensemble_.accum_acfp_.resize(ensemble_size_, 0.0);
      avg_ensemble_.coef_diff_integral_.resize(ensemble_size_, 0.0);
      avg_ensemble_.coef_visc_integral_.resize(ensemble_size_, 0.0);

      int expected_contribution = 0;
      for (int i = 0; i < ensemble_threads_num_; i++) {
        ensembles_[i].id_ = i;
        ensembles_[i].filename_ = "ensemble_" + std::to_string(i) + ".csv";
        ensembles_[i].cur_ = -i * ensemble_offset_;
        ensembles_[i].init_velocity_x_.resize(sys.particleNumber(),0.0);
        ensembles_[i].init_velocity_y_.resize(sys.particleNumber(),0.0);
        ensembles_[i].init_velocity_z_.resize(sys.particleNumber(),0.0);
        ensembles_[i].accum_acfv_.resize(ensemble_size_, 0.0);
        ensembles_[i].accum_acfp_.resize(ensemble_size_, 0.0);
        ensembles_[i].coef_diff_integral_.resize(ensemble_size_, 0.0);
        ensembles_[i].coef_visc_integral_.resize(ensemble_size_, 0.0);

        int ensemble_expected_contribution =
            settings.steps() / (ensemble_size_ + i * ensemble_offset_);
        expected_contribution += ensemble_expected_contribution;
      }
      if (expected_contribution < ensemble_count_) {
        std::ostringstream oss;
        oss << "With current setup accumulated ensembles will be: "
            << std::to_string(expected_contribution) << "/"
            << std::to_string(ensemble_count_);
        std::cerr << oss.str() << std::endl;
        exit(0);
      }
    }
  };
  // Запрещаем копирование
  EnsembleManager(const EnsembleManager &) = delete;
  EnsembleManager &operator=(const EnsembleManager &) = delete;

  ~EnsembleManager() = default;

  inline const bool enabled() const { return toggle_; }
  inline const bool completed() const { return finished_; };

  int getEnsembleThreadsNum() const { return ensemble_threads_num_; }
  int getEnsembleSize() const { return ensemble_size_; }
  int getEnsembleOffset() const { return ensemble_offset_; }
  int getEnsembleCount() const { return ensemble_count_; }
  int getCurrentEnsCount() const { return current_ens_count; }
  const std::vector<Ensemble> &getEnsembles() const { return ensembles_; }
  const std::vector<double> &getEnsemblesAccumAcfv() const {
    return avg_ensemble_.accum_acfv_;
  }
  const std::vector<double> &getEnsemblesAccumAcfp() const {
    return avg_ensemble_.accum_acfp_;
  }
  const std::vector<double> &getEnsemblesCoefDiffIntegral() const {
    return avg_ensemble_.coef_diff_integral_;
  }
  const std::vector<double> &getEnsemblesCoefViscIntegral() const {
    return avg_ensemble_.coef_visc_integral_;
  }
  std::string getData() const {
    std::ostringstream oss;
    oss << "Using ensembles to calculate transport coefficients"
        << "\n\tКоличество параллельных потоков ансамблей: "
        << ensemble_threads_num_
        << "\n\tКоличество шагов в ансамбле: " << ensemble_size_
        << "\n\tОтступ между ансамблями (шагов): " << ensemble_offset_
        << "\n\tКоличество рассчитанных ансамблей: " << ensemble_count_
        << "\n\tEnsembles Size: " << ensembles_.size();

    for (const auto &ensemble : ensembles_) {
      oss << "\n\t\tID Ансамбля: " << ensemble.id_
          << "\n\t\tТекущий шаг: " << ensemble.cur_
          << "\n\t\tInit Velocity Size: " << ensemble.init_velocity_x_.size()
          << "\n\t\tAccum ACFV Size: " << ensemble.accum_acfv_.size()
          << "\n\t\tCoef Diff Integral Size: "
          << ensemble.coef_diff_integral_.size()

          << "\n\t\tAccum ACFP Size: " << ensemble.accum_acfp_.size()
          << "\n\t\tCoef Visc Integral Size: "
          << ensemble.coef_visc_integral_.size() << "\n";
    }
    return oss.str();
  }

  void accumulateGreenCubo(System &sys) {
    if (current_ens_count < ensemble_count_) {
      const int blockSize =
          ensemble_threads_num_ /
          threadPool_.getThreads(); // или любая удобная величина
      std::cout << "Block size: " << blockSize << std::endl;
      std::vector<std::future<int>> futures;

      for (int start = 0; start < ensemble_threads_num_; start += blockSize) {
        int end = std::min(start + blockSize, ensemble_threads_num_);
        // enqueue задачу на вычисление для блока [start, end)
        futures.push_back(threadPool_.enqueue([this, &sys, start, end]() {
          int thread_finished_count = 0;
          if (start > ensembles_.size() || end > ensembles_.size()) {
            std::cerr << "Index overflow: " << start << "/" << end
                      << " >= " << ensembles_.size() << std::endl;
            std::abort();
          }
          for (int i = start; i < end; i++) {
            Ensemble &ens = ensembles_[i];
            // У ансамблей есть "задержка" начала расчетов, если еще не
            // дошли то просто увеличиваем текущий шаг ансамбля
            if (ens.cur_ < 0) {
              ens.cur_++;
              continue;
            }
            // Если текущий шаг 0, то инициализируем ансамбль
            if (ens.cur_ == 0) {
              initEnsemble(thread_finished_count, ens, sys);
            }
            // Диффузия
            double acfv = gk.getACFV(sys, ens.init_velocity_x_,ens.init_velocity_y_, ens.init_velocity_z_);
            ens.accum_acfv_[ens.cur_] = acfv;
            ens.coef_diff_integral_[ens.cur_] =
                integral_diff_.add_partial(acfv) / (3 * sys.particleNumber());
            // Вязкость
            double acfp = gk.getACFP(sys, ens.init_pressure_xy_,ens.init_pressure_yz_, ens.init_pressure_zx_);
            ens.accum_acfp_[ens.cur_] = acfp;
            ens.coef_visc_integral_[ens.cur_] =
                integral_visc_.add_partial(acfp) * sys.dimensions().volume() /
                (3 * sys.temperature());

            writeEnsembleDetailed(ens);

            ens.cur_++;
            // Если ансамбль полностью заполнен, записываем в усредненные
            // значения и обнуляем
            if (ens.cur_ >= ensemble_size_) {
              flushEnsemble(ens);
              thread_finished_count++;
            }
          }
          return thread_finished_count;
        }));
      }

      for (auto &f : futures) {
        current_ens_count += f.get();
      }
    } else {
      finished_ = true;
      for (int i = 0; i < ensemble_size_; i++) {
        avg_ensemble_.coef_diff_integral_[i] /= ensemble_count_;
        avg_ensemble_.coef_visc_integral_[i] /= ensemble_count_;

        avg_ensemble_.accum_acfv_[i] /= ensemble_count_;
        avg_ensemble_.accum_acfp_[i] /= ensemble_count_;
      }
      writeAvgEnsembleDetailed(avg_ensemble_);
    }
  };
};

#endif
