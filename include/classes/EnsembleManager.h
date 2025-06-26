#ifndef ENSEMBLE_MANAGER_H
#define ENSEMBLE_MANAGER_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "nlohmann/json.hpp"

#include "helpers/IncrementalIntegrator.h" //Интегральная функция
#include "classes/Matrix3.h"
#include "classes/ThreadPool.h"
#include "classes/Vector3.h"
#include "core/Settings.h"
#include "core/System.h"
#include "macroparams/Einstein.h"
#include "macroparams/GreenKubo.h"
#include "output/OutputManager.h"

/**
 * @brief Структура для хранения ансамбля.
 *
 * Структура для хранения ансамбля.
 */
struct Ensemble {
  /**
   * @brief Идентификатор ансамбля.
   *
   * Идентификатор ансамбля.
   */
  int id_;

  /**
   * @brief Внутренний счетчик ансамбля.
   *
   * Внутренний счетчик ансамбля.
   */
  int cur_;

  /**
   * @brief Имя файла для записи данных ансамбля.
   *
   * Имя файла для записи данных ансамбля.
   */
  std::string filename_;

  /**
   * @brief Файл для записи данных ансамбля.
   *
   * Файл для записи данных ансамбля.
   */
  std::ofstream output_file_;

  /**
   * @brief Начальные компоненты скорости.
   *
   * Начальные компоненты скорости.
   */
  std::vector<Vector3<double>> init_velocity_{{0.0, 0.0, 0.0}};

  /**
   * @brief Автокорреляционные функции скорости.
   *
   * Автокорреляционные функции скорости.
   */
  std::vector<double> accum_acfv_{0.0};

  /**
   * @brief Интеграл от АКФ для коэффициента диффузии.
   *
   * Интеграл от АКФ для коэффициента диффузии.
   */
  std::vector<double> coef_diff_integral_{0.0};

  /**
   * @brief Начальные тензоры давления системы.
   *
   * Начальные тензоры давления системы.
   */
  Matrix3 init_pressure_tensors_{0.0};

  /**
   * @brief Автокорреляционные функции тензоров давления.
   */
  std::vector<double> accum_acfp_{0.0};

  /**
   * @brief Интеграл от АКФ для коэффициента вязкости.
   *
   * Интеграл от АКФ для коэффициента вязкости.
   */
  std::vector<double> coef_visc_integral_{0.0};
};

using json = nlohmann::json;

/**
 * @brief Класс для управления ансамблями.
 *
 * Класс для управления ансамблями.
 */
class EnsembleManager {
private:
  /**
   * @brief Пул потоков для параллельных вычислений.
   *
   * Пул потоков для параллельных вычислений.
   */
  ThreadPool &threadPool_;

  /**
   * @brief Менеджер вывода данных.
   *
   * Менеджер вывода данных.
   */
  OutputManager &outputManager_;

  /**
   * @brief Интеграл от АКФ для коэффициента диффузии.
   *
   * Интеграл от АКФ для коэффициента диффузии.
   */
  IncrementalIntegrator integral_diff_;

  /**
   * @brief Интеграл от АКФ для коэффициента вязкости.
   *
   * Интеграл от АКФ для коэффициента вязкости.
   */
  IncrementalIntegrator integral_visc_;

  /**
   * @brief Методы Грин-Кубо для транспортных коэффициентов.
   *
   * Методы Грин-Кубо для транспортных коэффициентов.
   */
  GreenKubo gk;

  /**
   * @brief Флаг включения ансамблей.
   *
   * Флаг включения ансамблей.
   */
  bool toggle_{false};

  /**
   * @brief Флаг завершения расчетов.
   *
   * Флаг завершения расчетов.
   */
  bool finished_{false};

  /**
   * @brief Количество шагов.
   *
   * Количество шагов.
   */
  int nsteps_{0};

  /**
   * @brief Количество параллельных "потоков" ансамблей.
   *
   * Количество параллельных "потоков" ансамблей.
   */
  int ensemble_threads_num_;

  /**
   * @brief Размер ансамбля (кол-во шагов).
   *
   * Размер ансамбля (кол-во шагов).
   */
  int ensemble_size_;

  /**
   * @brief "Отступ" между ансамблями.
   *
   * "Отступ" между ансамблями.
   */
  int ensemble_offset_;

  /**
   * @brief Сколько всего полных ансамблей собирать.
   *
   * Сколько всего полных ансамблей собирать.
   */
  int ensemble_count_;

  /**
   * @brief Текущее количество ансамблей.
   *
   * Текущее количество ансамблей.
   */
  int current_ens_count{0};

  /**
   * @brief Вектор ансамблей.
   *
   * Вектор ансамблей.
   */
  std::vector<Ensemble> ensembles_;

  /**
   * @brief Сборник ансамблей.
   */
  Ensemble avg_ensemble_;

  /**
   * @brief Очистка ансамбля.
   *
   * Очистка ансамбля.
   * @param ens - ансамбль.
   */
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

  /**
   * @brief Инициализация ансамбля.
   *
   * Инициализация ансамбля.
   * @param thread_finished_count - количество завершенных потоков.
   * @param ens - ансамбль.
   * @param sys - система.
   */
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
      ens.cur_ = INT32_MIN; // Просто очень большое отрицательное число до
                             // которого не дойдем
      return;
    }
    const int pn = sys.particleNumber();
    const std::vector<Particle> &particles = sys.particles();
    for (int j = 0; j < pn; j++) {
      ens.init_velocity_[j] = particles[j].velocity();
    }
    ens.init_pressure_tensors_ = sys.pressureTensors();
    ens.filename_ = "ensemble_" + std::to_string(ens.id_) + ".csv";
    ens.output_file_ = std::ofstream(
        outputManager_.getEnsembleDir() + "/" + ens.filename_, std::ios::app);
  }

  /**
   * @brief Запись детальных данных ансамбля в файл.
   *
   * Запись детальных данных ансамбля в файл.
   * @param ens - ансамбль.
   */
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

  /**
   * @brief Запись усредненных данных ансамблей в файл.
   *
   * Запись усредненных данных ансамблей в файл.
   * @param avg_ensemble - усредненный ансамбль.
   */
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
  /**
   * @brief Конструктор класса EnsembleManager.
   *
   * Конструктор класса EnsembleManager.
   * @param config - конфигурация.
   * @param settings - настройки.
   * @param sys - система.
   * @param threadPool - пул потоков.
   * @param outputManager - менеджер вывода данных.
   */
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
        ensembles_[i].init_velocity_.resize(sys.particleNumber(), {0, 0, 0});
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

  /**
   * @brief Деструктор по умолчанию.
   *
   * Деструктор по умолчанию.
   */
  ~EnsembleManager() = default;

  /**
   * @brief Проверка включения ансамблей.
   *
   * Проверка включения ансамблей.
   * @return true - ансамбли включены, false - ансамбли выключены.
   */
  inline const bool enabled() const { return toggle_; }

  /**
   * @brief Проверка завершения расчетов.
   *
   * Проверка завершения расчетов.
   * @return true - расчеты завершены, false - расчеты не завершены.
   */
  inline const bool completed() const { return finished_; };

  /**
   * @brief Получение количества параллельных потоков ансамблей.
   *
   * Получение количества параллельных потоков ансамблей.
   * @return количество параллельных потоков ансамблей.
   */
  int getEnsembleThreadsNum() const { return ensemble_threads_num_; }

  /**
   * @brief Получение размера ансамбля.
   *
   * Получение размера ансамбля.
   * @return размер ансамбля.
   */
  int getEnsembleSize() const { return ensemble_size_; }
  
  /**
   * @brief Получение "отступа" между ансамблями.
   *
   * Получение "отступа" между ансамблями.
   * @return "отступ" между ансамблями.
   */
  int getEnsembleOffset() const { return ensemble_offset_; }

  /**
   * @brief Получение количества ансамблей.
   *
   * Получение количества ансамблей.
   * @return количество ансамблей.
   */
  int getEnsembleCount() const { return ensemble_count_; }

  /**
   * @brief Получение завершенного количества ансамблей.
   *
   * Получение завершенного количества ансамблей.
   * @return завершенное количество ансамблей.
   */
  int getCurrentEnsCount() const { return current_ens_count; }
  
  /**
   * @brief Получение вектора ансамблей.
   *
   * Получение вектора ансамблей.
   * @return вектор ансамблей.
   */
  const std::vector<Ensemble> &getEnsembles() const { return ensembles_; }
  
  /**
   * @brief Получение вектора усредненных значений АКФ скорости.
   *
   * Получение вектора усредненных значений АКФ скорости.
   * @return вектор усредненных значений АКФ скорости.
   */
  const std::vector<double> &getEnsemblesAccumAcfv() const {
    return avg_ensemble_.accum_acfv_;
  }

  /**
   * @brief Получение вектора усредненных значений АКФ давления.
   *
   * Получение вектора усредненных значений АКФ давления.
   * @return вектор усредненных значений АКФ давления.
   */
  const std::vector<double> &getEnsemblesAccumAcfp() const {
    return avg_ensemble_.accum_acfp_;
  }

  /**
   * @brief Получение вектора усредненных значений интеграла от АКФ для коэффициента диффузии.
   *
   * Получение вектора усредненных значений интеграла от АКФ для коэффициента диффузии.
   * @return вектор усредненных значений интеграла от АКФ для коэффициента диффузии.
   */
  const std::vector<double> &getEnsemblesCoefDiffIntegral() const {
    return avg_ensemble_.coef_diff_integral_;
  }

  /**
   * @brief Получение вектора усредненных значений интеграла от АКФ для коэффициента вязкости.
   *
   * Получение вектора усредненных значений интеграла от АКФ для коэффициента вязкости.
   * @return вектор усредненных значений интеграла от АКФ для коэффициента вязкости.
   */
  const std::vector<double> &getEnsemblesCoefViscIntegral() const {
    return avg_ensemble_.coef_visc_integral_;
  }

  /**
   * @brief Получение базовой информации о менеджере ансамблей.
   *
   * Получение базовой информации о менеджере ансамблей.
   * @return строка с информацией о менеджере ансамблей.
   */
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
          << "\n\t\tInit Velocity Size: " << ensemble.init_velocity_.size()
          << "\n\t\tAccum ACFV Size: " << ensemble.accum_acfv_.size()
          << "\n\t\tCoef Diff Integral Size: "
          << ensemble.coef_diff_integral_.size()

          << "\n\t\tAccum ACFP Size: " << ensemble.accum_acfp_.size()
          << "\n\t\tCoef Visc Integral Size: "
          << ensemble.coef_visc_integral_.size() << "\n";
    }
    return oss.str();
  }

  /**
   * @brief Расчет транспортных коэффициентов методом Грин-Кубо.
   *
   * Расчет транспортных коэффициентов методом Грин-Кубо.
   * @param sys - система.
   */
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
            double acfv = gk.getACFV(sys, ens.init_velocity_);
            ens.accum_acfv_[ens.cur_] = acfv;
            ens.coef_diff_integral_[ens.cur_] =
                integral_diff_.add_partial(acfv) / (3 * sys.particleNumber());
            // Вязкость
            double acfp = gk.getACFP(sys, ens.init_pressure_tensors_);
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
  }

  // Запрещаем копирование
  EnsembleManager(const EnsembleManager &) = delete;
  EnsembleManager &operator=(const EnsembleManager &) = delete;
};

#endif
