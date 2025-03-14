#ifndef ENSEMBLE_MANAGER_H
#define ENSEMBLE_MANAGER_H

struct Ensemble
{
  int id_;
  int cur_; // Внутренний счетчик ансамбля

  std::vector<Vector3<double>>                  // Начальные компоненты скорости
      init_velocity_{{0.0, 0.0, 0.0}};          // (для расчетов коэффициента диффузии)
  std::vector<double> accum_acfv_{0.0};         // Автокорреляционные функции скорости
  std::vector<double> coef_diff_integral_{0.0}; // Интеграл от АКФ

  Matrix3 init_pressure_tensors_{0.0}; // Начальные тензоры давления системы
  std::vector<double> accum_acfp_{
      0.0};                                     // Автокорреляционные функции тензоров давления
  std::vector<double> coef_visc_integral_{0.0}; // Интеграл от АКФ
};

#include "classes/IncrementalIntegrator.h" //Интегральная функция
#include "macroparams/GreenKubo.h"
#include "macroparams/Einstein.h"

using json = nlohmann::json;

class EnsembleManager
{
private:
  bool finished_{false};
  bool toggle_{false};
  ThreadPool &threadPool_;
  IncrementalIntegrator integral_diff_;
  IncrementalIntegrator integral_visc_;
  GreenKubo gk;
  int nsteps_;
  int ensemble_threads_num_; // Количество параллельных "потоков" ансамблей
  int ensemble_size_;        // Размер ансамбля (кол-во шагов)
  int ensemble_offset_;      // Отступ между ансамблями
  int ensemble_count_;       // Сколько всего полных ансамблей собирать

  int current_ens_count{0};
  std::vector<Ensemble> ensembles_;

  std::vector<double> ensembles_accum_acfv_{0.0};
  std::vector<double> ensembles_accum_acfp_{0.0};

  std::vector<double> ensembles_coef_diff_integral_{0.0};
  std::vector<double> ensembles_coef_visc_integral_{0.0};

  inline void flushEnsemble(Ensemble &ens)
  {
    ens.cur_ = 0;
    for (int i = 0; i < ensemble_size_; i++)
    {
      // Записываем
      ensembles_accum_acfv_[i] += ens.accum_acfv_[i];
      ensembles_accum_acfp_[i] += ens.accum_acfp_[i];
      ensembles_coef_diff_integral_[i] += ens.coef_diff_integral_[i];
      ensembles_coef_visc_integral_[i] += ens.coef_visc_integral_[i];
      // Очищаем
      ens.accum_acfv_[i] = 0.0;
      ens.accum_acfp_[i] = 0.0;
      ens.coef_diff_integral_[i] = 0.0;
      ens.coef_visc_integral_[i] = 0.0;
    }
  };

  inline void initEnsemble(int thread_finished_count, Ensemble &ens,
                           System &sys)
  {
    // Если не успеваем рассчитать полный ансамбль, то и не считаем, чтобы не
    // забирать ресурсы
    ens.id_ += thread_finished_count + current_ens_count;
    if (sys.currentStep() + ensemble_size_ > nsteps_)
    {
      ens.cur_ = -INT32_MAX; // Просто очень большое отрицательное число до
                             // которого не дойдем
      return;
    }
    const int pn = sys.particleNumber();
    const std::vector<Particle> &particles = sys.particles();
    for (int j = 0; j < pn; j++)
    {
      ens.init_velocity_[j] = particles[j].velocity();
    }
    ens.init_pressure_tensors_ = sys.pressureTensors();
  }

public:
  EnsembleManager(json &config, Settings &settings, System &sys,
                  ThreadPool &threadPool)
      : integral_diff_(settings.dt()), integral_visc_(settings.dt()),
        threadPool_(threadPool)
  {
    toggle_ = config.value("toggle", false);
    if (toggle_)
    {
      nsteps_ = settings.steps();
      ensemble_threads_num_ = config.value("threads", 1);
      ensemble_size_ = config.value("size", 1000000);
      ensemble_offset_ = config.value("offset", 1000000);
      ensemble_count_ = config.value("quantity", 100);
      ensembles_.resize(ensemble_threads_num_);
      ensembles_accum_acfv_.resize(ensemble_size_, 0.0);
      ensembles_accum_acfp_.resize(ensemble_size_, 0.0);

      ensembles_coef_diff_integral_.resize(ensemble_size_);
      ensembles_coef_visc_integral_.resize(ensemble_size_);

      int expected_contribution = 0;
      for (int i = 0; i < ensemble_threads_num_; i++)
      {
        ensembles_[i].id_ = i;
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
      if (expected_contribution < ensemble_count_)
      {
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
  const std::vector<double> &getEnsemblesAccumAcfv() const
  {
    return ensembles_accum_acfv_;
  }
  const std::vector<double> &getEnsemblesAccumAcfp() const
  {
    return ensembles_accum_acfp_;
  }
  const std::vector<double> &getEnsemblesCoefDiffIntegral() const
  {
    return ensembles_coef_diff_integral_;
  }
  const std::vector<double> &getEnsemblesCoefViscIntegral() const
  {
    return ensembles_coef_visc_integral_;
  }
  std::string getData() const
  {
    std::ostringstream oss;
    oss << "Using ensembles to calculate transport coefficients"
        << "\n\tКоличество параллельных потоков ансамблей: "
        << ensemble_threads_num_
        << "\n\tКоличество шагов в ансамбле: " << ensemble_size_
        << "\n\tОтступ между ансамблями (шагов): " << ensemble_offset_
        << "\n\tКоличество рассчитанных ансамблей: " << ensemble_count_
        << "\n\tEnsembles Size: " << ensembles_.size();

    for (const auto &ensemble : ensembles_)
    {
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

  void accumulateGreenCubo(System &sys)
  {
    if (current_ens_count < ensemble_count_)
    {
      const int blockSize =
          ensemble_threads_num_ /
          threadPool_.getThreads(); // или любая удобная величина
      std::vector<std::future<int>> futures;

      for (int start = 0; start < ensemble_threads_num_; start += blockSize)
      {
        int end = std::min(start + blockSize, ensemble_threads_num_);
        // enqueue задачу на вычисление для блока [start, end)
        futures.push_back(threadPool_.enqueue([this, &sys, start, end]()
                                              {
          int thread_finished_count = 0;
          for (int i = start; i < end; i++) {
            Ensemble &ens = ensembles_[i];
            // У ансамблей есть "задержка" начала расчетов, если еще не
            // дошли то просто увеличиваем текущий шаг ансамбля
            if (ens.cur_ < 0) {
              ens.cur_++;
              continue;
            }
            // Если ансамбль полностью заполнен, записываем в усредненные
            // значения и обнуляем
            if (ens.cur_ >= ensemble_size_) {
              flushEnsemble(ens);
              thread_finished_count++;
            }
            // Если текущий шаг 0, то инициализируем ансамбль
            if (ens.cur_ == 0) {
              initEnsemble(thread_finished_count, ens, sys);
            }
            // Диффузия
            double acfv = gk.getACFV(sys, ens);
            ens.accum_acfv_[ens.cur_] = acfv;
            ens.coef_diff_integral_[ens.cur_] =
                integral_diff_.add_partial(acfv) / (3 * sys.particleNumber());
            // Вязкость
            double acfp = gk.getACFP(sys, ens);
            ens.accum_acfp_[ens.cur_] = acfp;
            ens.coef_visc_integral_[ens.cur_] =
                integral_visc_.add_partial(acfp) * sys.dimensions().volume() /
                (3 * sys.temperature());

            ens.cur_++;
          }
          return thread_finished_count; }));
      }

      for (auto &f : futures)
      {
        current_ens_count += f.get();
      }
    }
    else
    {
      finished_ = true;
      for (int i = 0; i < ensemble_size_; i++)
      {
        ensembles_accum_acfp_[i] /= ensemble_count_;
        ensembles_accum_acfv_[i] /= ensemble_count_;
        ensembles_coef_diff_integral_[i] /= ensemble_count_;
        ensembles_coef_visc_integral_[i] /= ensemble_count_;
      }
    }
  };
};

#endif