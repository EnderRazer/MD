#ifndef INCREMENTAL_INTEGRATOR_H
#define INCREMENTAL_INTEGRATOR_H

/**
 * @brief Класс для интегрирования по времени.
 *
 * Класс для интегрирования по времени.
 */
class IncrementalIntegrator {
private:
  /**
   * @brief Интеграл.
   *
   * Интеграл.
   */
  double accumulated_integral_ = 0.0;

  /**
   * @brief Шаг интегрирования.
   *
   * Шаг интегрирования.
   */
  double h_;

  /**
   * @brief Последнее значение.
   *
   * Последнее значение.
   */
  double last_value_ = 0.0;

  
public:
  /**
   * @brief Конструктор класса IncrementalIntegrator.
   *
   * Конструктор класса IncrementalIntegrator.
   * @param step_size - шаг интегрирования.
   */
  IncrementalIntegrator(double step_size) : h_(step_size) {}

  /**
   * @brief Добавление частичного интеграла.
   *
   * Добавление частичного интеграла.
   * @param value - значение.
   * @return интеграл.
   */
  double add_partial(const double value) {
    double local_integral = (last_value_ + value) * h_ / 2.0;

    // Update state
    accumulated_integral_ += local_integral;
    last_value_ = value; // Store last function value

    return accumulated_integral_;
  }

  /**
   * @brief Получение общего интеграла.
   *
   * Получение общего интеграла.
   * @return интеграл.
   */
  double get_total_integral() const { return accumulated_integral_; }
};
#endif