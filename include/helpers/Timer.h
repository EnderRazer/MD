#ifndef TIMER_H
#define TIMER_H

#include <chrono>

/**
 * @brief Класс для измерения времени выполнения кода.
 * @details Класс предоставляет интерфейс для измерения времени выполнения кода.
 * Он позволяет запускать и останавливать таймер, а также получать время,
 * прошедшее с момента запуска таймера.
 */
class Timer {
public:
  /**
   * @brief Конструктор класса Timer.
   * @param mode - режим работы таймера.
   * @details Конструктор создает таймер с заданным режимом работы.
   */
  Timer(int mode) : running(false) { this->mode = mode; }

  /**
   * @brief Запуск таймера.
   * @details Метод запускает таймер.
   */
  inline void start() {
    if (!running) {
      start_time = std::chrono::high_resolution_clock::now();
      this->running = true;
    }
  }

  /**
   * @brief Остановка таймера.
   * @details Метод останавливает таймер.
   */
  inline void stop() {
    if (running) {
      end_time = std::chrono::high_resolution_clock::now();
      this->running = false;
    }
  }

  /**
   * @brief Получение времени, прошедшего с момента запуска таймера.
   * @details Метод возвращает время, прошедшее с момента запуска таймера.
   * @return Время, прошедшее с момента запуска таймера.
   */
  inline double elapsed() const {
    if (running) {
      auto current_time = std::chrono::high_resolution_clock::now();
      if (mode == 1)
        return std::chrono::duration<double, std::milli>(current_time -
                                                         start_time)
            .count();

      return std::chrono::duration<double>(current_time - start_time).count();
    } else {
      if (mode == 1)
        return std::chrono::duration<double, std::milli>(end_time - start_time)
            .count();

      return std::chrono::duration<double>(end_time - start_time).count();
    }
  }

private:
  /**
   * @brief Флаг работы таймера.
   */
  bool running;

  /**
   * @brief Режим работы таймера.
   */
  int mode;

  /**
   * @brief Время начала работы таймера.
   */
  std::chrono::high_resolution_clock::time_point start_time;

  /**
   * @brief Время окончания работы таймера.
   */
  std::chrono::high_resolution_clock::time_point end_time;
};

#endif // TIMER_H