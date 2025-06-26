#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <condition_variable>
#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>

/**
 * @brief Класс для управления пулом потоков.
 * @details Класс предоставляет интерфейс для управления пулом потоков.
 * Он позволяет добавлять задачи в очередь и получать результаты выполнения
 * задач.
 */
class ThreadPool {
public:
  /**
   * @brief Конструктор класса ThreadPool.
   * @param threads - количество потоков в пуле.
   * @details Конструктор создает пул потоков с заданным количеством потоков.
   * По умолчанию создается количество потоков = числу аппаратных ядер (если
   * возможно определить).
   */
  explicit ThreadPool(size_t threads = std::thread::hardware_concurrency())
      : stop_(false) {
    if (threads == 0) {
      throw std::runtime_error(
          "ThreadPool: количество потоков не может быть 0");
    }
    threads_ = threads;

    // Создаём заданное число потоков, каждый поток запускает worker()
    for (size_t i = 0; i < threads; ++i) {
      workers_.emplace_back([this] {
        // Пока не попросили остановиться, извлекаем задачи из очереди и
        // выполняем
        for (;;) {
          std::function<void()> task;
          {
            // Уникальный лок для работы с очередью
            std::unique_lock<std::mutex> lock(queue_mutex_);
            // Ждём, пока в очереди появится задача, либо поступит сигнал на
            // завершение
            condition_.wait(lock, [this] { return stop_ || !tasks_.empty(); });
            // Если попросили остановиться и задач больше нет — выходим из цикла
            if (stop_ && tasks_.empty()) {
              return;
            }
            // Извлекаем задачу из очереди
            task = std::move(tasks_.front());
            tasks_.pop();
          }
          // Выполняем задачу
          task();
        }
      });
    }
  }

  /**
   * @brief Деструктор класса ThreadPool.
   * @details Деструктор останавливает все потоки в пуле.
   */
  ~ThreadPool() {
    {
      // Ставим флаг stop_ = true
      std::unique_lock<std::mutex> lock(queue_mutex_);
      stop_ = true;
    }
    // Будим все потоки, чтобы они могли выйти из worker()
    condition_.notify_all();
    // Ждём их завершения
    for (std::thread &worker : workers_) {
      worker.join();
    }
  }

  /**
   * @brief Добавление новой задачи в пул.
   * @param f - функция, которая будет выполнена.
   * @param args - аргументы функции.
   * @return future, позволяющий получить результат задачи.
   */
  template <class F, class... Args>
  auto enqueue(F &&f, Args &&...args)
      -> std::future<typename std::invoke_result<F, Args...>::type> {
    using return_type = typename std::invoke_result<F, Args...>::type;

    // Упаковываем функцию и аргументы в packaged_task
    auto task_ptr = std::make_shared<std::packaged_task<return_type()>>(
        std::bind(std::forward<F>(f), std::forward<Args>(args)...));

    std::future<return_type> res = task_ptr->get_future();
    {
      std::unique_lock<std::mutex> lock(queue_mutex_);

      // Если пул остановлен, бросаем исключение
      if (stop_) {
        throw std::runtime_error("enqueue on stopped ThreadPool");
      }
      // Добавляем "лямбду-вызов" в очередь задач
      tasks_.emplace([task_ptr]() { (*task_ptr)(); });
    }
    // Будим один поток
    condition_.notify_one();
    return res;
  }

  /**
   * @brief Получение количества потоков в пуле.
   * @return Количество потоков в пуле.
   */
  const int getThreads() const { return threads_; };

  /**
   * @brief Получение данных о пуле потоков.
   * @return Данные о пуле потоков.
   */
  const std::string getData() const {
    std::ostringstream oss;
    oss << "Thread pool data:" << "\n\tThreads: " << threads_ << std::endl;
    return oss.str();
  }

  // Запрещаем копировать пул потоков (можно разрешить, но придётся аккуратно
  // реализовывать).
  ThreadPool(const ThreadPool &) = delete;
  ThreadPool &operator=(const ThreadPool &) = delete;

private:
  /**
   * @brief Количество потоков в пуле.
   */
  int threads_{0};

  /**
   * @brief Наши рабочие потоки.
   */
  std::vector<std::thread> workers_;

  /**
   * @brief Очередь задач (каждая задача — это std::function<void()>)
   */
  std::queue<std::function<void()>> tasks_;

  /**
   * @brief Синхронизация.
   */
  std::mutex queue_mutex_;

  /**
   * @brief Условная переменная.
   */
  std::condition_variable condition_;

  /**
   * @brief Флаг остановки пула.
   */
  bool stop_;
};

#endif // THREAD_POOL_H