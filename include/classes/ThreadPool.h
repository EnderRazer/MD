#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <condition_variable>
#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>

class ThreadPool {
public:
  // Конструктор: создаёт пул с заданным количеством потоков.
  // По умолчанию создаётся количество потоков = числу аппаратных ядер (если
  // возможно определить).
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

  // Запрещаем копировать пул потоков (можно разрешить, но придётся аккуратно
  // реализовывать).
  ThreadPool(const ThreadPool &) = delete;
  ThreadPool &operator=(const ThreadPool &) = delete;

  // Деструктор — останавливаем все потоки
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

  // Функция добавления (enqueue) новой задачи в пул.
  // Возвращает future, позволяющий получить результат задачи.
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
  const int getThreads() const { return threads_; };
  const std::string getData() const {
    std::ostringstream oss;
    oss << "Thread pool data:" << "\n\tThreads: " << threads_ << std::endl;
    return oss.str();
  }

private:
  int threads_{0};
  // Наши рабочие потоки
  std::vector<std::thread> workers_;
  // Очередь задач (каждая задача — это std::function<void()>)
  std::queue<std::function<void()>> tasks_;

  // Синхронизация
  std::mutex queue_mutex_;
  std::condition_variable condition_;
  // Флаг остановки пула
  bool stop_;
};

#endif // THREAD_POOL_H