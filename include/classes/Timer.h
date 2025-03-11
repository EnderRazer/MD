#ifndef TIMER_H
#define TIMER_H

#include <chrono>

class Timer {
public:
  Timer(int mode) : running(false) { this->mode = mode; }

  inline void start() {
    if (!running) {
      start_time = std::chrono::high_resolution_clock::now();
      this->running = true;
    }
  }

  inline void stop() {
    if (running) {
      end_time = std::chrono::high_resolution_clock::now();
      this->running = false;
    }
  }

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
  bool running;
  int mode;
  std::chrono::high_resolution_clock::time_point start_time;
  std::chrono::high_resolution_clock::time_point end_time;
};

#endif // TIMER_H