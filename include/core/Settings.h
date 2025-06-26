#ifndef SETTINGS_H
#define SETTINGS_H

#include "nlohmann/json.hpp"
#include <sstream>

using json = nlohmann::json;

/**
 * @brief Константы.
 * @details Константы для расчета.
 */
struct Constants {
  const double D{3};                   // Число степеней свободы
  const double KBOLTZMAN{1.380648528}; // Постоянная Больцмана
};

/**
 * @brief Класс для работы с настройками расчета.
 * @details Класс предоставляет методы для работы с настройками расчета.
 */
class Settings {
private:
  /**
   * @brief Константы.
   */
  Constants constants_;

  /**
   * @brief Сид.
   */
  int seed_{0};

  /**
   * @brief Число потоков.
   */
  int threads_{0};

  /**
   * @brief Флаг отладки.
   */
  bool debug_{false};

  /**
   * @brief Флаг профилирования.
   */
  bool profiling_{false};

  /**
   * @brief Флаг периодических граничных условий.
   */
  bool pbc_{false};

  /**
   * @brief Число шагов.
   */
  int nsteps_{0};

  /**
   * @brief Шаг по времени.
   */
  double delta_t_{0.0};

  /**
   * @brief Начальная температура.
   */
  double start_temp_{0.0};

  /**
   * @brief Название материала.
   */
  std::string material_name_{"Undefined"};

  /**
   * @brief Тип структуры материала.
   */
  std::string struct_type_{"Undefined"};

  /**
   * @brief Пользовательская структура материала.
   */
  json custom_layout_;

  /**
   * @brief Масса материала.
   */
  double mass_{0.0};

  
public:
  
  Settings() = default;
  ~Settings() = default;

  // Запрещаем копирование
  Settings(const Settings &) = delete;
  Settings &operator=(const Settings &) = delete;

  /**
   * @brief Конструктор.
   * @details Конструктор класса Settings.
   */
  Settings(json &config) {
    if (!config.contains("settings"))
      throw std::invalid_argument("В файле конфигурации нет настроек системы");
    json settings = config["settings"];
    // Важные переменные
    if (!settings.contains("pbc"))
      throw std::invalid_argument("В файле конфигурации нет настроек ПГУ");
    pbc_ = settings.value("pbc", true);

    if (!settings.contains("start_temp"))
      throw std::invalid_argument(
          "В файле конфигурации нет настройки начальной температуры системы");
    start_temp_ = settings.value("start_temp", 2.7315);

    if (!settings.contains("nsteps"))
      throw std::invalid_argument(
          "В файле конфигурации нет настройки количества шагов");
    nsteps_ = settings.value("nsteps", 1);

    if (!settings.contains("delta_t"))
      throw std::invalid_argument(
          "В файле конфигурации нет настройки дельты по времени");
    delta_t_ = settings.value("delta_t", 0.001);

    if (!config.contains("material"))
      throw std::invalid_argument("В файле конфигурации нет настроек вещества");
    json material = config["material"];

    if (!material.contains("struct_type"))
      throw std::invalid_argument(
          "В файле конфигурации нет настроек структуры вещества");
    struct_type_ = material.value("struct_type", "Undefined");

    if (!material.contains("mass"))
      throw std::invalid_argument(
          "В файле конфигурации нет настроек массы вещества");
    mass_ = material.value("mass", 0.0);

    if (struct_type_ == "Custom")
      if (!config.contains("custom"))
        throw std::invalid_argument(
            "В файле конфигурации нет настроек пользовательской структуры");
    custom_layout_ = config.value("custom", json::object());

    // Неважные переменные
    threads_ = settings.value("threads", 1);
    // omp_set_num_threads(threads_);
    debug_ = settings.value("debug", false);
    profiling_ = settings.value("profiling", false);

    srand(time(0));
    seed_ = rand();
  };
  
  // Getter and Setter for threads
  inline const int threads() const { return threads_; }
  inline void setThreads(int t) { threads_ = t; }

  // Getter and Setter for seed
  inline const int seed() const { return seed_; }

  // Getter and Setter for dedug (debug flag)
  inline const bool isDebug() const { return debug_; }
  inline void setDebug(bool d) { debug_ = d; }

  // Getter and Setter for profiling
  inline const bool isProfiling() const { return profiling_; }
  inline void setProfiling(bool p) { profiling_ = p; }

  // Getter and Setter for periodic boundary conditions (pbc)
  inline const bool hasPbc() const { return pbc_; }
  inline void setPbc(bool value) { pbc_ = value; }

  // Getter and Setter for start_temp
  inline const double startTemp() const { return start_temp_; }
  inline void setStartTemp(double temp) { start_temp_ = temp; }

  // Getter and Setter for nsteps
  inline const int steps() const { return nsteps_; }
  inline void setSteps(int steps) { nsteps_ = steps; }

  // Getter and Setter for delta_t
  inline const double dt() const { return delta_t_; }
  inline void setDt(double delta) { delta_t_ = delta; }

  // Getter and Setter for material name
  inline const std::string materialName() const { return material_name_; }
  inline void setMaterialName(std::string name) { material_name_ = name; }

  // Getter and Setter for material struct_type
  inline const std::string structType() const { return struct_type_; }
  inline void setStructType(std::string type) { struct_type_ = type; }

  // Getter and Setter for material mass
  inline const double mass() const { return mass_; }
  inline void setMass(double massa) { mass_ = massa; }

  // Getter for constants
  inline const Constants constants() const { return constants_; }

  inline const json customLayout() const { return custom_layout_; }

  /**
   * @brief Получение базовой информации о настройках.
   * @return Базовая информация о настройках.
   */
  std::string getData() const {
    std::ostringstream oss;
    oss.precision(16);
    oss << "Settings data:"
        << "\n\tDebug mode: " << debug_ << "\n\tProfiling: " << profiling_
        << "\n\tPBC mode: " << pbc_ << "\n\tThreads number: " << threads_
        << "\n\tdt: " << delta_t_ << "\n\tSteps count: " << nsteps_
        << "\n\tStarting temperature: " << start_temp_ << "\n";
    return oss.str();
  }
};

#endif // SETTINGS_H
