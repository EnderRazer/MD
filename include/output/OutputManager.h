#ifndef OUTPUT_MANAGER_H
#define OUTPUT_MANAGER_H

#include <filesystem>
#include <fstream>
#include <ostream>
#include <sstream>
#include <string>

#include "nlohmann/json.hpp"

#include "classes/Energy.h"
#include "core/Settings.h"
#include "core/System.h"

namespace fs = std::filesystem;
using json = nlohmann::json;

class Ensemble;
class EnsembleManager;

/**
 * @brief Класс для управления выходными данными.
 * @details Класс предоставляет методы для записи данных в файлы и потоки.
 */
class OutputManager {
private:
  /**
   * @brief Путь к основной директории.
   */
  std::string output_dir_;

  /**
   * @brief Путь к директории с шагами.
   */
  std::string step_dir_;

  /**
   * @brief Путь к директории с ансамблями.
   */
  std::string ens_dir_;

  /**
   * @brief Частота записи данных.
   */
  int frequency_;

public:
  /**
   * @brief Конструктор.
   * @param config - конфигурация.
   * @param settings - настройки.
   */
  OutputManager(const json &config, Settings &settings) {
    output_dir_ = config.value("main_directory", "data") + "_" +
                  std::to_string(settings.seed());
    if (!fs::exists(output_dir_)) {
      fs::create_directory(output_dir_);
    }
    step_dir_ =
        output_dir_ + "/" + config.value("steps_directory", "steps_data");
    if (!fs::exists(step_dir_)) {
      fs::create_directory(step_dir_);
    }
    ens_dir_ = output_dir_ + "/" +
               config.value("ensembles_directory", "ensembles_data");
    if (!fs::exists(ens_dir_)) {
      fs::create_directory(ens_dir_);
    }
    frequency_ = config.value("frequency", 1);
  }

  // Запрещаем копирование
  OutputManager(const OutputManager &) = delete;
  OutputManager &operator=(const OutputManager &) = delete;

  inline const int frequency() const { return frequency_; }
  inline const std::string &getEnsembleDir() const { return ens_dir_; }
  inline const std::string &getStepDir() const { return step_dir_; }
  inline const std::string &getOutputDir() const { return output_dir_; }

  /**
   * @brief Запись данных в файл.
   * @param filename - имя файла.
   * @param headers - заголовки.
   * @param data - данные.
   * @param append - флаг добавления данных.
   * @param delimiter - разделитель.
   * @param precision - точность.
   */
  template <typename T>
  void writeDataToFile(const std::string &filename,
                       const std::vector<std::string> &headers,
                       const std::vector<std::vector<T>> &data,
                       bool append = false, const std::string &delimiter = ";",
                       int precision = 16) const {
    // Check if file exists (for header writing decision)
    bool file_exists = fs::exists(filename);

    // Open file in appropriate mode
    std::ofstream file(filename, append ? std::ios::app : std::ios::out);

    writeToStream(file, headers, data, !file_exists || !append, delimiter,
                  precision);
  }

  /**
   * @brief Запись данных в поток.
   * @param stream - поток.
   * @param headers - заголовки.
   * @param data - данные.
   * @param write_headers - флаг записи заголовков.
   * @param delimiter - разделитель.
   * @param precision - точность.
   */
  template <typename T>
  void writeDataToStream(std::ofstream &stream,
                         const std::vector<std::string> &headers,
                         const std::vector<std::vector<T>> &data,
                         bool write_headers = true,
                         const std::string &delimiter = ";",
                         int precision = 16) const {
    writeToStream(stream, headers, data, write_headers, delimiter, precision);
  }

  /**
   * @brief Запись данных в поток.
   * @param stream - поток.
   * @param headers - заголовки.
   * @param data - данные.
   * @param write_headers - флаг записи заголовков.
   * @param delimiter - разделитель.
   * @param precision - точность.
   */
  template <typename T, typename StreamType>
  void writeToStream(StreamType &stream, const std::vector<std::string> &headers,
                const std::vector<std::vector<T>> &data, bool write_headers,
                const std::string &delimiter, int precision) const {
    // Set precision for floating point numbers
    stream << std::fixed << std::setprecision(precision);

    // Write header if requested
    if (write_headers) {
      for (size_t i = 0; i < headers.size(); ++i) {
        stream << headers[i];
        if (i < headers.size() - 1) {
          stream << delimiter;
        }
      }
      stream << "\n";
    }

    // Write data rows
    for (const auto &row : data) {
      for (size_t i = 0; i < row.size(); ++i) {
        stream << row[i];
        if (i < row.size() - 1) {
          stream << delimiter;
        }
      }
      stream << "\n";
    }
  }

  /**
   * @brief Запись строки в файл.
   * @param filename - имя файла.
   * @param headers - заголовки.
   * @param row - строка.
   * @param append - флаг добавления данных.
   * @param delimiter - разделитель.
   * @param precision - точность.
   */
  template <typename T>
  void writeRowToFile(const std::string &filename,
                      const std::vector<std::string> &headers,
                      const std::vector<T> &row, bool append = false,
                      const std::string &delimiter = ";",
                      int precision = 16) const {
    std::vector<std::vector<T>> data = {row};
    writeDataToFile(filename, headers, data, append, delimiter, precision);
  }

  /**
   * @brief Запись строки в поток.
   * @param stream - поток.
   * @param headers - заголовки.
   * @param row - строка.
   * @param write_headers - флаг записи заголовков.
   * @param delimiter - разделитель.
   * @param precision - точность.
   */
  template <typename T>
  void writeRowToStream(std::ofstream &stream,
                        const std::vector<std::string> &headers,
                        const std::vector<T> &row, bool write_headers = true,
                        const std::string &delimiter = ";",
                        int precision = 16) const {
    std::vector<std::vector<T>> data = {row};
    writeDataToStream(stream, headers, data, write_headers, delimiter,
                      precision);
  }

  /**
   * @brief Запись конфигурации расчетов в файл.
   * @param config - конфигурация.
   */
  void writeModellingProperties(json &config) const {
    std::string filename = output_dir_ + "/launch_config.json";
    std::ofstream file(filename, std::ios::app);
    file << config.dump(4) << std::endl;
    file.close();
  }

  /**
   * @brief Запись свойств системы в файл.
   * @param sys - система.
   */
  void writeSystemProperties(System &sys) const {
    std::string filename = output_dir_ + "/properties.csv";

    std::vector<std::string> headers = {"Step",      "Temperature", "Pressure",
                                        "KineticE",  "PotentialE",  "ThermalE",
                                        "InternalE", "TotalE",      "Pulse"};

    std::vector<double> data = {
        static_cast<double>(sys.currentStep()),
        sys.temperature(),
        sys.pressure(),
        sys.energies().get(Energy::EnergyType::Kinetic),
        sys.energies().get(Energy::EnergyType::Potential),
        sys.energies().get(Energy::EnergyType::Thermodynamic),
        sys.energies().get(Energy::EnergyType::Internal),
        sys.energies().get(Energy::EnergyType::Full),
        sys.pulse()};

    writeRowToFile(filename, headers, data, true);
  }

  /**
   * @brief Запись средних свойств системы в файл.
   * @param sys - система.
   */
  void writeAvgSystemProperties(System &sys) const {
    std::string filename = output_dir_ + "/properties_avg.csv";

    std::vector<std::string> headers = {"Step",      "Temperature", "Pressure",
                                        "KineticE",  "PotentialE",  "ThermalE",
                                        "InternalE", "TotalE",      "Pulse"};

    std::vector<double> data = {
        static_cast<double>(sys.currentStep()),
        sys.temperatureAvg(),
        sys.pressureAvg(),
        sys.energiesAvg().get(Energy::EnergyType::Kinetic),
        sys.energiesAvg().get(Energy::EnergyType::Potential),
        sys.energiesAvg().get(Energy::EnergyType::Thermodynamic),
        sys.energiesAvg().get(Energy::EnergyType::Internal),
        sys.energiesAvg().get(Energy::EnergyType::Full),
        sys.pulseAvg()};

    writeRowToFile(filename, headers, data, true);
  }

  /**
   * @brief Запись пошаговых данных о частицах в файл.
   * @param sys - система.
   */
  void writeStepData(System &sys) const {
    std::string filename = step_dir_ + "/step_" +
                           std::to_string(sys.currentStep()) + "_particles.csv";

    std::vector<std::string> headers = {
        "n",          "Cx",       "Cy",       "Cz",        "Vx",
        "Vy",         "Vz",       "Fx",       "Fy",        "Fz",
        "PotentialE", "KineticE", "ThermalE", "InternalE", "TotalE"};

    std::vector<std::vector<double>> data;
    int id = 0;

    for (const auto &particle : sys.particles()) {
      const auto &id = particle.id();
      const auto &coord = particle.coord();
      const auto &velocity = particle.velocity();
      const auto &force = particle.force();
      const auto &energies = particle.energies();

      std::vector<double> row = {
          static_cast<double>(id),
          coord.x(),
          coord.y(),
          coord.z(),
          velocity.x(),
          velocity.y(),
          velocity.z(),
          force.x(),
          force.y(),
          force.z(),
          energies.get(Energy::EnergyType::Potential),
          energies.get(Energy::EnergyType::Kinetic),
          energies.get(Energy::EnergyType::Thermodynamic),
          energies.get(Energy::EnergyType::Internal),
          energies.get(Energy::EnergyType::Full)};

      data.push_back(row);
    }

    writeDataToFile(filename, headers, data);
  }

  /**
   * @brief Получение базовой информации о выходном менеджере.
   * @return Базовая информация о выходном менеджере.
   */
  const std::string getData() const {
    std::ostringstream oss;
    oss << "Output manager data:"
        << "\n\tMain directory: " << output_dir_ << "\n\tSteps directory"
        << step_dir_ << "\n\tEnsembles directory: " << ens_dir_ << std::endl;
    return oss.str();
  }
};

#endif // OUTPUT_MANAGER_H
