#ifndef BACKUP_MANAGER_H
#define BACKUP_MANAGER_H

#include "nlohmann/json.hpp"
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <stdio.h>

#include "classes/Particle.h"
#include "classes/Vector3.h"
#include "core/System.h"

using json = nlohmann::json;
/**
 * @brief Класс для управления резервными копиями данных.
 *
 * Класс BackupManager отвечает за создание и восстановление резервных копий
 * данных системы. Он позволяет включать и отключать резервное копирование,
 * устанавливать частоту создания резервных копий и восстанавливать данные
 * из резервных копий.
 */
class BackupManager {
public:

  /**
   * @brief Получение состояния включения резервного копирования.
   *
   * Получение состояния включения резервного копирования.
   * @return true, если резервное копирование включено, false в противном случае.
   */
  inline bool enabled() const { return toggle_; }

  /**
   * @brief Получение частоты создания резервных копий.
   *
   * Получение частоты создания резервных копий.
   * @return частота создания резервных копий.
   */
  inline int frequency() const { return backup_frequency_; }

  /**
   * @brief Создание резервной копии.
   *
   * Создание резервной копии.
   * @param sys - система, данные которой будут записаны в резервную копию.
   * @return true, если резервная копия создана успешно, false в противном случае.
   */
  bool createBackup(System &sys) const {
    std::ostringstream oss;
    oss << base_folder_ << "/"
        << "backup_" << sys.currentStep() << ".bak";
    cleanupBackup(sys.currentStep());
    return writeRecordsToFile(oss.str(), sys);
  }

  /**
   * @brief Восстановление данных из резервной копии.
   *
   * Восстановление данных из резервной копии.
   * @param sys - система, данные которой будут восстановлены из резервной копии.
   * @return true, если данные восстановлены успешно, false в противном случае.
   */
  bool restoreBackup(System &sys) {
    std::string input;
    std::cout << "Введите путь до файла бекапа: ";
    std::cin >> input;
    return readRecordsFromFile(
        input + "backup_" + std::to_string(backup_restore_step_) + ".bak", sys);
  };

  /**
   * @brief Получение базовой информации о резервном копировании.
   *
   * Получение базовой информации о резервном копировании.
   * @return строка с информацией о резервном копировании.
   */
  std::string getData() const {
    std::ostringstream oss;
    oss << "Backupdata"
        << "\n\tBackup enabled: " << toggle_
        << "\n\tBackup restore step: " << backup_restore_step_
        << "\n\tBackup frequency: " << backup_frequency_ << "\n";
    return oss.str();
  }

  /**
   * @brief Конструктор класса BackupManager.
   *
   * Конструктор класса BackupManager принимает конфигурацию и настройки системы.
   * 
   */
  BackupManager(json config, Settings &settings) {
    base_folder_ = "data_" + std::to_string(settings.seed()) + "/backups";
    if (!std::filesystem::exists(base_folder_)) {
      if (!std::filesystem::create_directories(base_folder_)) {
        std::cerr << "Failed to create backup directory: " << base_folder_
                  << std::endl;
      }
    }
    toggle_ = config.value("toggle", false);
    backup_frequency_ = config.value("frequency", 100);
    backup_restore_step_ = config.value("restore_step", 0);
  };
  /**
   * @brief Деструктор класса BackupManager.
   *
   * Деструктор класса BackupManager по умолчанию.
   */
  ~BackupManager() = default;

private:
  /**
   * @brief Базовый путь для резервных копий.
   *
   * Базовый путь для резервных копий данных системы.
   */
  std::string base_folder_{""};

  /**
   * @brief Флаг включения резервного копирования.
   *
   * Флаг, указывающий, включено ли резервное копирование.
   */
  bool toggle_{false};

  /**
   * @brief Частота создания резервных копий.
   *
   * Частота создания резервных копий в шагах.
   */
  int backup_frequency_{100};

  /**
   * @brief Шаг восстановления резервной копии.
   *
   * Шаг, на котором будет произведено восстановление резервной копии.
   */
  int backup_restore_step_{0};

  /**
   * @brief Вспомогательная функция для записи записей в двоичный файл.
   *
   * Вспомогательная функция для записи записей в двоичный файл.
   *
   * @param filename Имя файла для записи.
   * @param sys Система, данные которой будут записаны.
   * @return true, если запись прошла успешно, false в противном случае.
   */
  bool writeRecordsToFile(const std::string &filename, System &sys) const {
    // Open file in binary mode, std::ios::out ensures we overwrite
    std::ofstream outFile(filename, std::ios::binary | std::ios::out);
    if (!outFile.is_open()) {
      std::cerr << "Could not open file for writing: " << filename << std::endl;
      return false;
    }

    // OPTIONAL: Write a header or “magic number” to ensure file validity
    // Example:
    const uint32_t magicNumber = sys.currentStep();
    outFile.write(reinterpret_cast<const char *>(&magicNumber),
                  sizeof(magicNumber));

    // Write number of records
    size_t recordCount = sys.particleNumber();
    outFile.write(reinterpret_cast<const char *>(&recordCount),
                  sizeof(recordCount));

    // Write each record in binary
    // NOTE: This assumes 'Record' is trivially copyable (POD).
    // If not, you'd serialize field by field.
    for (const Particle &p : sys.particles()) {
      int id = p.id();
      double mass = p.mass();
      outFile.write(reinterpret_cast<const char *>(&id), sizeof(int));
      outFile.write(reinterpret_cast<const char *>(&mass), sizeof(double));
      outFile.write(reinterpret_cast<const char *>(&p.coord()),
                    sizeof(Vector3<double>));
      outFile.write(reinterpret_cast<const char *>(&p.velocity()),
                    sizeof(Vector3<double>));
    }

    outFile.close();
    return true;
  }

  /**
   * @brief Вспомогательная функция для чтения записей в двоичном формате.
   *
   * Вспомогательная функция для чтения записей в двоичном формате.
   *
   * @param filename Имя файла для чтения.
   * @param sys Система, данные которой будут прочитаны.
   * @return true, если чтение прошло успешно, false в противном случае.
   */
  bool readRecordsFromFile(const std::string &filename, System &sys) {
    std::ifstream inFile(filename, std::ios::binary | std::ios::in);
    if (!inFile.is_open()) {
      throw std::invalid_argument("Could not open file for reading");
      return false;
    }

    // Read and validate magic number
    uint32_t magicNumber;
    inFile.read(reinterpret_cast<char *>(&magicNumber), sizeof(magicNumber));
    if (magicNumber != backup_restore_step_) {
      std::cerr << "Invalid backup file (bad restore step)." << std::endl;
      return false;
    }

    // Read number of records
    size_t recordCount = 0;
    inFile.read(reinterpret_cast<char *>(&recordCount), sizeof(recordCount));
    std::cout << "Reading " << recordCount << " records." << std::endl;
    sys.particles().clear();
    sys.particles().reserve(recordCount);

    // Get current position after reading header
    std::streampos headerEnd = inFile.tellg();

    // Calculate file size
    inFile.seekg(0, std::ios::end);
    std::streampos fileSize = inFile.tellg();

    // Reset to position after header
    inFile.seekg(headerEnd);

    // Calculate expected remaining size for old and new formats
    size_t recordSizeOldFormat = sizeof(double) + sizeof(Vector3<double>) * 2;
    size_t recordSizeNewFormat =
        sizeof(int) + sizeof(double) + sizeof(Vector3<double>) * 2;
    size_t expectedOldSize = recordCount * recordSizeOldFormat;
    size_t expectedNewSize = recordCount * recordSizeNewFormat;

    // Get actual remaining size
    std::streamoff remainingSize = fileSize - headerEnd;

    // Determine if we're using old format
    bool isOldFormat = (static_cast<size_t>(remainingSize) == expectedOldSize);

    std::cout << "Detected " << (isOldFormat ? "old" : "new")
              << " backup format." << std::endl;

    // Read in each Record
    for (size_t i = 0; i < recordCount; ++i) {
      int id = i;
      double mass;
      Vector3<double> coord;
      Vector3<double> velocity;

      if (!isOldFormat) {
        // New format has explicit ID
        inFile.read(reinterpret_cast<char *>(&id), sizeof(int));
      }
      inFile.read(reinterpret_cast<char *>(&mass), sizeof(double));
      inFile.read(reinterpret_cast<char *>(&coord), sizeof(Vector3<double>));
      inFile.read(reinterpret_cast<char *>(&velocity), sizeof(Vector3<double>));

      Particle p;
      p.setId(id);
      p.setMass(mass);
      p.setCoord(coord);
      p.setVelocity(velocity);
      sys.particles().push_back(std::move(p)); // Use move to avoid extra copy
    }
    inFile.close();
    return true;
  }

  /**
   * @brief Удаление старых резервных копий, оставляя только последние 5.
   *
   * Удаление старых резервных копий, оставляя только последние 5.
   *
   * @param step Шаг, на котором будет произведено удаление резервных копий.
   */
  void cleanupBackup(int step) const {
    if (step <= 5 * backup_frequency_) {
      return; // No need to delete backups yet
    }

    int old_step = step - 5 * backup_frequency_;
    // Preserve milestone backups (divisible by 100000 or 1000000)
    if (old_step % 100000 == 0 || old_step % 1000000 == 0) {
      return;
    }

    char filename[256]; // Ensure the buffer is large enough
    std::snprintf(filename, sizeof(filename), "%s/backup_%d.bak",
                  base_folder_.c_str(), old_step);

    std::filesystem::path backupPath(filename);
    if (std::filesystem::exists(backupPath)) {
      try {
        std::filesystem::remove(backupPath);
      } catch (const std::filesystem::filesystem_error &e) {
        std::cerr << "Error deleting backup: " << e.what() << std::endl;
      }
    }
  }
};

#endif // BACKUP_MANAGER_H
