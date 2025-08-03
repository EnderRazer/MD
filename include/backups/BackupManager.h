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

#include "classes/Particles.h"
#include "core/System.h"

using json = nlohmann::json;
// The BackupManager class handles saving/restoring these Records
class BackupManager {
public:
  inline bool enabled() const { return toggle_; }
  inline int frequency() const { return backup_frequency_; }

  // Write backup data to file
  // Returns true on success, false on failure
  bool createBackup(System &sys) const {
    std::ostringstream oss;
    oss << base_folder_ << "/"
        << "backup_" << sys.currentStep() << ".bak";
    cleanupBackup(sys.currentStep());
    return writeRecordsToFile(oss.str(), sys);
  }
  // Restore data from backup file
  // Returns true on success, false on failure
  bool restoreBackup(System &sys) {
    std::string input;
    std::cout << "Введите путь до файла бекапа: ";
    std::cin >> input;
    return readRecordsFromFile(
        input + "backup_" + std::to_string(backup_restore_step_) + ".bak", sys);
  };

  std::string getData() const {
    std::ostringstream oss;
    oss << "Backupdata"
        << "\n\tBackup enabled: " << toggle_
        << "\n\tBackup restore step: " << backup_restore_step_
        << "\n\tBackup frequency: " << backup_frequency_ << "\n";
    return oss.str();
  }

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
  ~BackupManager() = default;

private:
  std::string base_folder_{""};
  bool toggle_{false};
  int backup_frequency_{100};
  int backup_restore_step_{0};
  // Helper to write records in binary
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
    Particles &p = sys.particles();  
    for(int i = 0; i < p.size(); i++){
      outFile.write(reinterpret_cast<const char *>(&i), sizeof(int));
      outFile.write(reinterpret_cast<const char *>(&p.mass_[i]), sizeof(double));
      outFile.write(reinterpret_cast<const char *>(&p.coord_x_[i]),sizeof(double));
      outFile.write(reinterpret_cast<const char *>(&p.coord_y_[i]),sizeof(double));
      outFile.write(reinterpret_cast<const char *>(&p.coord_z_[i]),sizeof(double));
      outFile.write(reinterpret_cast<const char *>(&p.velocity_x_[i]),sizeof(double));
      outFile.write(reinterpret_cast<const char *>(&p.velocity_y_[i]),sizeof(double));
      outFile.write(reinterpret_cast<const char *>(&p.velocity_z_[i]),sizeof(double));
    }

    outFile.close();
    return true;
  }

  // Helper to read records in binary
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
    sys.particles().resize(recordCount);

    // Get current position after reading header
    std::streampos headerEnd = inFile.tellg();

    // Calculate file size
    inFile.seekg(0, std::ios::end);
    std::streampos fileSize = inFile.tellg();

    // Reset to position after header
    inFile.seekg(headerEnd);

    // Calculate expected remaining size for old and new formats
    //size_t recordSizeOldFormat = sizeof(double) + sizeof(Vector3<double>) * 2;
    //size_t recordSizeNewFormat =
    //    sizeof(int) + sizeof(double) + sizeof(Vector3<double>) * 2;
    //size_t expectedOldSize = recordCount * recordSizeOldFormat;
    //size_t expectedNewSize = recordCount * recordSizeNewFormat;

    // Get actual remaining size
    //std::streamoff remainingSize = fileSize - headerEnd;

    // Determine if we're using old format
    //bool isOldFormat = (static_cast<size_t>(remainingSize) == expectedOldSize);

    //std::cout << "Detected " << (isOldFormat ? "old" : "new")
    //          << " backup format." << std::endl;

    // Read in each Record
    Particles &particles = sys.particles();
    for (size_t i = 0; i < recordCount; ++i) {
      int id = i;
      double mass;
      double coord_x,coord_y,coord_z;
      double velocity_x,velocity_y,velocity_z;

      
      inFile.read(reinterpret_cast<char *>(&id), sizeof(int));
      inFile.read(reinterpret_cast<char *>(&mass), sizeof(double));
      inFile.read(reinterpret_cast<char *>(&coord_x), sizeof(double));
      inFile.read(reinterpret_cast<char *>(&coord_y), sizeof(double));
      inFile.read(reinterpret_cast<char *>(&coord_z), sizeof(double));
      inFile.read(reinterpret_cast<char *>(&velocity_x), sizeof(double));
      inFile.read(reinterpret_cast<char *>(&velocity_y), sizeof(double));
      inFile.read(reinterpret_cast<char *>(&velocity_z), sizeof(double));

      particles.mass_[i] = mass;
      
      particles.coord_x_[i] = coord_x;
      particles.coord_y_[i] = coord_y;
      particles.coord_z_[i] = coord_z;

      particles.velocity_x_[i] = velocity_x;
      particles.velocity_y_[i] = velocity_y;
      particles.velocity_z_[i] = velocity_z;
    }
    inFile.close();
    return true;
  }

  // Delete old backups, keeping only the most recent 5
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
