#ifndef BACKUP_MANAGER_H
#define BACKUP_MANAGER_H
#include <filesystem>
using json = nlohmann::json;
// The BackupManager class handles saving/restoring these Records
class BackupManager
{
public:
  inline bool enabled() const { return toggle_; }
  inline int frequency() const { return backup_frequency_; }
  // Write backup data to file
  // Returns true on success, false on failure
  bool createBackup(System &sys) const
  {
    std::ostringstream oss;
    oss << base_folder_ << "/" << "backup_" << sys.currentStep() << ".bak";
    return writeRecordsToFile(oss.str(), sys);
  }
  // Restore data from backup file
  // Returns true on success, false on failure
  bool restoreBackup(System &sys)
  {
    std::string input;
    std::cout << "Введите путь до файла бекапа: ";
    std::cin >> input;
    return readRecordsFromFile(
        input + "backup_" + std::to_string(backup_restore_step_) + ".bak",
        sys);
  };

  std::string getData() const
  {
    std::ostringstream oss;
    oss << "Backupdata" << "\n\tBackup enabled: " << toggle_
        << "\n\tBackup restore step: " << backup_restore_step_
        << "\n\tBackup frequency: " << backup_frequency_ << "\n";
    return oss.str();
  }

  BackupManager(json config, Settings &settings)
  {
    base_folder_ = "data_" + std::to_string(settings.seed()) + "/backups";
    if (!std::filesystem::exists(base_folder_))
    {
      if (!std::filesystem::create_directories(base_folder_))
      {
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
  bool writeRecordsToFile(const std::string &filename, System &sys) const
  {
    // Open file in binary mode, std::ios::out ensures we overwrite
    std::ofstream outFile(filename, std::ios::binary | std::ios::out);
    if (!outFile.is_open())
    {
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
    for (const Particle &p : sys.particles())
    {
      double mass = p.getMass();
      outFile.write(reinterpret_cast<const char *>(&mass), sizeof(double));
      outFile.write(reinterpret_cast<const char *>(&p.coord()),
                    sizeof(Vector3<double>));
      outFile.write(reinterpret_cast<const char *>(&p.velocity()),
                    sizeof(Vector3<double>));
    }

    outFile.close();
    return true;
  }

  // Helper to read records in binary
  bool readRecordsFromFile(const std::string &filename, System &sys)
  {
    std::ifstream inFile(filename, std::ios::binary | std::ios::in);
    if (!inFile.is_open())
    {
      throw std::invalid_argument("Could not open file for reading");
      return false;
    }

    // Read and validate magic number
    uint32_t magicNumber;
    inFile.read(reinterpret_cast<char *>(&magicNumber), sizeof(magicNumber));
    if (magicNumber != backup_restore_step_)
    {
      std::cerr << "Invalid backup file (bad restore step)." << std::endl;
      return false;
    }

    // Read number of records
    size_t recordCount = 0;
    inFile.read(reinterpret_cast<char *>(&recordCount), sizeof(recordCount));
    std::cout << "Reading " << recordCount << " records." << std::endl;
    sys.particles().clear();
    sys.particles().reserve(recordCount);
    // Read in each Record
    for (size_t i = 0; i < recordCount; ++i)
    {
      double mass;
      Vector3<double> coord;
      Vector3<double> velocity;

      inFile.read(reinterpret_cast<char *>(&mass), sizeof(double));
      inFile.read(reinterpret_cast<char *>(&coord), sizeof(Vector3<double>));
      inFile.read(reinterpret_cast<char *>(&velocity), sizeof(Vector3<double>));

      Particle p;
      p.setMass(mass);
      p.setCoord(coord);
      p.setVelocity(velocity);
      sys.particles().push_back(std::move(p)); // Use move to avoid extra copy
    }
    inFile.close();
    return true;
  }
};

#endif // BACKUP_MANAGER_H