#ifndef OUTPUT_MANAGER_H
#define OUTPUT_MANAGER_H

namespace fs = std::filesystem;
using json = nlohmann::json;

class OutputManager {
private:
  std::string output_dir_;
  std::string step_dir_;
  std::string ens_dir_;

public:
  OutputManager(const std::string &output_dir, Settings &settings)
      : output_dir_(output_dir + "_" + std::to_string(settings.seed())) {
    if (!fs::exists(output_dir_)) {
      fs::create_directory(output_dir_);
    }
    step_dir_ = output_dir_ + "/steps";
    if (!fs::exists(step_dir_)) {
      fs::create_directory(step_dir_);
    }
    ens_dir_ = output_dir_ + "/ensembles";
    if (!fs::exists(ens_dir_)) {
      fs::create_directory(ens_dir_);
    }
  }

  // Запрещаем копирование
  OutputManager(const OutputManager &) = delete;
  OutputManager &operator=(const OutputManager &) = delete;

  void writeModellingProperties(json &config) const {
    std::string filename = output_dir_ + "/launch_config.json";
    std::ofstream file(filename, std::ios::app);
    file << config.dump(4) << std::endl;
    file.close();
  }
  // Write accumulated properties (temperature, pressure, energies, pulse) over
  // steps
  void writeSystemProperties(System &sys) const {
    std::string filename = output_dir_ + "/properties.dat";
    bool file_exists = fs::exists(filename);

    std::ofstream file(filename, std::ios::app);
    file << std::fixed
         << std::setprecision(16); // Fixed-point notation with high precision
    // Define column width
    constexpr int width = 24;
    constexpr int id_width = 8;
    // Write header if file does not exist
    if (!file_exists) {
      file << std::setw(id_width) << "Step" << std::setw(width) << "Temperature"
           << std::setw(width) << "Pressure" << std::setw(width) << "TempAvg"
           << std::setw(width) << "PressAvg" << std::setw(width) << "KineticE"
           << std::setw(width) << "PotentialE" << std::setw(width) << "TotalE"
           << std::setw(width) << "Pulse"
           << "\n";
    }
    // Write formatted data
    file << std::setw(id_width) << sys.currentStep() << std::setw(width)
         << sys.temperature() << std::setw(width) << sys.pressure()
         << std::setw(width) << sys.temperatureAvg() << std::setw(width)
         << sys.pressureAvg() << std::setw(width)
         << sys.energies().get(Energy::EnergyType::Kinetic) << std::setw(width)
         << sys.energies().get(Energy::EnergyType::Potential)
         << std::setw(width) << sys.energies().get(Energy::EnergyType::Full)
         << std::setw(width) << sys.pulse() << "\n";
    file.close();
  }

  // Create a directory for each step and write particle data
  void writeStepData(System &sys) const {
    writeParticlesDetailed(sys, step_dir_ + "/step_" +
                                    std::to_string(sys.currentStep()) +
                                    "_particles.dat");
  }

  void writeEnsemblesDetailed(const EnsembleManager &manager) const {
    std::vector<Ensemble> ensembles = manager.getEnsembles();

    for (int i = 0; i < ensembles.size(); i++) {
      if (ensembles[i].cur_ - 1 < 0)
        return;
      std::string filename =
          ens_dir_ + "/ensemble_" + std::to_string(ensembles[i].id_);
      bool file_exists = fs::exists(filename);
      std::ofstream file(filename, std::ios::app);
      file << std::fixed
           << std::setprecision(16); // Fixed-point notation with high precision
      // Define column width
      constexpr int width = 24;
      constexpr int id_width = 6;

      // Header
      if (!file_exists)
        file << std::setw(id_width) << "n" << std::setw(width) << "ACFV"
             << std::setw(width) << "CDiff" << std::setw(width) << "ACFP"
             << std::setw(width) << "CVisc"
             << "\n";

      file << std::setw(id_width) << ensembles[i].cur_ - 1 << std::setw(width)
           << ensembles[i].accum_acfv_[ensembles[i].cur_ - 1]
           << std::setw(width)
           << ensembles[i].coef_diff_integral_[ensembles[i].cur_ - 1]
           << std::setw(width)
           << ensembles[i].accum_acfp_[ensembles[i].cur_ - 1]
           << std::setw(width)
           << ensembles[i].coef_visc_integral_[ensembles[i].cur_ - 1] << "\n";
      file.close();
    }
  }

  void writeAvgEnsembleDetailed(const EnsembleManager &manager) const {
    std::string filename = ens_dir_ + "/ensemble_avg";
    bool file_exists = fs::exists(filename);
    std::ofstream file(filename);
    file << std::fixed
         << std::setprecision(16); // Fixed-point notation with high precision
    // Define column width
    constexpr int width = 24;
    constexpr int id_width = 6;

    // Header
    if (!file_exists)
      file << std::setw(id_width) << "n" << std::setw(width) << "ACFV"
           << std::setw(width) << "CDiff" << std::setw(width) << "ACFP"
           << std::setw(width) << "CVisc"
           << "\n";

    std::vector<double> ensembles_accum_acfv_ = manager.getEnsemblesAccumAcfv();
    std::vector<double> ensembles_accum_acfp_ = manager.getEnsemblesAccumAcfp();

    std::vector<double> ensembles_coef_diff_integral_ =
        manager.getEnsemblesCoefDiffIntegral();
    std::vector<double> ensembles_coef_visc_integral_ =
        manager.getEnsemblesCoefViscIntegral();

    for (int i = 0; i < manager.getEnsembleSize(); i++) {
      file << std::setw(id_width) << i << std::setw(width)
           << ensembles_accum_acfv_[i] << std::setw(width)
           << ensembles_coef_diff_integral_[i] << std::setw(width)
           << ensembles_accum_acfp_[i] << std::setw(width)
           << ensembles_coef_visc_integral_[i] << "\n";
    }
    file.close();
  }

  const std::string getData() const {
    std::ostringstream oss;
    oss << "Output manager data:" << "\n\tMain directory: " << output_dir_
        << "\n\tSteps directory" << step_dir_
        << "\n\tEnsembles directory: " << ens_dir_ << std::endl;
    return oss.str();
  }

private:
  void writeParticlesDetailed(System &sys, const std::string &filename) const {
    std::ofstream file(filename);
    file << std::fixed
         << std::setprecision(16); // Fixed-point notation with high precision
    // Define column width
    constexpr int width = 24;
    constexpr int id_width = 6;

    // Header
    file << std::setw(id_width) << "n" << std::setw(width) << "Cx"
         << std::setw(width) << "Cy" << std::setw(width) << "Cz"
         << std::setw(width) << "Vx" << std::setw(width) << "Vy"
         << std::setw(width) << "Vz" << std::setw(width) << "Fx"
         << std::setw(width) << "Fy" << std::setw(width) << "Fz"
         << std::setw(width) << "PotentialE" << std::setw(width) << "KineticE"
         << std::setw(width) << "TotalE"
         << "\n";
    int id = 0;
    for (const auto &particle : sys.particles()) {
      const auto &coord = particle.coord();
      const auto &velocity = particle.velocity();
      const auto &force = particle.force();
      const auto &energies = particle.energies();
      file << std::setw(id_width) << id++ << std::setw(width) << coord.x()
           << std::setw(width) << coord.y() << std::setw(width) << coord.z()
           << std::setw(width) << velocity.x() << std::setw(width)
           << velocity.y() << std::setw(width) << velocity.z()
           << std::setw(width) << force.x() << std::setw(width) << force.y()
           << std::setw(width) << force.z() << std::setw(width)
           << energies.get(Energy::EnergyType::Potential) << std::setw(width)
           << energies.get(Energy::EnergyType::Kinetic) << std::setw(width)
           << energies.get(Energy::EnergyType::Full) << "\n";
    }
    file.close();
  }
};

#endif // OUTPUT_MANAGER_H