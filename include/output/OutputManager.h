#ifndef OUTPUT_MANAGER_H
#define OUTPUT_MANAGER_H

namespace fs = std::filesystem;
using json = nlohmann::json;

class OutputManager
{
private:
  std::string output_dir_;
  std::string step_dir_;
  std::string ens_dir_;

  int frequency_;

public:
  OutputManager(const json &config, Settings &settings)
  {
    output_dir_ = config.value("main_directory", "data") + "_" + std::to_string(settings.seed());
    if (!fs::exists(output_dir_))
    {
      fs::create_directory(output_dir_);
    }
    step_dir_ = output_dir_ + "/" + config.value("steps_directory", "steps_data");
    if (!fs::exists(step_dir_))
    {
      fs::create_directory(step_dir_);
    }
    ens_dir_ = output_dir_ + "/" + config.value("ensembles_directory", "ensembles_data");
    if (!fs::exists(ens_dir_))
    {
      fs::create_directory(ens_dir_);
    }
    frequency_ = config.value("frequency", 1);
  }

  // Запрещаем копирование
  OutputManager(const OutputManager &) = delete;
  OutputManager &operator=(const OutputManager &) = delete;

  inline const int frequency() const { return frequency_; }

  void writeModellingProperties(json &config) const
  {
    std::string filename = output_dir_ + "/launch_config.json";
    std::ofstream file(filename, std::ios::app);
    file << config.dump(4) << std::endl;
    file.close();
  }
  // Write accumulated properties (temperature, pressure, energies, pulse) over
  // steps
  void writeSystemProperties(System &sys) const
  {
    std::string filename = output_dir_ + "/properties.csv";
    bool file_exists = fs::exists(filename);

    std::ofstream file(filename, std::ios::app);
    file << std::fixed
         << std::setprecision(16); // Fixed-point notation with high precision

    // Write header if file does not exist
    if (!file_exists)
    {
      file << "Step;" << "Temperature;"
           << "Pressure;" << "TempAvg;"
           << "PressAvg;" << "KineticE;"
           << "PotentialE;" << "TotalE;"
           << "Pulse;"
           << "\n";
    }
    // Write formatted data
    file << sys.currentStep() << ";"
         << sys.temperature() << ";"
         << sys.pressure() << ";"
         << sys.temperatureAvg() << ";"
         << sys.pressureAvg() << ";"
         << sys.energies().get(Energy::EnergyType::Kinetic) << ";"
         << sys.energies().get(Energy::EnergyType::Potential) << ";"
         << sys.energies().get(Energy::EnergyType::Full) << ";"
         << sys.pulse() << ";"
         << "\n";
    file.close();
  }

  // Create a directory for each step and write particle data
  void writeStepData(System &sys) const
  {
    writeParticlesDetailed(sys, step_dir_ + "/step_" +
                                    std::to_string(sys.currentStep()) +
                                    "_particles.csv");
  }

  void writeEnsemblesDetailed(const EnsembleManager &manager) const
  {
    std::vector<Ensemble> ensembles = manager.getEnsembles();

    for (Ensemble &ens : ensembles)
    {
      if (ens.cur_ - 1 < 0)
        continue;
      std::string filename =
          ens_dir_ + "/ensemble_" + std::to_string(ens.id_) + ".csv";
      bool file_exists = fs::exists(filename);
      std::ofstream file(filename, std::ios::app);
      file << std::fixed
           << std::setprecision(16); // Fixed-point notation with high precision

      // Header
      if (!file_exists)
        file << "n;" << "ACFV;" << "CDiff;" << "ACFP;" << "CVisc;" << "\n";

      file << ens.cur_ - 1 << ";"
           << ens.accum_acfv_[ens.cur_ - 1] << ";"
           << ens.coef_diff_integral_[ens.cur_ - 1] << ";"
           << ens.accum_acfp_[ens.cur_ - 1] << ";"
           << ens.coef_visc_integral_[ens.cur_ - 1] << ";"
           << "\n";
      file.close();
    }
  }

  void writeAvgEnsembleDetailed(const EnsembleManager &manager) const
  {
    std::string filename = ens_dir_ + "/ensemble_avg.csv";
    bool file_exists = fs::exists(filename);
    std::ofstream file(filename);
    file << std::fixed
         << std::setprecision(16); // Fixed-point notation with high precision

    // Header
    if (!file_exists)
      file << "n;" << "ACFV;"
           << "CDiff;" << "ACFP;"
           << "CVisc;"
           << "\n";

    std::vector<double> ensembles_accum_acfv_ = manager.getEnsemblesAccumAcfv();
    std::vector<double> ensembles_accum_acfp_ = manager.getEnsemblesAccumAcfp();

    std::vector<double> ensembles_coef_diff_integral_ =
        manager.getEnsemblesCoefDiffIntegral();
    std::vector<double> ensembles_coef_visc_integral_ =
        manager.getEnsemblesCoefViscIntegral();

    for (int i = 0; i < manager.getEnsembleSize(); i++)
    {
      file << i << ";"
           << ensembles_accum_acfv_[i] << ";"
           << ensembles_coef_diff_integral_[i] << ";"
           << ensembles_accum_acfp_[i] << ";"
           << ensembles_coef_visc_integral_[i] << ";"
           << "\n";
    }
    file.close();
  }

  const std::string getData() const
  {
    std::ostringstream oss;
    oss << "Output manager data:" << "\n\tMain directory: " << output_dir_
        << "\n\tSteps directory" << step_dir_
        << "\n\tEnsembles directory: " << ens_dir_ << std::endl;
    return oss.str();
  }

private:
  void writeParticlesDetailed(System &sys, const std::string &filename) const
  {
    std::ofstream file(filename);
    file << std::fixed
         << std::setprecision(16); // Fixed-point notation with high precision

    // Header
    file << "n;"
         << "Cx;" << "Cy;" << "Cz;"
         << "Vx;" << "Vy;" << "Vz;"
         << "Fx;" << "Fy;" << "Fz;"
         << "PotentialE;" << "KineticE;" << "TotalE;"
         << "\n";
    int id = 0;
    for (const auto &particle : sys.particles())
    {
      const auto &coord = particle.coord();
      const auto &velocity = particle.velocity();
      const auto &force = particle.force();
      const auto &energies = particle.energies();
      file << id++ << ";"
           << coord.x() << ";"
           << coord.y() << ";"
           << coord.z() << ";"
           << velocity.x() << ";"
           << velocity.y() << ";"
           << velocity.z() << ";"
           << force.x() << ";"
           << force.y() << ";"
           << force.z() << ";"
           << energies.get(Energy::EnergyType::Potential) << ";"
           << energies.get(Energy::EnergyType::Kinetic) << ";"
           << energies.get(Energy::EnergyType::Full) << ";"
           << "\n";
    }
    file.close();
  }
};

#endif // OUTPUT_MANAGER_H