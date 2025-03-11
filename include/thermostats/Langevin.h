//===================================================================================
// Реализация термостата Ланжевена
//===================================================================================
#ifndef THERMOSTAT_LNGVN
#define THERMOSTAT_LNGVN
using json = nlohmann::json;

class ThermostatLangevine : public Thermostat {
public:
  const ThermostatType type_{ThermostatType::LANGEVIN};

private:
  int seed_{0};
  bool toggle_{false};
  double tau_t_{0.0};
  double pref_temperature_{0.0};
  Vector3<double> Vs_{}; // Скорость потока

  double sigma_lang_{0.0}; // Константа для термостата Ланжевена

  double eps_{0.0};                   // Случайная величина (0;1)
  std::vector<Vector3<double>> Fr_{}; // Случайная сила
  std::vector<Vector3<double>> Ft_{}; // Сила трения
public:
  ~ThermostatLangevine() noexcept override = default;
  ThermostatLangevine(json &config, Settings &settings, int particle_number) {
    toggle_ = config.value("toggle", false);
    pref_temperature_ = config.value("pref_temp", 0.0);
    tau_t_ = config.value("tau", 0.0);
    Vs_ = config["stream_velocity"].get<std::vector<double>>();
    sigma_lang_ = sqrt((2 * settings.constants().KBOLTZMAN * pref_temperature_ *
                        settings.mass()) /
                       (tau_t_ * settings.dt()));
    seed_ = settings.seed();
    Fr_.reserve(particle_number);
    Ft_.reserve(particle_number);
  }

  // Запрещаем копирование
  ThermostatLangevine(const ThermostatLangevine &) = delete;
  ThermostatLangevine &operator=(const ThermostatLangevine &) = delete;

  void LGVN_generateFr(System &sys) {
    int pn = sys.particleNumber();
    for (int i = 0; i < pn / 2; i++) {

      eps_ = gasdev(seed_);
      Fr_[i].x() = sigma_lang_ * eps_;
      Fr_[pn / 2 + i].x() = -Fr_[i].x();

      eps_ = gasdev(seed_);
      Fr_[i].y() = sigma_lang_ * eps_;
      Fr_[pn / 2 + i].y() = -Fr_[i].y();

      eps_ = gasdev(seed_);
      Fr_[i].z() = sigma_lang_ * eps_;
      Fr_[pn / 2 + i].z() = -Fr_[i].z();
    }
  };

  void LGVN_generateFt(System &sys) {
    int pn = sys.particleNumber();
    std::vector<Particle> &particles = sys.particles();
    for (int i = 0; i < pn; i++) {
      Ft_[i] =
          -particles[i].getMass() * (particles[i].velocity() - Vs_) / tau_t_;
    }
  };
  void applyTemperatureControl(System &sys) override {
    LGVN_generateFr(sys);
    LGVN_generateFt(sys);
    int pn = sys.particleNumber();
    std::vector<Particle> &particles = sys.particles();
    for (int i = 0; i < pn; i++) {
      particles[i].force() += 2 * (Fr_[i] + Ft_[i]);
    }
  };

  std::string getData() const override {
    std::ostringstream oss;
    oss.precision(16);
    oss << "Using Langevin thermostat!" << "\n\tToggled: " << toggle_
        << "\n\tPrefered temperature: " << pref_temperature_
        << "\n\tTau_t: " << tau_t_ << "\n\tStream velocity: " << Vs_
        << "\n\tSIGMA_LANG: " << sigma_lang_ << "\n";
    return oss.str();
  }
  inline const ThermostatType getThermostatType() const override {
    return type_;
  }
};
#endif