class Barostat {
public:
  enum class BarostatType { BERENDSEN };
  virtual ~Barostat() = default;
  virtual void applyPressureControl(System &sys) = 0;
  virtual std::string getData() const = 0;

  inline virtual const BarostatType getBarostatType() const = 0;
};
