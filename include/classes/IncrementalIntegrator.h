#ifndef INCREMENTAL_INTEGRATOR_H
#define INCREMENTAL_INTEGRATOR_H
class IncrementalIntegrator {
private:
  double accumulated_integral_ = 0.0;
  double h_;
  double last_value_ = 0.0;

public:
  explicit IncrementalIntegrator(double step_size) : h_(step_size) {}

  double add_partial(const double value) {
    double local_integral = (last_value_ + value) * h_ / 2.0;

    // Update state
    accumulated_integral_ += local_integral;
    last_value_ = value; // Store last function value

    return accumulated_integral_;
  }

  double get_total_integral() const { return accumulated_integral_; }
};
#endif