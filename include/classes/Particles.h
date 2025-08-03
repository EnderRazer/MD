#ifndef PARTICLE_H
#define PARTICLE_H

#include <cstddef>
#include <iostream>
#include <vector>

#include "AlignedAllocator.h"

using AlignedVector = std::vector<double, AlignedAllocator<double, 64>>;

struct ForceCalcValues {
  int interaction_count{0};
  double rVec_x{0}, rVec_y{0}, rVec_z{0};
  
  double force_x{0},force_y{0},force_z{0};

  double e_pot{0.0};

  double virial_xx{0},virial_xy{0},virial_xz{0};
  double virial_yx{0},virial_yy{0},virial_yz{0};
  double virial_zx{0},virial_zy{0},virial_zz{0};
};

class Particles {
public:
  size_t N{};
  alignas(64) std::vector<double, AlignedAllocator<double, 64>>  mass_{}; // Particle mass

  // Position (coordinate)
  alignas(64) std::vector<double, AlignedAllocator<double, 64>> coord_x_, coord_y_, coord_z_;   

  // Velocity
  alignas(64) std::vector<double, AlignedAllocator<double, 64>>  velocity_x_{}, velocity_y_{}, velocity_z_{}; 

  // Force acting on the particle
  alignas(64) std::vector<double, AlignedAllocator<double, 64>>  force_x_{}, force_y_{}, force_z_{};    

  // Virials (3x3 matrix)
  alignas(64) std::vector<double, AlignedAllocator<double, 64>>  virial_xx_{}, virial_xy_{}, virial_xz_{};
  alignas(64) std::vector<double, AlignedAllocator<double, 64>>  virial_yx_{}, virial_yy_{}, virial_yz_{};
  alignas(64) std::vector<double, AlignedAllocator<double, 64>>  virial_zx_{}, virial_zy_{}, virial_zz_{};

  // Energies
  alignas(64) std::vector<double, AlignedAllocator<double, 64>>  e_pot_{},e_kin_{},e_int_{},e_term_{},e_full_{};

  // Pulse of the particle 
  alignas(64) std::vector<double, AlignedAllocator<double, 64>>  pulse_{}; 

  // EAM Properties
  // Electron density of the particle
  alignas(64) std::vector<double, AlignedAllocator<double, 64>>  electron_density_{}; 

  // Pair potential energy of the particle
  alignas(64) std::vector<double, AlignedAllocator<double, 64>>  pair_potential_{};   

  // OM of particles
  alignas(64) std::vector<double, AlignedAllocator<double, 64>>  ce_{};
  alignas(64) std::vector<double, AlignedAllocator<double, 64>>  d_ce_{};
public:
  // Default constructor: Initializes everything to zero.
  Particles() = default;
  ~Particles() = default;

  Particles(size_t n) : N(n) {
    mass_.resize(N);
    coord_x_.resize(N); coord_y_.resize(N); coord_z_.resize(N);
    velocity_x_.resize(N); velocity_y_.resize(N); velocity_z_.resize(N);
    force_x_.resize(N); force_y_.resize(N); force_z_.resize(N);

    virial_xx_.resize(N); virial_xy_.resize(N); virial_xz_.resize(N);
    virial_yx_.resize(N); virial_yy_.resize(N); virial_yz_.resize(N);
    virial_zx_.resize(N); virial_zy_.resize(N); virial_zz_.resize(N);

    e_kin_.resize(N), e_pot_.resize(N), e_int_.resize(N), e_term_.resize(N), e_full_.resize(N);
    
    pulse_.resize(N);

    electron_density_.resize(N);
    pair_potential_.resize(N);
    ce_.resize(N);
    d_ce_.resize(N);
  };

  void resize(size_t n) {
    std::cout << "Calling resize with n = " << n << std::endl;
    N = n;
    mass_.resize(n);
    coord_x_.resize(n); coord_y_.resize(n); coord_z_.resize(n);
    velocity_x_.resize(n); velocity_y_.resize(n); velocity_z_.resize(n);
    force_x_.resize(n); force_y_.resize(n); force_z_.resize(n);

    virial_xx_.resize(n); virial_xy_.resize(n); virial_xz_.resize(n);
    virial_yx_.resize(n); virial_yy_.resize(n); virial_yz_.resize(n);
    virial_zx_.resize(n); virial_zy_.resize(n); virial_zz_.resize(n);

    e_kin_.resize(n), e_pot_.resize(n), e_int_.resize(n), e_term_.resize(n), e_full_.resize(n);
    pulse_.resize(n);

    electron_density_.resize(n);
    pair_potential_.resize(n);
    ce_.resize(N);
    d_ce_.resize(N);
  }

  inline void applyForceInteraction(const size_t index, const ForceCalcValues &result) {
    force_x_[index] = result.force_x;
    force_y_[index] = result.force_y;
    force_z_[index] = result.force_z;

    virial_xx_[index] = result.virial_xx;
    virial_xy_[index] = result.virial_xy;
    virial_xz_[index] = result.virial_xz;

    virial_yx_[index] = result.virial_yx;
    virial_yy_[index] = result.virial_yy;
    virial_yz_[index] = result.virial_yz;

    virial_zx_[index] = result.virial_zx;
    virial_zy_[index] = result.virial_zy;
    virial_zz_[index] = result.virial_zz;

    e_pot_[index] = result.e_pot;
  }

  inline void updateEnergy(double &vcm_x,double &vcm_y, double vcm_z) {
    double halfMass = 0.0;
    double velocitySquared = 0.0;
    double relativeVelocitySquared = 0.0;


    double e_kin = 0.0;
    double e_term = 0.0;

    for(int i=0; i<N;i++){
      halfMass = mass_[i]*0.5;
      velocitySquared = (velocity_x_[i]*velocity_x_[i]) + (velocity_y_[i]*velocity_y_[i]) + (velocity_z_[i]*velocity_z_[i]);
      relativeVelocitySquared = ((velocity_x_[i]-vcm_x)*(velocity_x_[i]-vcm_x)) 
                                  + ((velocity_y_[i]-vcm_y)*(velocity_y_[i]-vcm_y)) 
                                  + ((velocity_z_[i]-vcm_z)*(velocity_z_[i]-vcm_z));
      e_kin = halfMass*velocitySquared;
      e_term = halfMass*relativeVelocitySquared;
      
      e_kin_[i] = e_kin;
      e_term_[i] = e_term;
      e_int_[i] = e_pot_[i]+e_term;
      e_full_[i] = e_pot_[i]+e_kin;
    }
  }

  inline void updatePulse() {
    double velocitySquared = 0.0;
    for(int i=0; i<N;i++){
      velocitySquared = (velocity_x_[i]*velocity_x_[i]) + (velocity_y_[i]*velocity_y_[i]) + (velocity_z_[i]*velocity_z_[i]);
      pulse_[i] = mass_[i]*velocitySquared;
    }
  }

  void computeCoordinatesSIMD(double dt) {
    double* __restrict x = coord_x_.data();
    double* __restrict y = coord_y_.data();
    double* __restrict z = coord_z_.data();
    const double* __restrict vx = velocity_x_.data();
    const double* __restrict vy = velocity_y_.data();
    const double* __restrict vz = velocity_z_.data();

  #pragma clang loop vectorize(enable)
    for (size_t i = 0; i < N; ++i) {
      x[i] += vx[i] * dt;
      y[i] += vy[i] * dt;
      z[i] += vz[i] * dt;
    }
  }

  void computeVelocitiesSIMD(double mt) {
    double* __restrict vx = velocity_x_.data();
    double* __restrict vy = velocity_y_.data();
    double* __restrict vz = velocity_z_.data();

    const double* __restrict fx = velocity_x_.data();
    const double* __restrict fy = velocity_y_.data();
    const double* __restrict fz = velocity_z_.data();


  #pragma clang loop vectorize(enable)
    for (size_t i = 0; i < N; ++i) {
      vx[i] += fx[i]*mt;
      vy[i] += fy[i]*mt;
      vz[i] += fz[i]*mt;
    }
  }

  // Вспомогательная функция Morton
  inline uint64_t expandBits(uint32_t v) {
      uint64_t x = v & 0x1fffff;
      x = (x | x << 32) & 0x1f00000000ffff;
      x = (x | x << 16) & 0x1f0000ff0000ff;
      x = (x | x << 8)  & 0x100f00f00f00f00f;
      x = (x | x << 4)  & 0x10c30c30c30c30c3;
      x = (x | x << 2)  & 0x1249249249249249;
      return x;
  }
  
  inline uint64_t morton3D(uint32_t x, uint32_t y, uint32_t z) {
      return expandBits(x) | (expandBits(y) << 1) | (expandBits(z) << 2);
  }

  void reorderByMorton() {
    if (N == 0) return;

    // 1. Генерируем ключи Morton
    std::vector<uint64_t> keys(N);
    std::vector<size_t> indices(N);
    for (size_t i = 0; i < N; ++i) indices[i] = i;

    // Находим границы для нормализации координат
    double minX = *std::min_element(coord_x_.begin(), coord_x_.end());
    double minY = *std::min_element(coord_y_.begin(), coord_y_.end());
    double minZ = *std::min_element(coord_z_.begin(), coord_z_.end());
    double maxX = *std::max_element(coord_x_.begin(), coord_x_.end());
    double maxY = *std::max_element(coord_y_.begin(), coord_y_.end());
    double maxZ = *std::max_element(coord_z_.begin(), coord_z_.end());

    double scaleX = (maxX - minX > 1e-12) ? (1023.0 / (maxX - minX)) : 1.0;
    double scaleY = (maxY - minY > 1e-12) ? (1023.0 / (maxY - minY)) : 1.0;
    double scaleZ = (maxZ - minZ > 1e-12) ? (1023.0 / (maxZ - minZ)) : 1.0;

    for (size_t i = 0; i < N; ++i) {
        uint32_t xi = (uint32_t)((coord_x_[i] - minX) * scaleX);
        uint32_t yi = (uint32_t)((coord_y_[i] - minY) * scaleY);
        uint32_t zi = (uint32_t)((coord_z_[i] - minZ) * scaleZ);
        keys[i] = morton3D(xi, yi, zi);
    }

    // 2. Сортируем индексы по ключам
    std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
        return keys[a] < keys[b];
    });

    // 3. Перепаковываем все массивы
    auto reorderVec = [&](auto &vec) {
        using T = typename std::decay<decltype(vec)>::type::value_type;
        std::vector<T, AlignedAllocator<T, 64>> tmp(N);
        for (size_t i = 0; i < N; ++i) tmp[i] = vec[indices[i]];
        vec.swap(tmp);
    };

    reorderVec(mass_);
    reorderVec(coord_x_); reorderVec(coord_y_); reorderVec(coord_z_);
    reorderVec(velocity_x_); reorderVec(velocity_y_); reorderVec(velocity_z_);
    reorderVec(force_x_); reorderVec(force_y_); reorderVec(force_z_);
    reorderVec(virial_xx_); reorderVec(virial_xy_); reorderVec(virial_xz_);
    reorderVec(virial_yx_); reorderVec(virial_yy_); reorderVec(virial_yz_);
    reorderVec(virial_zx_); reorderVec(virial_zy_); reorderVec(virial_zz_);
    reorderVec(e_pot_); reorderVec(e_kin_); reorderVec(e_int_);
    reorderVec(e_term_); reorderVec(e_full_);
    reorderVec(pulse_);
    reorderVec(electron_density_);
    reorderVec(pair_potential_);
  }
  //Геттеры
  constexpr size_t size() const { return N; }
  /*
  //Масса
  inline double& mass(size_t i) { return mass_[i]; }
  inline const double& mass(size_t i) const { return mass_[i]; }
  //Координаты
  inline double& coordX(size_t i) { return coord_x_[i]; }
  inline const double& coordX(size_t i) const { return coord_x_[i]; }
  inline double& coordY(size_t i) { return coord_y_[i]; }
  inline const double& coordY(size_t i) const { return coord_y_[i]; }
  inline double& coordZ(size_t i) { return coord_z_[i]; }
  inline const double& coordZ(size_t i) const { return coord_z_[i]; }
  //Скорости
  inline double& velocityX(size_t i) { return velocity_x_[i]; }
  inline const double& velocityX(size_t i) const { return velocity_x_[i]; }
  inline double& velocityY(size_t i) { return velocity_y_[i]; }
  inline const double& velocityY(size_t i) const { return velocity_y_[i]; }
  inline double& velocityZ(size_t i) { return velocity_z_[i]; }
  inline const double& velocityZ(size_t i) const { return velocity_z_[i]; }
  //Силы
  inline double& forceX(size_t i) { return force_x_[i]; }
  inline const double& forceX(size_t i) const { return force_x_[i]; }
  inline double& forceY(size_t i) { return force_y_[i]; }
  inline const double& forceY(size_t i) const { return force_y_[i]; }
  inline double& forceZ(size_t i) { return force_z_[i]; }
  inline const double& forceZ(size_t i) const { return force_z_[i]; }

  //Вириалы
  inline double& virialXX(size_t i) { return virial_xx_[i]; }
  inline const double& virialXX(size_t i) const { return virial_xx_[i]; }
  inline double& virialXY(size_t i) { return virial_xy_[i]; }
  inline const double& virialXY(size_t i) const { return virial_xy_[i]; }
  inline double& virialXZ(size_t i) { return virial_xz_[i]; }
  inline const double& virialXZ(size_t i) const { return virial_xz_[i]; }

  inline double& virialYX(size_t i) { return virial_yx_[i]; }
  inline const double& virialYX(size_t i) const { return virial_yx_[i]; }
  inline double& virialYY(size_t i) { return virial_yy_[i]; }
  inline const double& virialYY(size_t i) const { return virial_yy_[i]; }
  inline double& virialYZ(size_t i) { return virial_yz_[i]; }
  inline const double& virialYZ(size_t i) const { return virial_yz_[i]; }

  inline double& virialZX(size_t i) { return virial_zx_[i]; }
  inline const double& virialZX(size_t i) const { return virial_zx_[i]; }
  inline double& virialZY(size_t i) { return virial_zy_[i]; }
  inline const double& virialZY(size_t i) const { return virial_zy_[i]; }
  inline double& virialZZ(size_t i) { return virial_zz_[i]; }
  inline const double& virialZZ(size_t i) const { return virial_zz_[i]; }

  // Energies
  inline double& ePot(size_t i) {return e_pot_[i]; }
  inline const double& ePot(size_t i) const {return e_pot_[i]; }

  inline double& eKin(size_t i) {return e_kin_[i]; }
  inline const double& eKin(size_t i) const {return e_kin_[i]; }

  inline double& eTerm(size_t i) {return e_term_[i]; }
  inline const double& eTerm(size_t i) const {return e_term_[i]; }

  inline double& eInt(size_t i) {return e_int_[i]; }
  inline const double& eInt(size_t i) const {return e_int_[i]; }

  inline double& eFull(size_t i) {return e_full_[i]; }
  inline const double& eFull(size_t i) const {return e_full_[i]; }

  // Pulse of the particle 
  inline double& pulse(size_t i) {return pulse_[i]; }
  inline const double& pulse(size_t i) const {return pulse_[i]; }

  // EAM Properties
  // Electron density of the particle
  inline double& electronDensity(size_t i) {return electron_density_[i]; }
  inline const double& electronDensity(size_t i) const {return electron_density_[i]; }
  // Pair potential energy of the particle
  inline double& pairPart(size_t i) {return pair_potential_[i]; }
  inline const double& pairPart(size_t i) const {return pair_potential_[i]; }  
  */
  // Запрещаем копирование
  Particles(const Particles &) = delete;
  Particles &operator=(const Particles &) = delete;
};
#endif // PARTICLE_H
