#ifndef VECTOR3_H
#define VECTOR3_H

template <typename T> class Vector3 {
  private:
    T components[3];

  public:
    // Default constructor (zero vector)
    Vector3() : components{T(0), T(0), T(0)} {}
    Vector3(const std::vector<T> &vec) {
        if (vec.size() != 3) {
            throw std::invalid_argument("Vector3 can only be initialized from a std::vector with exactly 3 elements.");
        }
        components[0] = vec[0];
        components[1] = vec[1];
        components[2] = vec[2];
    }
    // Constructor with given values.
    Vector3(T x, T y, T z) : components{x, y, z} {}

    // Accessors for x, y, z (non-const and const versions)
    T &x() { return components[0]; }
    T &y() { return components[1]; }
    T &z() { return components[2]; }

    const T &x() const { return components[0]; }
    const T &y() const { return components[1]; }
    const T &z() const { return components[2]; }

    // Overload operator[] for safe indexed access
    T &operator[](size_t index) {
        if (index >= 3)
            throw std::out_of_range("Index out of range in Vector3");
        return components[index];
    }

    const T &operator[](size_t index) const {
        if (index >= 3)
            throw std::out_of_range("Index out of range in Vector3");
        return components[index];
    }

    // Returns the squared magnitude of the vector: x*x + y*y + z*z.
    T lengthSquared() const {
        return components[0] * components[0] + components[1] * components[1] + components[2] * components[2];
    }

    // Returns the magnitude (length) of the vector: sqrt(x*x + y*y + z*z).
    T length() const { return std::sqrt(lengthSquared()); }

    Vector3<T> operator*(T scalar) const {
        return Vector3<T>(components[0] * scalar, components[1] * scalar, components[2] * scalar);
    }
    Vector3<T> operator+(T scalar) const {
        return Vector3<T>(components[0] + scalar, components[1] + scalar, components[2] + scalar);
    }
    Vector3<T> operator-(T scalar) const {
        return Vector3<T>(components[0] - scalar, components[1] - scalar, components[2] - scalar);
    }

    friend Vector3<T> operator*(T scalar, const Vector3<T> &vec) { return vec * scalar; }

    Vector3<T> &operator*=(T scalar) {
        components[0] *= scalar;
        components[1] *= scalar;
        components[2] *= scalar;
        return *this;
    }

    Vector3<T> operator/(T scalar) const {
        return Vector3<T>(components[0] / scalar, components[1] / scalar, components[2] / scalar);
    }

    Vector3<T> operator+(const Vector3<T> &other) const {
        return Vector3<T>(components[0] + other.components[0], components[1] + other.components[1],
                          components[2] + other.components[2]);
    }
    Vector3<T> operator+=(const Vector3<T> &other) {
        components[0] += other.components[0];
        components[1] += other.components[1];
        components[2] += other.components[2];
        return *this;
    }
    Vector3<T> operator-(const Vector3<T> &other) const {
        return Vector3<T>(components[0] - other.components[0], components[1] - other.components[1],
                          components[2] - other.components[2]);
    }
    // Returns a string representation of the vector.
    friend std::ostream &operator<<(std::ostream &os, const Vector3<T> &vec) {
        os << "(" << vec.components[0] << ", " << vec.components[1] << ", " << vec.components[2] << ")";
        return os;
    }
};

#endif