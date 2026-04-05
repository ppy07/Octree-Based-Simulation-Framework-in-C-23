#pragma once

#include <array>
#include <cstddef>
#include <concepts>
#include <initializer_list>
#include <cmath>
#include <type_traits>

namespace oktal {

// =======================
// Numeric concept
// =======================
template <typename T>
concept Numeric = std::integral<T> || std::floating_point<T>;

// =======================
// Vec class
// =======================
template <Numeric T, std::size_t DIM>
class Vec {
public:
    // Default constructor (clang-tidy: rely on default member init)
    constexpr Vec() noexcept = default;

    // Constant constructor
    constexpr explicit Vec(const T& value) noexcept {
        _v.fill(value);
    }

    // Initializer-list constructor
    constexpr Vec(std::initializer_list<T> list) noexcept {
        std::size_t i = 0;
        for (const T& value : list) {
            if (i >= DIM) {
                break;
            }
            _v.at(i) = value;
            ++i;
        }
    }

    // Size
    [[nodiscard]] constexpr std::size_t size() const noexcept {
        return DIM;
    }

    // Raw data access
    constexpr T* data() noexcept {
        return _v.data();
    }

    [[nodiscard]] constexpr const T* data() const noexcept {
        return _v.data();
    }

    // Element access (clang-tidy: avoid runtime subscripting)
    constexpr T& operator[](std::size_t i) {
        return _v.at(i);
    }

    constexpr const T& operator[](std::size_t i) const {
        return _v.at(i);
    }

    // Iterators
    [[nodiscard]] constexpr auto begin() noexcept {
        return _v.begin();
    }

    [[nodiscard]] constexpr auto end() noexcept {
        return _v.end();
    }

    [[nodiscard]] constexpr auto begin() const noexcept {
        return _v.begin();
    }

    [[nodiscard]] constexpr auto end() const noexcept {
        return _v.end();
    }

    // Equality
    constexpr bool operator==(const Vec& other) const noexcept {
        for (std::size_t i = 0; i < DIM; ++i) {
            if ((*this)[i] != other[i]) {
                return false;
            }
        }
        return true;
    }

    constexpr bool operator!=(const Vec& other) const noexcept {
        return !(*this == other);
    }

    // Arithmetic operators
    constexpr Vec operator+(const Vec& other) const {
        Vec result{};
        for (std::size_t i = 0; i < DIM; ++i) {
            result[i] = (*this)[i] + other[i];
        }
        return result;
    }

    constexpr Vec operator-(const Vec& other) const {
        Vec result{};
        for (std::size_t i = 0; i < DIM; ++i) {
            result[i] = (*this)[i] - other[i];
        }
        return result;
    }

    constexpr Vec operator*(const T& scalar) const {
        Vec result{};
        for (std::size_t i = 0; i < DIM; ++i) {
            result[i] = (*this)[i] * scalar;
        }
        return result;
    }

    constexpr Vec operator/(const T& scalar) const {
        Vec result{};
        for (std::size_t i = 0; i < DIM; ++i) {
            result[i] = (*this)[i] / scalar;
        }
        return result;
    }

    friend constexpr Vec operator*(const T& scalar, const Vec& v) {
        return v * scalar;
    }

    // Augmented assignments
    constexpr Vec& operator+=(const Vec& rhs) {
        for (std::size_t i = 0; i < DIM; ++i) {
            (*this)[i] += rhs[i];
        }
        return *this;
    }

    constexpr Vec& operator-=(const Vec& rhs) {
        for (std::size_t i = 0; i < DIM; ++i) {
            (*this)[i] -= rhs[i];
        }
        return *this;
    }

    template <typename S>
    requires std::is_arithmetic_v<S>
    constexpr Vec& operator*=(S scalar) {
        for (std::size_t i = 0; i < DIM; ++i) {
            (*this)[i] *= static_cast<T>(scalar);
        }
        return *this;
    }

    template <typename S>
    requires std::is_arithmetic_v<S>
    constexpr Vec& operator/=(S scalar) {
        for (std::size_t i = 0; i < DIM; ++i) {
            (*this)[i] /= static_cast<T>(scalar);
        }
        return *this;
    }

    // Magnitudes
    [[nodiscard]] constexpr T sqrMagnitude() const noexcept {
        T sum{};
        for (std::size_t i = 0; i < DIM; ++i) {
            const T v = (*this)[i];
            sum += v * v;
        }
        return sum;
    }

    template <typename U = T>
    requires std::floating_point<U>
    [[nodiscard]] U magnitude() const {
        return std::sqrt(static_cast<U>(sqrMagnitude()));
    }

    // Unary minus (disabled for unsigned types)
    template <typename U = T>
    requires (!std::unsigned_integral<U>)
    constexpr Vec operator-() const {
        Vec result{};
        for (std::size_t i = 0; i < DIM; ++i) {
            result[i] = -(*this)[i];
        }
        return result;
    }

    // Type-converting constructor (FIX: no access to other._v)
    template <typename S>
    requires std::convertible_to<S, T>
    constexpr explicit Vec(const Vec<S, DIM>& other) noexcept {
        for (std::size_t i = 0; i < DIM; ++i) {
            _v.at(i) = static_cast<T>(other[i]);
        }
    }

private:
    std::array<T, DIM> _v{}; // clang-tidy wants default member init
};

// =======================
// Common aliases
// =======================
using Vec3F = Vec<float, 3>;
using Vec3D = Vec<double, 3>;

} // namespace oktal
