#pragma once

#include <type_traits>
#include <concepts>

#include <oktal/geometry/Vec.hpp>

namespace oktal {

template <typename T = double>
requires std::is_arithmetic_v<T>
class Box {
public:
    using vector_type = Vec<T, 3>;

    // constructors
    Box() = default;

    Box(const vector_type& min, const vector_type& max)
        : minCorner_{min}, maxCorner_{max}
    {
    }

    // ======================
    // getters / setters
    // ======================

    [[nodiscard]] vector_type& minCorner() noexcept
    {
        return minCorner_;
    }

    [[nodiscard]] const vector_type& minCorner() const noexcept
    {
        return minCorner_;
    }

    [[nodiscard]] vector_type& maxCorner() noexcept
    {
        return maxCorner_;
    }

    [[nodiscard]] const vector_type& maxCorner() const noexcept
    {
        return maxCorner_;
    }

    // ======================
    // observers
    // ======================

    // center: floating-point only
    template <typename U = T>
    requires std::floating_point<U>
    [[nodiscard]] vector_type center() const noexcept
    {
        return (minCorner_ + maxCorner_) / T(2);
    }

    // extents: works for all arithmetic types
    [[nodiscard]] vector_type extents() const noexcept
    {
        return maxCorner_ - minCorner_;
    }

    // volume: works for all arithmetic types
    [[nodiscard]] T volume() const noexcept
    {
        const auto e = extents();
        return e[0] * e[1] * e[2];
    }

private:
    vector_type minCorner_{};
    vector_type maxCorner_{};
};

} // namespace oktal
