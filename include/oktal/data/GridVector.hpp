#pragma once

#include <concepts>
#include <cstddef>
#include <experimental/mdspan>
#include <ranges>
#include <type_traits>
#include <utility>

#include <oktal/octree/CellGrid.hpp>

namespace oktal {

// ============================================================
// GridVectorView alias
// ============================================================

template <typename T, std::size_t Q>
requires std::semiregular<std::remove_const_t<T>> && (Q > 0uz)
using GridVectorView =
    std::conditional_t<
        (Q == 1uz),
        std::experimental::mdspan<
            T,
            std::experimental::extents<std::size_t, std::experimental::dynamic_extent>,
            std::experimental::layout_left>,
        std::experimental::mdspan<
            T,
            std::experimental::extents<std::size_t, std::experimental::dynamic_extent, Q>,
            std::experimental::layout_left>>;

// ============================================================
// GridVector class
// ============================================================

template <std::semiregular T, std::size_t Q>
requires (Q > 0uz)
class GridVector {
public:
    using value_type = T;
    static constexpr std::size_t NUM_COMPONENTS = Q;

    // ------------------------------------------------------------
    // Constructor / Destructor
    // ------------------------------------------------------------

    explicit GridVector(const CellGrid& cells)
        : cells_{&cells}
        , allocSize_{cells.size() * Q}
        , data_{(allocSize_ != 0uz) ? new T[allocSize_] : nullptr}
    {
    }

    ~GridVector()
    {
        delete[] data_;
    }

    // ------------------------------------------------------------
    // Copy semantics
    // ------------------------------------------------------------

    GridVector(const GridVector& other)
        : cells_{other.cells_}
        , allocSize_{other.allocSize_}
        , data_{(allocSize_ != 0uz) ? new T[allocSize_] : nullptr}
    {
        if (allocSize_ != 0uz) {
            // ranges::copy wants an OUTPUT ITERATOR, not an output range
            std::ranges::copy(std::span<const T>(other.data_, allocSize_), data_);
        }
    }

    GridVector& operator=(const GridVector& other)
    {
        if (this == &other) {
            return *this;
        }

        cells_ = other.cells_;

        if (allocSize_ != other.allocSize_) {
            delete[] data_;
            allocSize_ = other.allocSize_;
            data_ = (allocSize_ != 0uz) ? new T[allocSize_] : nullptr;
        }

        if (allocSize_ != 0uz) {
            std::ranges::copy(std::span<const T>(other.data_, allocSize_), data_);
        }

        return *this;
    }

    // ------------------------------------------------------------
    // Move semantics
    // ------------------------------------------------------------

    GridVector(GridVector&& other) noexcept
        : cells_{std::exchange(other.cells_, nullptr)}
        , allocSize_{std::exchange(other.allocSize_, 0uz)}
        , data_{std::exchange(other.data_, nullptr)}
    {
    }

    GridVector& operator=(GridVector&& other) noexcept
    {
        if (this == &other) {
            return *this;
        }

        delete[] data_;
        cells_ = std::exchange(other.cells_, nullptr);
        allocSize_ = std::exchange(other.allocSize_, 0uz);
        data_ = std::exchange(other.data_, nullptr);

        return *this;
    }

    // ------------------------------------------------------------
    // Observers
    // ------------------------------------------------------------

    [[nodiscard]] std::size_t allocSize() const noexcept { return allocSize_; }

    [[nodiscard]] T* data() noexcept { return data_; }
    [[nodiscard]] const T* data() const noexcept { return data_; }

    // ------------------------------------------------------------
    // Views
    // ------------------------------------------------------------

    [[nodiscard]] auto view() noexcept
    {
        const std::size_t numCells = allocSize_ / Q;

        if constexpr (Q == 1uz) {
            return GridVectorView<T, Q>{data_, numCells};
        } else {
            // 2D mdspan needs (cells, components)
            return GridVectorView<T, Q>{data_, numCells, Q};
        }
    }

    [[nodiscard]] auto view() const noexcept
    {
        const std::size_t numCells = allocSize_ / Q;

        if constexpr (Q == 1uz) {
            return GridVectorView<const T, Q>{data_, numCells};
        } else {
            return GridVectorView<const T, Q>{data_, numCells, Q};
        }
    }

    [[nodiscard]] auto const_view() const noexcept { return view(); }

    // ------------------------------------------------------------
    // Implicit conversion to views
    // ------------------------------------------------------------

    operator GridVectorView<T, Q>() noexcept { return view(); }
    operator GridVectorView<const T, Q>() const noexcept { return view(); }

    // ------------------------------------------------------------
    // operator[] (C++23 multidimensional subscript)
    // ------------------------------------------------------------

    template <std::size_t R = Q>
    requires (R == 1uz)
    T& operator[](std::size_t cellIdx) noexcept
    {
        return view()[cellIdx];
    }

    template <std::size_t R = Q>
    requires (R == 1uz)
    const T& operator[](std::size_t cellIdx) const noexcept
    {
        return view()[cellIdx];
    }

    template <std::size_t R = Q>
    requires (R > 1uz)
    T& operator[](std::size_t cellIdx, std::size_t compIdx) noexcept
    {
        return view()[cellIdx, compIdx];
    }

    template <std::size_t R = Q>
    requires (R > 1uz)
    const T& operator[](std::size_t cellIdx, std::size_t compIdx) const noexcept
    {
        return view()[cellIdx, compIdx];
    }

private:
    const CellGrid* cells_{nullptr};

    // declaration order matters for -Wreorder
    std::size_t allocSize_{0uz};
    T* data_{nullptr};
};

} // namespace oktal
