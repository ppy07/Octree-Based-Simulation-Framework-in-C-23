#pragma once

#include "oktal/octree/CellOctree.hpp"
#include "oktal/geometry/Vec.hpp"
#include "oktal/geometry/Box.hpp"
#include <cassert>
#include <array>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <variant>
#include <vector>
#include <iterator>
#include <ranges>

namespace oktal {

struct NoPeriodicity {
    constexpr NoPeriodicity() = default;

    std::optional<Vec<std::size_t,3>>
    operator()(const Vec<std::ptrdiff_t,3>& coords,
            std::size_t level) const
    {
        const std::size_t N = 1uz << level;

        for (std::size_t i = 0; i < 3uz; ++i) {
            const std::ptrdiff_t ci = coords[i];

            if (ci < 0) {
                return std::nullopt;
            }
            const auto uci = static_cast<std::size_t>(ci);
            if (uci >= N) {
                return std::nullopt;
            }
        }

        return Vec<std::size_t,3>{
            static_cast<std::size_t>(coords[0]),
            static_cast<std::size_t>(coords[1]),
            static_cast<std::size_t>(coords[2])
        };
    }

};

class Torus {
public:
    explicit Torus(std::array<bool,3> periodic)
        : periodic_(periodic) {}

    std::optional<Vec<std::size_t,3>>
    operator()(const Vec<std::ptrdiff_t,3>& coords,
            std::size_t level) const
    {
        const std::size_t N = 1uz << level;
        const auto Np = static_cast<std::ptrdiff_t>(N);

        Vec<std::size_t,3> out{};

        for (std::size_t i = 0; i < 3uz; ++i) {
            const std::ptrdiff_t ci = coords[i];

            if (periodic_.at(i)) {
                std::ptrdiff_t v = ci % Np;
                if (v < 0) {v += Np;}
                out[i] = static_cast<std::size_t>(v);
            } else {
                if (ci < 0 || ci >= Np) {
                    return std::nullopt;
                }
                out[i] = static_cast<std::size_t>(ci);
            }
        }

        return out;
    }


private:
    std::array<bool,3> periodic_;
};

class CellGridBuilder;

class CellGrid {
public:
    static constexpr std::size_t NOT_ENUMERATED =
        std::numeric_limits<std::size_t>::max();
    static constexpr std::size_t NO_NEIGHBOR = NOT_ENUMERATED;

    class CellView {
    public:
        CellView() noexcept = default;
        CellView(const CellGrid* g, std::size_t i) noexcept
            : grid_(g), idx_(i) {}

        [[nodiscard]] std::size_t enumerationIndex() const noexcept {
            return idx_;
        }

        [[nodiscard]] bool isValid() const noexcept {
            return (grid_ != nullptr) && (idx_ != NOT_ENUMERATED);
        }

        explicit operator bool() const noexcept { return isValid(); }
        operator std::size_t() const noexcept { return idx_; }

        [[nodiscard]] MortonIndex mortonIndex() const {
            return grid_->morton_indices_[idx_];
        }

        [[nodiscard]] std::size_t level() const {
            return mortonIndex().level();
        }
        [[nodiscard]] const CellOctree& octree() const
        {
            assert(grid_ != nullptr);
            return grid_->octree();
        }


        [[nodiscard]] Vec<double,3> center() const
        {
            const auto cell = octree().getCell(mortonIndex());
            assert(cell.has_value());
            return cell->center();
        }

        [[nodiscard]] Box<double> boundingBox() const
        {
            const auto cell = octree().getCell(mortonIndex());
            assert(cell.has_value());
            return cell->boundingBox();
        }


        [[nodiscard]] CellView neighbor(Vec<std::ptrdiff_t,3> dir) const {
            const auto& adj = grid_->neighborIndices(dir);
            const auto n = adj[idx_];
            if (n == NO_NEIGHBOR){ return {};}
            return {grid_, n};
        }

    private:
        const CellGrid* grid_{nullptr};
        std::size_t idx_{NOT_ENUMERATED};
    };

    class Iterator {
    public:
        using iterator_concept  = std::forward_iterator_tag;
        using iterator_category = std::forward_iterator_tag;
        using value_type        = CellView;
        using difference_type   = std::ptrdiff_t;

        Iterator() = default;
        Iterator(const CellGrid* g, std::size_t i) : grid_(g), idx_(i) {}

        value_type operator*() const { return {grid_, idx_}; }

        Iterator& operator++() { ++idx_; return *this; }
        Iterator operator++(int) { auto t = *this; ++(*this); return t; }

        friend bool operator==(const Iterator& a, const Iterator& b) {
            return a.grid_ == b.grid_ && a.idx_ == b.idx_;
        }
        friend bool operator!=(const Iterator& a, const Iterator& b) {
            return !(a == b);
        }

    private:
        const CellGrid* grid_{nullptr};
        std::size_t idx_{0};
    };

    static CellGridBuilder create(std::shared_ptr<const CellOctree> ot);

    [[nodiscard]] const CellOctree& octree() const noexcept {
        return *octree_;
    }

    [[nodiscard]] std::span<const MortonIndex>
    mortonIndices() const noexcept {
        return morton_indices_;
    }

    [[nodiscard]] std::size_t
    getEnumerationIndex(std::size_t streamIdx) const noexcept {
        if (streamIdx >= stream_to_enum_.size()) {return NOT_ENUMERATED;}
        return stream_to_enum_[streamIdx];
    }
    [[nodiscard]] CellView operator[](std::size_t i) const {
        return {this, i};
    }

    [[nodiscard]] std::size_t
    getEnumerationIndex(const CellOctree::CellView& cell) const noexcept {
        return getEnumerationIndex(cell.streamIndex());
    }

    [[nodiscard]] std::span<const std::size_t>
    neighborIndices(Vec<std::ptrdiff_t,3> dir) const;

    [[nodiscard]] Iterator begin() const { return {this, 0}; }
    [[nodiscard]] Iterator end()   const { return {this, morton_indices_.size()}; }
    [[nodiscard]] std::size_t size() const noexcept {
        return morton_indices_.size();
    }

private:
    friend class CellGridBuilder;
    explicit CellGrid(std::shared_ptr<const CellOctree> ot)
        : octree_(std::move(ot)) {}

    std::shared_ptr<const CellOctree> octree_;
    std::vector<MortonIndex> morton_indices_;
    std::vector<std::size_t> stream_to_enum_;

    std::vector<Vec<std::ptrdiff_t,3>> directions_;
    std::vector<std::vector<std::size_t>> adjacency_;
};

class CellGridBuilder {
public:
    explicit CellGridBuilder(std::shared_ptr<const CellOctree> ot);

    CellGridBuilder& levels(std::span<const std::size_t> levels);
    CellGridBuilder& levels(std::initializer_list<std::size_t> levels);

    CellGridBuilder& neighborhood(std::span<const Vec<std::ptrdiff_t,3>> directions);
    CellGridBuilder& neighborhood(std::initializer_list<Vec<std::ptrdiff_t,3>> directions);

    CellGridBuilder& periodicityMapper(const NoPeriodicity& mapper);
    CellGridBuilder& periodicityMapper(const Torus& mapper);

    [[nodiscard]] CellGrid build() const;

private:
    std::shared_ptr<const CellOctree> octree_;
    std::vector<std::size_t> levels_;
    std::vector<Vec<std::ptrdiff_t,3>> neighborhood_;
    std::variant<NoPeriodicity, Torus> periodicity_{NoPeriodicity{}};
};

} // namespace oktal