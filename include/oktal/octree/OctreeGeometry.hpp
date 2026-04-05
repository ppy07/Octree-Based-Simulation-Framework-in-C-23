#pragma once

#include <cstddef>

#include <oktal/geometry/Vec.hpp>
#include <oktal/geometry/Box.hpp>
#include <oktal/octree/MortonIndex.hpp>

namespace oktal {

class OctreeGeometry {
public:
    using vector_type = Vec<double, 3>;
    using box_type    = Box<double>;

    // unit cube [0,1]^3
    OctreeGeometry() = default;

    OctreeGeometry(const vector_type& origin, double sidelength) noexcept
        : origin_{origin}, sidelength_{sidelength} {}

    [[nodiscard]] vector_type origin() const noexcept {
        return origin_;
    }

    [[nodiscard]] double sidelength() const noexcept {
        return sidelength_;
    }

    [[nodiscard]] double dx(std::size_t level) const noexcept {
        return sidelength_ / static_cast<double>(std::size_t{1} << level);
    }

    [[nodiscard]] vector_type cellExtents(std::size_t level) const noexcept {
        const double h = dx(level);
        return {h, h, h};
    }

    [[nodiscard]] vector_type cellMinCorner(const MortonIndex& m) const {
        const auto coords = m.gridCoordinates();
        const double h = dx(m.level());
        return origin_ + vector_type{
            static_cast<double>(coords[0]) * h,
            static_cast<double>(coords[1]) * h,
            static_cast<double>(coords[2]) * h
        };
    }

    [[nodiscard]] vector_type cellMaxCorner(const MortonIndex& m) const {
        const double h = dx(m.level());
        return cellMinCorner(m) + vector_type{h, h, h};
    }

    [[nodiscard]] vector_type cellCenter(const MortonIndex& m) const {
        return (cellMinCorner(m) + cellMaxCorner(m)) / 2.0;
    }

    [[nodiscard]] box_type cellBoundingBox(const MortonIndex& m) const {
        return {cellMinCorner(m), cellMaxCorner(m)};
    }

private:
    vector_type origin_{0.0, 0.0, 0.0};
    double sidelength_{1.0};
};

} // namespace oktal
