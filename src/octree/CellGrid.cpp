#include "oktal/octree/CellGrid.hpp"

#include <algorithm>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace oktal {

CellGridBuilder::CellGridBuilder(std::shared_ptr<const CellOctree> ot)
    : octree_(std::move(ot)) {}

CellGridBuilder& CellGridBuilder::levels(std::span<const std::size_t> l) {
    levels_.assign(l.begin(), l.end());
    return *this;
}

CellGridBuilder& CellGridBuilder::levels(std::initializer_list<std::size_t> l) {
    levels_.assign(l.begin(), l.end());
    return *this;
}

CellGridBuilder&
CellGridBuilder::neighborhood(std::span<const Vec<std::ptrdiff_t, 3>> d) {
    neighborhood_.assign(d.begin(), d.end());
    return *this;
}

CellGridBuilder&
CellGridBuilder::neighborhood(std::initializer_list<Vec<std::ptrdiff_t, 3>> d) {
    neighborhood_.assign(d.begin(), d.end());
    return *this;
}

CellGridBuilder& CellGridBuilder::periodicityMapper(const NoPeriodicity& mapper) {
    periodicity_ = mapper;
    return *this;
}

CellGridBuilder& CellGridBuilder::periodicityMapper(const Torus& mapper) {
    periodicity_ = mapper;
    return *this;
}

CellGrid CellGridBuilder::build() const {
    CellGrid grid{octree_};

    std::vector<std::size_t> lvls = levels_;
    if (lvls.empty()) {
        for (std::size_t l = 0; l < octree_->numberOfLevels(); ++l) {
            lvls.push_back(l);
        }
    }

    for (auto l : lvls) {
        for (auto cell : octree_->horizontalRange(l)) {
            grid.morton_indices_.push_back(cell.mortonIndex());
        }
    }

    grid.stream_to_enum_.assign(
        octree_->numberOfNodes(), CellGrid::NOT_ENUMERATED);

    for (std::size_t i = 0; i < grid.morton_indices_.size(); ++i) {
        const auto c = octree_->getCell(grid.morton_indices_[i]);
        if (!c) {
            continue;
        }
        grid.stream_to_enum_[c->streamIndex()] = i;
    }

    if (neighborhood_.empty()) {
        return grid;
    }

    grid.directions_ = neighborhood_;
    grid.adjacency_.resize(neighborhood_.size(),
        std::vector<std::size_t>(grid.size(), CellGrid::NO_NEIGHBOR));

    for (std::size_t d = 0; d < neighborhood_.size(); ++d) {
        const auto& off = neighborhood_[d];

        for (std::size_t i = 0; i < grid.size(); ++i) {
            const auto& m = grid.morton_indices_[i];
            const auto level = m.level();
            const auto base = m.gridCoordinates();

            const Vec<std::ptrdiff_t, 3> target{
                static_cast<std::ptrdiff_t>(base[0]) + off[0],
                static_cast<std::ptrdiff_t>(base[1]) + off[1],
                static_cast<std::ptrdiff_t>(base[2]) + off[2]
            };

            const std::optional<Vec<std::size_t, 3>> mapped =
                std::visit([&](const auto& p) {
                    return p(target, level);
                }, periodicity_);

            if (!mapped) {
                continue;
            }

            const MortonIndex nm =
                MortonIndex::fromGridCoordinates(level, *mapped);

            const auto c = octree_->getCell(nm);

            // Ensure we found a valid cell AND it's at the same level.
            // If getCell returns an ancestor, it's not a neighbor in this grid.
            if (!c || c->mortonIndex().level() != level) {
                continue;
            }

            grid.adjacency_[d][i] =
                grid.getEnumerationIndex(c->streamIndex());
        }
    }

    return grid;
}

CellGridBuilder CellGrid::create(std::shared_ptr<const CellOctree> ot) {
    return CellGridBuilder{std::move(ot)};
}

std::span<const std::size_t>
CellGrid::neighborIndices(Vec<std::ptrdiff_t, 3> dir) const {
    for (std::size_t i = 0; i < directions_.size(); ++i) {
        if (directions_[i] == dir) {
            return adjacency_[i];
        }
    }
    throw std::out_of_range("CellGrid::neighborIndices");
}

} // namespace oktal
