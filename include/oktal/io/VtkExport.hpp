#pragma once

#include <filesystem>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include <oktal/octree/CellOctree.hpp>
#include <oktal/octree/CellGrid.hpp>
#include <oktal/data/GridVector.hpp> // GridVectorView

// We need SnapshotHtgFile in the header because writeGridVector is templated.
#include <advpt/htgfile/VtkHtgFile.hpp>

namespace oktal::io::vtk {

void exportOctree(
    const CellOctree& octree,
    const std::filesystem::path& filepath
);

// ============================================================
// Export CellGrid + grid vectors (chainable)
// ============================================================

class CellGridExport {
public:
    CellGridExport(const CellGrid& grid, const std::filesystem::path& filepath);

    // ------------------------------------------------------------
    // Scalar overload: span (Task 5b)
    // ------------------------------------------------------------
    template <typename T>
    CellGridExport& writeGridVector(std::string_view name, std::span<const T> values)
    {
        if (values.size() != grid_->size()) {
            throw std::invalid_argument("writeGridVector: size mismatch");
        }

        const auto& ot = grid_->octree();

        // VTK expects one value PER STREAM NODE.
        std::vector<T> cellData(ot.numberOfNodes(), T{0});

        // Enumerated cells -> fill at streamIndex
        for (std::size_t i = 0; i < grid_->size(); ++i) {
            const MortonIndex m = grid_->mortonIndices()[i];
            const auto cv = ot.getCell(m);
            if (!cv) {
                continue;
            }
            const std::size_t s = cv->streamIndex();
            if (s < cellData.size()) {
                cellData[s] = values[i];
            }
        }

        file_.writeCellData(std::string(name),
                            std::span<const T>(cellData.data(), cellData.size()));
        return *this;
    }

    // Convenience overload (tests often pass std::vector directly)
    template <typename T>
    CellGridExport& writeGridVector(std::string_view name, const std::vector<T>& values)
    {
        return writeGridVector<T>(name, std::span<const T>(values.data(), values.size()));
    }

    // ------------------------------------------------------------
    // Task 6b: GridVectorView overload (scalar + vector)
    // Signature required by sheet/tests:
    // template <typename T, size_t Q>
    // auto &writeGridVector(const std::string &id, GridVectorView<const T, Q> vector)
    // ------------------------------------------------------------
    template <typename T, std::size_t Q>
    CellGridExport& writeGridVector(const std::string& id, ::oktal::GridVectorView<const T, Q> vec)
    {
        const auto& ot = grid_->octree();
        const std::size_t numNodes = ot.numberOfNodes();

        using E = std::remove_const_t<T>;

        // -------------------------
        // Q == 1: scalar
        // -------------------------
        if constexpr (Q == 1uz) {
            if (vec.extent(0) != grid_->size()) {
                throw std::invalid_argument("writeGridVector(GridVectorView): size mismatch");
            }

            std::vector<E> cellData(numNodes, E{0});

            for (std::size_t i = 0; i < grid_->size(); ++i) {
                const MortonIndex m = grid_->mortonIndices()[i];
                const auto cv = ot.getCell(m);
                if (!cv) {
                    continue;
                }
                const std::size_t s = cv->streamIndex();
                if (s < cellData.size()) {
                    cellData[s] = vec[i];
                }
            }

            file_.writeCellData(id, std::span<const E>(cellData.data(), cellData.size()));
            return *this;
        }

        // -------------------------
        // Q > 1: vector field
        // Convert SoA -> AoS because VTK expects layout_right (AoS)
        // -------------------------
        else {
            if (vec.extent(0) != grid_->size()) {
                throw std::invalid_argument("writeGridVector(GridVectorView): size mismatch");
            }

            // AoS storage: [node][component] contiguous
            std::vector<E> cellData(numNodes * Q, E{0});

            for (std::size_t i = 0; i < grid_->size(); ++i) {
                const MortonIndex m = grid_->mortonIndices()[i];
                const auto cv = ot.getCell(m);
                if (!cv) {
                    continue;
                }
                const std::size_t s = cv->streamIndex();
                if (s >= numNodes) {
                    continue;
                }

                for (std::size_t q = 0; q < Q; ++q) {
                    cellData[s * Q + q] = vec[i, q];
                }
            }

            // IMPORTANT: MdCellDataView must be NON-CONST element type
            // otherwise HighFive rejects it (const double not supported).
            const advpt::htgfile::MdCellDataView<E, 2> out(cellData.data(), numNodes, Q);
            file_.writeCellData(id, out);

            return *this;
        }
    }

private:
    const CellGrid* grid_{nullptr};
    advpt::htgfile::SnapshotHtgFile file_;
};

[[nodiscard]] CellGridExport exportCellGrid(
    const CellGrid& grid,
    const std::filesystem::path& filepath
);

} // namespace oktal::io::vtk
