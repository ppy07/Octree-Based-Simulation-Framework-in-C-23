#include "oktal/io/VtkExport.hpp"

#include <advpt/htgfile/VtkHtgFile.hpp>
#include <cstdint>
#include <vector>

namespace oktal::io::vtk {

namespace {

/**
 * Packs a bit vector into uint8_t bytes, MSB-first.
 * VTK requires that the first node's bit is the most significant bit of the first byte.
 */
std::vector<std::uint8_t> packBits(const std::vector<bool>& bits)
{
    std::vector<std::uint8_t> out;
    if (bits.empty()) {
        return out;
    }

    std::uint8_t currentByte = 0;
    int bitCount = 0;

    for (const bool b : bits) {
        currentByte <<= 1;
        if (b) {
            currentByte |= 1;
        }
        bitCount++;

        if (bitCount == 8) {
            out.push_back(currentByte);
            currentByte = 0;
            bitCount = 0;
        }
    }

    // Pad the last byte by shifting remaining bits to the left (MSB side)
    if (bitCount > 0) {
        currentByte <<= (8 - bitCount);
        out.push_back(currentByte);
    }

    return out;
}

} // namespace

void exportOctree(const CellOctree& octree, const std::filesystem::path& filepath)
{
    using advpt::htgfile::HyperTree;
    using advpt::htgfile::SnapshotHtgFile;

    HyperTree ht;

    // 1. Geometry
    const auto& geo = octree.geometry();
    ht.xCoords = {geo.origin()[0], geo.origin()[0] + geo.sidelength()};
    ht.yCoords = {geo.origin()[1], geo.origin()[1] + geo.sidelength()};
    ht.zCoords = {geo.origin()[2], geo.origin()[2] + geo.sidelength()};

    // 2. Topology
    const std::size_t levels = octree.numberOfLevels();
    for (std::size_t l = 0; l < levels; ++l) {
        ht.nodesPerDepth.push_back(static_cast<std::int64_t>(octree.numberOfNodes(l)));
    }

    std::vector<bool> descriptorBits;
    std::vector<bool> maskBits;
    std::vector<int32_t> levelData;

    // 3. Traversal: BFS level by level
    for (std::size_t l = 0; l < levels; ++l) {
        for (const auto& node : octree.nodesStream(l)) {
            const bool maskValue = node.isPhantom() && !node.isRefined();
            maskBits.push_back(maskValue);

            if (l + 1 < levels) {
                descriptorBits.push_back(node.isRefined());
            }

            levelData.push_back(static_cast<int32_t>(l));
        }
    }

    ht.descriptor = packBits(descriptorBits);
    ht.mask = packBits(maskBits);

    // 4. File I/O
    auto file = SnapshotHtgFile::create(filepath, ht);

    // 5. level attribute
    file.writeCellData("level", std::span<const int32_t>(levelData.data(), levelData.size()));
}

// ============================================================
// Task 5b: exportCellGrid
// ============================================================

CellGridExport::CellGridExport(const CellGrid& grid, const std::filesystem::path& filepath)
    : grid_(&grid)
    , file_([&]() {
        using advpt::htgfile::HyperTree;
        using advpt::htgfile::SnapshotHtgFile;

        const auto& octree = grid.octree();

        HyperTree ht;

        // Geometry
        const auto& geo = octree.geometry();
        ht.xCoords = {geo.origin()[0], geo.origin()[0] + geo.sidelength()};
        ht.yCoords = {geo.origin()[1], geo.origin()[1] + geo.sidelength()};
        ht.zCoords = {geo.origin()[2], geo.origin()[2] + geo.sidelength()};

        // Nodes per depth
        const std::size_t levels = octree.numberOfLevels();
        for (std::size_t l = 0; l < levels; ++l) {
            ht.nodesPerDepth.push_back(static_cast<std::int64_t>(octree.numberOfNodes(l)));
        }

        std::vector<bool> descriptorBits;
        std::vector<bool> maskBits;
        std::vector<int32_t> levelData;

        for (std::size_t l = 0; l < levels; ++l) {
            for (const auto& node : octree.nodesStream(l)) {
                const bool maskValue = node.isPhantom() && !node.isRefined();
                maskBits.push_back(maskValue);

                if (l + 1 < levels) {
                    descriptorBits.push_back(node.isRefined());
                }

                levelData.push_back(static_cast<int32_t>(l));
            }
        }

        ht.descriptor = packBits(descriptorBits);
        ht.mask = packBits(maskBits);

        auto file = SnapshotHtgFile::create(filepath, ht);

        // always write level array
        file.writeCellData("level",
                           std::span<const int32_t>(levelData.data(), levelData.size()));

        return file;
    }())
{
}

CellGridExport exportCellGrid(const CellGrid& grid, const std::filesystem::path& filepath)
{
    return {grid, filepath};
}

} // namespace oktal::io::vtk
