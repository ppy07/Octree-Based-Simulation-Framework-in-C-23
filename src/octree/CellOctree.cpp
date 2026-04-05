#include "oktal/octree/CellOctree.hpp"

#include <stdexcept>
#include <string_view>
#include <vector>
#include <memory>   // <-- added
#include <string>   // <-- added

namespace oktal {

namespace {

// ------------------------------------------------------------
// Split descriptor "....|....|...." into views per level
// ------------------------------------------------------------
std::vector<std::string_view> splitLevels(std::string_view desc)
{
    std::vector<std::string_view> levels;

    std::size_t pos = 0;
    while (true) {
        const std::size_t next = desc.find('|', pos);
        levels.push_back(desc.substr(pos, next - pos));

        if (next == std::string_view::npos) {
            break;
        }
        pos = next + 1;
    }
    return levels;
}

// ------------------------------------------------------------
// Validate characters: . R P X
// ------------------------------------------------------------
void validateChars(std::string_view level)
{
    for (const char c : level) {
        if (c == '.' || c == 'R' || c == 'P' || c == 'X') {
            continue;
        }
        throw std::invalid_argument("invalid descriptor character");
    }
}

// ------------------------------------------------------------
// Count refined nodes in a level (R or X)
// ------------------------------------------------------------
std::size_t countRefined(std::string_view level)
{
    std::size_t refined = 0;
    for (const char c : level) {
        if (c == 'R' || c == 'X') {
            ++refined;
        }
    }
    return refined;
}

// ------------------------------------------------------------
// Validate descriptor structure:
//  - non-empty
//  - root exists and has size 1
//  - next level size = refinedCount*8
// ------------------------------------------------------------
void validateStructure(const std::vector<std::string_view>& levels)
{
    if (levels.empty()) {
        throw std::invalid_argument("empty descriptor");
    }
    if (levels.front().size() != 1) {
        throw std::invalid_argument("invalid root");
    }

    for (std::size_t l = 0; l < levels.size(); ++l) {
        validateChars(levels[l]);

        if (l + 1 < levels.size()) {
            const std::size_t refined = countRefined(levels[l]);
            const std::size_t expectedNext = refined * 8;

            if (levels[l + 1].size() != expectedNext) {
                throw std::invalid_argument("invalid refinement structure");
            }
        }
    }
}

} // namespace

// ======================
// Node
// ======================
CellOctree::Node::Node(bool refined, bool phantom, std::size_t childIdx) noexcept
    : data_(childIdx & indexMask)
{
    setRefined(refined);
    setPhantom(phantom);
}

bool CellOctree::Node::isRefined() const noexcept
{
    return (data_ & refinedMask) != 0;
}

bool CellOctree::Node::isPhantom() const noexcept
{
    return (data_ & phantomMask) != 0;
}

void CellOctree::Node::setRefined(bool v) noexcept
{
    if (v) {
        data_ |= refinedMask;
    } else {
        data_ &= ~refinedMask;
    }
}

void CellOctree::Node::setPhantom(bool v) noexcept
{
    if (v) {
        data_ |= phantomMask;
    } else {
        data_ &= ~phantomMask;
    }
}

std::size_t CellOctree::Node::childrenStartIndex() const noexcept
{
    return data_ & indexMask;
}

void CellOctree::Node::setChildrenStartIndex(std::size_t idx) noexcept
{
    data_ = (data_ & (refinedMask | phantomMask)) | (idx & indexMask);
}

std::size_t CellOctree::Node::childIndex(std::size_t branch) const noexcept
{
    return childrenStartIndex() + branch;
}

// ======================
// Constructors
// ======================
CellOctree::CellOctree()
{
    nodes_.emplace_back(); // root
    levelStart_.push_back(0);
    levelCount_.push_back(1);
}

CellOctree::CellOctree(const OctreeGeometry& geometry)
    : geometry_(geometry)
{
    nodes_.emplace_back(); // root
    levelStart_.push_back(0);
    levelCount_.push_back(1);
}

// ======================
// Geometry
// ======================
const OctreeGeometry& CellOctree::geometry() const noexcept
{
    return geometry_;
}

// ======================
// Nodes stream
// ======================
std::size_t CellOctree::numberOfNodes() const noexcept
{
    return nodes_.size();
}

std::size_t CellOctree::numberOfLevels() const noexcept
{
    return levelCount_.size();
}

std::size_t CellOctree::numberOfNodes(std::size_t level) const
{
    if (level >= levelCount_.size()) {
        return 0;
    }
    return levelCount_[level];
}

std::span<const CellOctree::Node> CellOctree::nodesStream() const noexcept
{
    return {nodes_.data(), nodes_.size()};
}

std::span<const CellOctree::Node> CellOctree::nodesStream(std::size_t level) const
{
    if (level >= levelCount_.size()) {
        return {};
    }

    const std::span<const Node> all(nodes_);
    return all.subspan(levelStart_[level], levelCount_[level]);
}

// ======================
// Descriptor
// ======================
CellOctree CellOctree::fromDescriptor(std::string_view desc)
{
    const auto levels = splitLevels(desc);
    validateStructure(levels);

    CellOctree tree;
    tree.nodes_.clear();
    tree.levelStart_.clear();
    tree.levelCount_.clear();

    std::size_t global = 0;
    for (const auto lvl : levels) {
        tree.levelStart_.push_back(global);
        tree.levelCount_.push_back(lvl.size());

        for (const char c : lvl) {
            bool refined = false;
            bool phantom = false;

            if (c == 'R') {
                refined = true;
            }
            if (c == 'P') {
                phantom = true;
            }
            if (c == 'X') {
                refined = true;
                phantom = true;
            }

            tree.nodes_.emplace_back(refined, phantom);
            ++global;
        }
    }

    for (std::size_t l = 0; l + 1 < tree.levelCount_.size(); ++l) {
        std::size_t next = tree.levelStart_[l + 1];

        for (std::size_t i = 0; i < tree.levelCount_[l]; ++i) {
            auto& node = tree.nodes_[tree.levelStart_[l] + i];

            if (node.isRefined()) {
                node.setChildrenStartIndex(next);
                next += 8;
            }
        }
    }

    return tree;
}

// ======================
// CellView
// ======================
CellOctree::CellView::CellView(
    const CellOctree* tree,
    std::size_t streamIdx,
    MortonIndex morton)
    : tree_(tree), streamIdx_(streamIdx), morton_(morton)
{
}

MortonIndex CellOctree::CellView::mortonIndex() const noexcept
{
    return morton_;
}

bool CellOctree::CellView::isRoot() const noexcept
{
    return morton_.isRoot();
}

bool CellOctree::CellView::isRefined() const noexcept
{
    return tree_->nodes_[streamIdx_].isRefined();
}

std::size_t CellOctree::CellView::level() const noexcept
{
    return morton_.level();
}

std::size_t CellOctree::CellView::streamIndex() const noexcept
{
    return streamIdx_;
}

Vec<double, 3> CellOctree::CellView::center() const
{
    return boundingBox().center();
}

Box<double> CellOctree::CellView::boundingBox() const
{
    return tree_->geometry_.cellBoundingBox(morton_);
}

// ======================
// Cell queries
// ======================
std::optional<CellOctree::CellView> CellOctree::getCell(MortonIndex m) const
{
    if (m.level() >= numberOfLevels()) {
        return std::nullopt;
    }

    std::size_t idx = 0;
    const auto path = m.getPath();

    for (std::size_t l = 0; l < m.level(); ++l) {
        if (!nodes_[idx].isRefined()) {
            return std::nullopt;
        }
        idx = nodes_[idx].childIndex(path[l]);
    }

    if (nodes_[idx].isPhantom()) {
        return std::nullopt;
    }

    return CellView(this, idx, m);
}

bool CellOctree::cellExists(MortonIndex m) const
{
    return getCell(m).has_value();
}

std::optional<CellOctree::CellView> CellOctree::getRootCell() const
{
    return getCell(MortonIndex{});
}

std::optional<CellOctree::CellView> CellOctree::getCell(morton_bits_t bits) const
{
    return getCell(MortonIndex(bits));
}

bool CellOctree::cellExists(morton_bits_t bits) const
{
    return cellExists(MortonIndex(bits));
}

std::optional<CellOctree::CellView>
CellOctree::getCell(std::initializer_list<morton_bits_t> bits) const
{
    if (bits.size() != 1) {
        return std::nullopt;
    }
    return getCell(*bits.begin());
}

bool CellOctree::cellExists(std::initializer_list<morton_bits_t> bits) const
{
    return getCell(bits).has_value();
}

// ======================
// Utilities (Task 5b): Uniform grid factory
// ======================
std::shared_ptr<const CellOctree>
CellOctree::createUniformGrid(std::size_t level)
{
    return createUniformGrid(OctreeGeometry{}, level);
}

std::shared_ptr<const CellOctree>
CellOctree::createUniformGrid(OctreeGeometry geom, std::size_t level)
{
    // Build a descriptor:
    // - levels 0..level-1 : all nodes are refined & phantom -> 'X'
    // - level 'level'     : all nodes are leaves & non-phantom -> '.'
    //
    // Special case level==0 => just one leaf '.'.
    std::string desc;

    auto pow8 = [](std::size_t l) -> std::size_t {
        std::size_t v = 1;
        for (std::size_t i = 0; i < l; ++i) {v *= 8;
        }
        return v;
    };

    if (level == 0) {
        desc = ".";
    } else {
        for (std::size_t l = 0; l < level; ++l) {
            if (l > 0) {
                desc.push_back('|');
            }
            desc.append(pow8(l), 'X');
        }
        desc.push_back('|');
        desc.append(pow8(level), '.');
    }

    CellOctree tree = CellOctree::fromDescriptor(desc);
    tree.geometry_ = (geom);

    return std::make_shared<const CellOctree>(std::move(tree));
}

} // namespace oktal
