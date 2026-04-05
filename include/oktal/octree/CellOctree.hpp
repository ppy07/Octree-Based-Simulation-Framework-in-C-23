#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <concepts>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include <oktal/geometry/Box.hpp>
#include <oktal/geometry/Vec.hpp>
#include <oktal/octree/MortonIndex.hpp>
#include <oktal/octree/OctreeGeometry.hpp>

namespace oktal {

class CellOctree {
public:
    // ======================
    // Constructors
    // ======================
    CellOctree();
    explicit CellOctree(const OctreeGeometry& geometry);

    // ======================
    // Geometry
    // ======================
    [[nodiscard]] const OctreeGeometry& geometry() const noexcept;

    // ======================
    // Node
    // ======================
    class Node {
    public:
        Node() noexcept = default;
        Node(bool refined, bool phantom, std::size_t childIdx = 0) noexcept;

        [[nodiscard]] bool isRefined() const noexcept;
        [[nodiscard]] bool isPhantom() const noexcept;

        void setRefined(bool v) noexcept;
        void setPhantom(bool v) noexcept;

        [[nodiscard]] std::size_t childrenStartIndex() const noexcept;
        void setChildrenStartIndex(std::size_t idx) noexcept;

        [[nodiscard]] std::size_t childIndex(std::size_t branch) const noexcept;

    private:
        using storage_t = std::size_t;

        static constexpr storage_t refinedMask =
            storage_t{1} << (sizeof(storage_t) * 8 - 1);
        static constexpr storage_t phantomMask =
            storage_t{1} << (sizeof(storage_t) * 8 - 2);
        static constexpr storage_t indexMask =
            ~(refinedMask | phantomMask);

        storage_t data_ = 0;
    };

    // ======================
    // CellView
    // ======================
    class CellView {
    public:
        [[nodiscard]] MortonIndex mortonIndex() const noexcept;
        [[nodiscard]] bool isRoot() const noexcept;
        [[nodiscard]] bool isRefined() const noexcept;
        [[nodiscard]] std::size_t level() const noexcept;
        [[nodiscard]] std::size_t streamIndex() const noexcept;

        [[nodiscard]] Vec<double, 3> center() const;
        [[nodiscard]] Box<double> boundingBox() const;

    private:
        friend class CellOctree;

        CellView(const CellOctree* tree,
                 std::size_t streamIdx,
                 MortonIndex morton);

        const CellOctree* tree_{nullptr};
        std::size_t streamIdx_{0uz};
        MortonIndex morton_;
    };

    // ======================
    // Nodes stream
    // ======================
    [[nodiscard]] std::size_t numberOfNodes() const noexcept;
    [[nodiscard]] std::size_t numberOfLevels() const noexcept;
    [[nodiscard]] std::size_t numberOfNodes(std::size_t level) const;

    [[nodiscard]] std::span<const Node> nodesStream() const noexcept;
    [[nodiscard]] std::span<const Node> nodesStream(std::size_t level) const;

    // ======================
    // Descriptor
    // ======================
    [[nodiscard]] static CellOctree fromDescriptor(std::string_view desc);

    // ======================
    // Utilities (Task b - Uniform grid)
    // ======================
    [[nodiscard]] static std::shared_ptr<const CellOctree>
    createUniformGrid(std::size_t level);

    [[nodiscard]] static std::shared_ptr<const CellOctree>
    createUniformGrid(OctreeGeometry geom, std::size_t level);

    // ======================
    // Cell Queries
    // ======================
    [[nodiscard]] std::optional<CellView> getCell(MortonIndex m) const;
    [[nodiscard]] bool cellExists(MortonIndex m) const;

    [[nodiscard]] std::optional<CellView> getRootCell() const;

    [[nodiscard]] std::optional<CellView> getCell(morton_bits_t bits) const;
    [[nodiscard]] bool cellExists(morton_bits_t bits) const;

    [[nodiscard]] std::optional<CellView> getCell(std::initializer_list<morton_bits_t> bits) const;
    [[nodiscard]] bool cellExists(std::initializer_list<morton_bits_t> bits) const;

    // ============================================================
    // ===================== TASK 4b: Ranges =======================
    // ============================================================

    [[nodiscard]] auto preOrderDepthFirstRange() const;
    [[nodiscard]] auto horizontalRange(std::size_t level) const;

private:
    std::vector<Node> nodes_;
    std::vector<std::size_t> levelStart_;
    std::vector<std::size_t> levelCount_;
    OctreeGeometry geometry_;
};

// ============================================================
// ===================== TASK 4a: OctreeCursor =================
// ============================================================

class OctreeCursor {
public:
    OctreeCursor() noexcept : octree_(nullptr) {}

    explicit OctreeCursor(const CellOctree& ot)
        : octree_(&ot), path_{0uz} {}

    OctreeCursor(const CellOctree& ot, std::span<const std::size_t> path)
        : octree_(&ot), path_(path.begin(), path.end()) {}

    [[nodiscard]] const CellOctree* octree() const noexcept { return octree_; }

    [[nodiscard]] std::span<const std::size_t> path() const noexcept
    {
        return {path_.data(), path_.size()};
    }

    [[nodiscard]] bool empty() const noexcept { return octree_ == nullptr; }

    [[nodiscard]] bool end() const noexcept
    {
        return (octree_ != nullptr) && path_.empty();
    }

    [[nodiscard]] std::size_t currentLevel() const noexcept
    {
        return path_.size() - 1;
    }

    [[nodiscard]] std::size_t currentStreamIndex() const noexcept
    {
        return path_.back();
    }

    [[nodiscard]] bool firstSibling() const noexcept
    {
        if (currentLevel() == 0) {
            return true;
        }
        const auto parentIdx = path_[path_.size() - 2];
        const auto base = octree_->nodesStream()[parentIdx].childrenStartIndex();
        return currentStreamIndex() == base;
    }

    [[nodiscard]] bool lastSibling() const noexcept
    {
        if (currentLevel() == 0) {
            return true;
        }
        const auto parentIdx = path_[path_.size() - 2];
        const auto base = octree_->nodesStream()[parentIdx].childrenStartIndex();
        return currentStreamIndex() == base + 7;
    }

    [[nodiscard]] MortonIndex mortonIndex() const noexcept
    {
        if (end()) {
            return MortonIndex{};
        }
        if (currentLevel() == 0) {
            return MortonIndex{};
        }

        std::vector<morton_bits_t> digits;
        digits.reserve(currentLevel());

        for (std::size_t l = 1; l < path_.size(); ++l) {
            const std::size_t parent = path_[l - 1];
            const std::size_t node   = path_[l];

            const std::size_t base =
                octree_->nodesStream()[parent].childrenStartIndex();

            const auto digit = static_cast<morton_bits_t>(node - base);
            digits.push_back(digit);
        }

        return MortonIndex::fromPath(digits);
    }

    [[nodiscard]] std::optional<CellOctree::CellView> currentCell() const
    {
        if (end()) {
            return std::nullopt;
        }
        return octree_->getCell(mortonIndex());
    }

    friend bool operator==(const OctreeCursor& a, const OctreeCursor& b) noexcept
    {
        return (a.octree_ == b.octree_) && (a.path_ == b.path_);
    }

    friend bool operator!=(const OctreeCursor& a, const OctreeCursor& b) noexcept
    {
        return !(a == b);
    }

    void toEnd() noexcept
    {
        if (octree_ != nullptr) {
            path_.clear();
        }
    }

    void ascend() noexcept
    {
        if (end()) {
            return;
        }

        if (path_.size() <= 1) {
            toEnd();
            return;
        }

        path_.pop_back();
    }

    void descend() noexcept
    {
        if (end()) {
            return;
        }

        const auto idx = currentStreamIndex();
        const auto node = octree_->nodesStream()[idx];

        if (!node.isRefined()) {
            return;
        }

        path_.push_back(node.childrenStartIndex());
    }

    void descend(std::size_t childIdx)
    {
        if (childIdx >= 8) {
            throw std::out_of_range("OctreeCursor::descend: invalid childIdx");
        }

        if (end()) {
            return;
        }

        const auto idx = currentStreamIndex();
        const auto node = octree_->nodesStream()[idx];

        if (!node.isRefined()) {
            return;
        }

        path_.push_back(node.childrenStartIndex() + childIdx);
    }

    void previousSibling() noexcept
    {
        if (end() || currentLevel() == 0) {
            return;
        }

        if (!firstSibling()) {
            path_.back() -= 1;
        }
    }

    void nextSibling() noexcept
    {
        if (end() || currentLevel() == 0) {
            return;
        }

        if (!lastSibling()) {
            path_.back() += 1;
        }
    }

    void toSibling(std::size_t siblingIdx)
    {
        if (currentLevel() == 0) {
            if (siblingIdx != 0) {
                throw std::out_of_range("OctreeCursor::toSibling: root has no siblings");
            }
            return;
        }

        if (siblingIdx >= 8) {
            throw std::out_of_range("OctreeCursor::toSibling: invalid siblingIdx");
        }

        const auto parentIdx = path_[path_.size() - 2];
        const auto base = octree_->nodesStream()[parentIdx].childrenStartIndex();
        path_.back() = base + siblingIdx;
    }

private:
    const CellOctree* octree_;
    std::vector<std::size_t> path_;
};

// ============================================================
// ===================== TASK 4b: Iterator Policy Concept ======
// ============================================================

template <typename TPolicy>
concept OctreeIteratorPolicy =
    std::semiregular<TPolicy> &&
    requires(const TPolicy pol, OctreeCursor& c) {
        { pol.advance(c) } -> std::same_as<void>;
    };

// ============================================================
// ===================== TASK 4b: Iterator Template ============
// ============================================================

template <OctreeIteratorPolicy TPolicy>
class OctreeIterator {
public:
    using iterator_concept  = std::forward_iterator_tag;
    using iterator_category = std::forward_iterator_tag;
    using value_type        = CellOctree::CellView;
    using difference_type   = std::ptrdiff_t;
    using reference         = value_type;

    OctreeIterator() = default;

    OctreeIterator(OctreeCursor cursor, TPolicy policy)
        : cursor_(std::move(cursor)), policy_(std::move(policy))
    {
        skipPhantoms();
    }

    reference operator*() const
    {
        return cursor_.currentCell().value(); // NOLINT
    }

    OctreeIterator& operator++()
    {
        policy_.advance(cursor_);
        skipPhantoms();
        return *this;
    }

    OctreeIterator operator++(int)
    {
        OctreeIterator tmp = *this;
        ++(*this);
        return tmp;
    }

    friend bool operator==(const OctreeIterator& a, const OctreeIterator& b)
    {
        return a.cursor_ == b.cursor_;
    }

    friend bool operator!=(const OctreeIterator& a, const OctreeIterator& b)
    {
        return !(a == b);
    }

private:
    void skipPhantoms()
    {
        while (!cursor_.end() && !cursor_.currentCell().has_value()) {
            policy_.advance(cursor_);
        }
    }

    OctreeCursor cursor_;
    TPolicy policy_{};
};

// ============================================================
// ===================== TASK 4b: Range Template ===============
// ============================================================

template <OctreeIteratorPolicy TPolicy>
class OctreeCellsRange {
public:
    using iterator = OctreeIterator<TPolicy>;

    OctreeCellsRange(OctreeCursor start, OctreeCursor end, TPolicy policy)
        : start_(std::move(start)), end_(std::move(end)), policy_(std::move(policy))
    {
    }

    [[nodiscard]] iterator begin() const { return iterator{start_, policy_}; }
    [[nodiscard]] iterator end()   const { return iterator{end_, policy_}; }

private:
    OctreeCursor start_;
    OctreeCursor end_;
    TPolicy policy_;
};

// ============================================================
// ===================== TASK 4b: Policies =====================
// ============================================================

class PreOrderDepthFirstPolicy {
public:
    PreOrderDepthFirstPolicy() = default;

    static void advance(OctreeCursor& c)
    {
        if (c.end()) {
            return;
        }

        const auto idx = c.currentStreamIndex();
        if (c.octree()->nodesStream()[idx].isRefined()) {
            c.descend();
            return;
        }

        while (true) {
            if (!c.lastSibling()) {
                c.nextSibling();
                return;
            }
            c.ascend();
            if (c.end()) {
                return;
            }
        }
    }
};

class HorizontalPolicy {
public:
    HorizontalPolicy() = default;
    explicit HorizontalPolicy(std::size_t level) : level_(level) {}

    void advance(OctreeCursor& c) const
    {
        if (c.end()) {
            return;
        }

        if (level_ == 0) {
            c.toEnd();
            return;
        }

        do {
            PreOrderDepthFirstPolicy::advance(c);
            if (c.end()) {
                return;
            }
        } while (c.currentLevel() != level_);
    }

private:
    std::size_t level_{0uz};
};

// ============================================================
// ===================== Range factories =======================
// ============================================================

inline auto CellOctree::preOrderDepthFirstRange() const
{
    const OctreeCursor start{*this};
    const OctreeCursor end{*this, {}};
    return OctreeCellsRange<PreOrderDepthFirstPolicy>(start, end, PreOrderDepthFirstPolicy{});
}

inline auto CellOctree::horizontalRange(std::size_t level) const
{
    const OctreeCursor end{*this, {}};

    if (level >= numberOfLevels()) {
        return OctreeCellsRange<HorizontalPolicy>(end, end, HorizontalPolicy{level});
    }

    if (level == 0) {
        const OctreeCursor start{*this};
        return OctreeCellsRange<HorizontalPolicy>(start, end, HorizontalPolicy{0});
    }

    OctreeCursor start{*this};
    while (!start.end() && start.currentLevel() != level) {
        PreOrderDepthFirstPolicy::advance(start);
    }

    if (start.end()) {
        return OctreeCellsRange<HorizontalPolicy>(end, end, HorizontalPolicy{level});
    }

    return OctreeCellsRange<HorizontalPolicy>(start, end, HorizontalPolicy{level});
}

} // namespace oktal
