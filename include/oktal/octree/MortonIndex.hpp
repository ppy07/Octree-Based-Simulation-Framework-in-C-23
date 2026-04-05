#pragma once

#include <cstdint>
#include <cstddef>
#include <vector>
#include <ranges>
#include <concepts>

#include <oktal/geometry/Vec.hpp>

namespace oktal {

// underlying bit type
using morton_bits_t = std::uint64_t;

class [[maybe_unused]] MortonIndex {
public:
    // maximum depth (3 bits per level)
    static constexpr std::size_t MAX_DEPTH =
        (sizeof(morton_bits_t) * 8) / 3;

    // ======================
    // Constructors & basics
    // ======================
    MortonIndex() noexcept;
    explicit MortonIndex(morton_bits_t bits) noexcept;

    [[nodiscard]] morton_bits_t getBits() const noexcept;

    // ======================
    // Path conversion
    // ======================
    static MortonIndex fromPath(const std::vector<morton_bits_t>& path);

    [[nodiscard]] std::vector<morton_bits_t> getPath() const;

    // generic range overload
    template <std::ranges::input_range R>
    requires std::convertible_to<std::ranges::range_value_t<R>, morton_bits_t>
    static MortonIndex fromPath(const R& range)
    {
        std::vector<morton_bits_t> tmp;
        for (const auto& v : range) {
            tmp.push_back(static_cast<morton_bits_t>(v));
        }
        return fromPath(tmp);
    }

    // ======================
    // Position queries
    // ======================
    [[nodiscard]] std::size_t level() const noexcept;
    [[nodiscard]] bool isRoot() const noexcept;
    [[nodiscard]] std::size_t siblingIndex() const noexcept;
    [[nodiscard]] bool isFirstSibling() const noexcept;
    [[nodiscard]] bool isLastSibling() const noexcept;

    // ======================
    // Traversal
    // ======================
    [[nodiscard]] MortonIndex parent() const noexcept;
    [[nodiscard]] MortonIndex safeParent() const;
    [[nodiscard]] MortonIndex child(morton_bits_t digit) const noexcept;
    [[nodiscard]] MortonIndex safeChild(morton_bits_t digit) const;

    // ======================
    // Partial order
    // ======================
    bool operator==(const MortonIndex& other) const noexcept;
    bool operator!=(const MortonIndex& other) const noexcept;
    bool operator<(const MortonIndex& other) const noexcept;
    bool operator>(const MortonIndex& other) const noexcept;
    bool operator<=(const MortonIndex& other) const noexcept;
    bool operator>=(const MortonIndex& other) const noexcept;

    // ======================
    // Grid coordinates
    // ======================
    [[nodiscard]] Vec<std::size_t, 3> gridCoordinates() const;

    static MortonIndex fromGridCoordinates(
        std::size_t level,
        const Vec<std::size_t, 3>& coords);

private:
    morton_bits_t _bits;
};

} // namespace oktal
