#include <oktal/octree/MortonIndex.hpp>
#include <oktal/geometry/Vec.hpp>

#include <bit>
#include <stdexcept>

namespace oktal {

// ======================
// Constructors & basics
// ======================

MortonIndex::MortonIndex() noexcept
    : _bits(1)
{
}

MortonIndex::MortonIndex(morton_bits_t bits) noexcept
    : _bits(bits)
{
}

morton_bits_t MortonIndex::getBits() const noexcept
{
    return _bits;
}

// ======================
// Path conversion
// ======================

MortonIndex MortonIndex::fromPath(
    const std::vector<morton_bits_t>& path)
{
    if (path.size() > MAX_DEPTH) {
        throw std::invalid_argument(
            "MortonIndex::fromPath: path too long");
    }

    morton_bits_t bits = morton_bits_t(1) << (3 * path.size());

    for (std::size_t i = 0; i < path.size(); ++i) {
        if (path[i] > 7) {
            throw std::invalid_argument(
                "MortonIndex::fromPath: invalid digit");
        }

        bits |= path[i] << (3 * (path.size() - 1 - i));
    }

    return MortonIndex(bits);
}

std::vector<morton_bits_t> MortonIndex::getPath() const
{
    if (_bits == morton_bits_t(1)) {
        return {};
    }

    const std::size_t highest_bit =
        static_cast<std::size_t>(std::bit_width(_bits)) - 1;

    const std::size_t depth = highest_bit / 3;
    std::vector<morton_bits_t> path(depth);

    for (std::size_t i = 0; i < depth; ++i) {
        path[i] =
            (_bits >> (3 * (depth - 1 - i))) & morton_bits_t(0b111);
    }

    return path;
}

// ======================
// Position queries
// ======================

std::size_t MortonIndex::level() const noexcept
{
    return (static_cast<std::size_t>(std::bit_width(_bits)) - 1) / 3;
}

bool MortonIndex::isRoot() const noexcept
{
    return _bits == morton_bits_t(1);
}

std::size_t MortonIndex::siblingIndex() const noexcept
{
    if (isRoot()) {
        return 0;
    }

    const auto path = getPath();
    return static_cast<std::size_t>(path.back());
}

bool MortonIndex::isFirstSibling() const noexcept
{
    return isRoot() || siblingIndex() == 0;
}

bool MortonIndex::isLastSibling() const noexcept
{
    return isRoot() || siblingIndex() == 7;
}

// ======================
// Traversal
// ======================

MortonIndex MortonIndex::parent() const noexcept
{
    return MortonIndex(_bits >> 3);
}

MortonIndex MortonIndex::safeParent() const
{
    if (isRoot()) {
        throw std::logic_error(
            "MortonIndex::safeParent: root has no parent");
    }

    return parent();
}

MortonIndex MortonIndex::child(morton_bits_t digit) const noexcept
{
    return MortonIndex((_bits << 3) | digit);
}

MortonIndex MortonIndex::safeChild(morton_bits_t digit) const
{
    if (level() >= MAX_DEPTH) {
        throw std::logic_error(
            "MortonIndex::safeChild: maximum depth exceeded");
    }

    return child(digit);
}

// ======================
// Partial order
// ======================

bool MortonIndex::operator==(const MortonIndex& other) const noexcept
{
    return _bits == other._bits;
}

bool MortonIndex::operator!=(const MortonIndex& other) const noexcept
{
    return !(*this == other);
}

bool MortonIndex::operator>(const MortonIndex& other) const noexcept
{
    if (*this == other) {
        return false;
    }

    const std::size_t myLevel = level();
    const std::size_t otherLevel = other.level();

    if (myLevel >= otherLevel) {
        return false;
    }

    const std::size_t shift = 3 * (otherLevel - myLevel);
    return _bits == (other._bits >> shift);
}

bool MortonIndex::operator<(const MortonIndex& other) const noexcept
{
    return other > *this;
}

bool MortonIndex::operator>=(const MortonIndex& other) const noexcept
{
    return (*this == other) || (*this > other);
}

bool MortonIndex::operator<=(const MortonIndex& other) const noexcept
{
    return (*this == other) || (*this < other);
}

// ======================
// Grid coordinates
// ======================

Vec<std::size_t, 3> MortonIndex::gridCoordinates() const
{
    std::size_t x = 0;
    std::size_t y = 0;
    std::size_t z = 0;

    const std::size_t L = level();

    for (std::size_t i = 0; i < L; ++i) {
        const morton_bits_t d =
            (_bits >> (3 * (L - 1 - i))) & morton_bits_t(0b111);

        const std::size_t bx = d & 1;
        const std::size_t by = (d >> 1) & 1;
        const std::size_t bz = (d >> 2) & 1;

        x = (x << 1) | bx;
        y = (y << 1) | by;
        z = (z << 1) | bz;
    }

    return Vec<std::size_t, 3>{x, y, z};
}

MortonIndex MortonIndex::fromGridCoordinates(
    std::size_t level,
    const Vec<std::size_t, 3>& coords)
{
    if (level > MAX_DEPTH) {
        throw std::invalid_argument(
            "MortonIndex::fromGridCoordinates: level too large");
    }

    morton_bits_t bits = morton_bits_t(1) << (3 * level);

    for (std::size_t i = 0; i < level; ++i) {
        const std::size_t shift = level - 1 - i;

        const morton_bits_t bx = (coords[0] >> shift) & 1;
        const morton_bits_t by = (coords[1] >> shift) & 1;
        const morton_bits_t bz = (coords[2] >> shift) & 1;

        const morton_bits_t d = bx | (by << 1) | (bz << 2);
        bits |= d << (3 * (level - 1 - i));
    }

    return MortonIndex(bits);
}

} // namespace oktal
