#include "oktal/geometry/PeriodicBox.hpp"

#include <algorithm>

namespace oktal {


std::array<double, 3>
PeriodicBox::mapIntoBox(std::array<double, 3> point) const {

  auto mapToInterval = [](double t, double lower, double upper) {
    const double intervalSize{upper - lower};
    const double tNormalized{(t - lower) / intervalSize};
    const double tInUnitInterval{tNormalized - std::floor(tNormalized)};
    return lower + tInUnitInterval * intervalSize;
  };

  for (std::size_t i = 0; i < point.size(); ++i) {
    if (periodicity_.at(i)) {
        point.at(i) = mapToInterval(point.at(i),
                                 minCorner_.at(i),
                                 maxCorner_.at(i));
    }
  }
  

  return point;
}

double PeriodicBox::sqrDistance(std::array<double, 3> pointA,
                                std::array<double, 3> pointB) const {
  auto periodicDistance = [&](std::size_t coord) {
    if (!periodicity_.at(coord)) {
      return pointA.at(coord) - pointB.at(coord);
    }

    const double lower{minCorner_.at(coord)};
    const double upper{maxCorner_.at(coord)};
    const double p0{pointA.at(coord)};
    const double p1{pointB.at(coord)};
    

    const double d0{std::abs(p0 - p1)};
    const double d1{std::abs(p0 - (lower - (upper - p1)))};
    const double d2{std::abs(p0 - (upper + (p1 - lower)))};

    return std::min({d0, d1, d2});
  };

  const double dx{periodicDistance(0)};
  const double dy{periodicDistance(1)};
  const double dz{periodicDistance(2)};

  return dx * dx + dy * dy + dz * dz;
}



} // namespace oktal
