#pragma once

#include "oktal/lbm/D3Q19.hpp"
#include "oktal/octree/CellGrid.hpp"

#include <cstddef>

namespace oktal::lbm {

namespace detail {

// Discrete equilibrium distribution for D3Q19 (equation (19) from the task page):
// f_i^eq(rho, u) = w_i * rho * ( 1 + (2 c_i^T u - u^T u)/(2 c_s^2) + ( (c_i^T u)^2 )/(2 c_s^4) )
inline double feq(std::size_t i, double rho, double ux, double uy, double uz)
{
  const double cs2 = D3Q19::C_S_SQ;
  const double inv2cs2 = 1.0 / (2.0 * cs2);
  const double inv2cs4 = 1.0 / (2.0 * cs2 * cs2);

  // Use .at() on std::array to satisfy clang-tidy (constant-array-index rule)
  const auto& c = D3Q19::CS.at(i);

  // Use constant indices on oktal::Vec to satisfy clang-tidy and because Vec has no .at()
  const double cu =
      static_cast<double>(c[0uz]) * ux +
      static_cast<double>(c[1uz]) * uy +
      static_cast<double>(c[2uz]) * uz;

  const double uu = ux * ux + uy * uy + uz * uz;

  const double term = 1.0 + (2.0 * cu - uu) * inv2cs2 + (cu * cu) * inv2cs4;

  // Use .at() on std::array weights
  return D3Q19::W.at(i) * rho * term;
}

} // namespace detail

// ============================================================
// InitializePdfs  (Eq. 20)
// Signature (Listing 19):
// void (CellGrid::CellView cell,
//       D3Q19LatticeView pdfs,
//       GridVectorView<const double, 1> rho,
//       GridVectorView<const double, 3> u);
// ============================================================
class InitializePdfs {
public:
  void operator()(CellGrid::CellView cell,
                  D3Q19LatticeView pdfs,
                  GridVectorView<const double, 1> rho,
                  GridVectorView<const double, 3> u) const
  {
    const std::size_t cellIdx{cell};

    const double r  = rho[cellIdx];
    const double ux = u[cellIdx, 0uz];
    const double uy = u[cellIdx, 1uz];
    const double uz = u[cellIdx, 2uz];

    for (std::size_t i = 0uz; i < D3Q19::Q; ++i) {
      pdfs[cellIdx, i] = detail::feq(i, r, ux, uy, uz);
    }
  }
};

// ============================================================
// ComputeMacroscopicQuantities  (Eq. 17,18)
// Signature (Listing 20):
// void (CellGrid::CellView cell,
//       D3Q19LatticeConstView pdfs,
//       GridVectorView<double, 1> rho,
//       GridVectorView<double, 3> u);
// ============================================================
class ComputeMacroscopicQuantities {
public:
  void operator()(CellGrid::CellView cell,
                  D3Q19LatticeConstView pdfs,
                  GridVectorView<double, 1> rho,
                  GridVectorView<double, 3> u) const
  {
    const std::size_t cellIdx{cell};

    double r  = 0.0;
    double mx = 0.0;
    double my = 0.0;
    double mz = 0.0;

    for (std::size_t i = 0uz; i < D3Q19::Q; ++i) {
      const double fi = pdfs[cellIdx, i];
      r += fi;

      const auto& c = D3Q19::CS.at(i);
      mx += static_cast<double>(c[0uz]) * fi;
      my += static_cast<double>(c[1uz]) * fi;
      mz += static_cast<double>(c[2uz]) * fi;
    }

    rho[cellIdx] = r;

    if (r != 0.0) {
      const double invR = 1.0 / r;
      u[cellIdx, 0uz] = mx * invR;
      u[cellIdx, 1uz] = my * invR;
      u[cellIdx, 2uz] = mz * invR;
    } else {
      // Defensive: avoid NaNs in degenerate states
      u[cellIdx, 0uz] = 0.0;
      u[cellIdx, 1uz] = 0.0;
      u[cellIdx, 2uz] = 0.0;
    }
  }
};

// ============================================================
// Collide  (Eq. 21)
// Signature (Listing 19):
// void (CellGrid::CellView cell,
//       D3Q19LatticeView pdfs,
//       GridVectorView<const double, 1> rho,
//       GridVectorView<const double, 3> u);
// ============================================================
class Collide {
public:
  explicit Collide(double omega) : omega_{omega} {}

  void operator()(CellGrid::CellView cell,
                  D3Q19LatticeView pdfs,
                  GridVectorView<const double, 1> rho,
                  GridVectorView<const double, 3> u) const
  {
    const std::size_t cellIdx{cell};

    const double r  = rho[cellIdx];
    const double ux = u[cellIdx, 0uz];
    const double uy = u[cellIdx, 1uz];
    const double uz = u[cellIdx, 2uz];

    for (std::size_t i = 0uz; i < D3Q19::Q; ++i) {
      const double fi  = pdfs[cellIdx, i];
      const double feq = detail::feq(i, r, ux, uy, uz);
      pdfs[cellIdx, i] = fi + omega_ * (feq - fi);
    }
  }

private:
  double omega_;
};

// ============================================================
// Stream  (Eq. 22)
// Two-array approach: pdfsDst gets the streamed result from pdfsSrc.
// Signature (Listing 21):
// void (CellGrid::CellView cell,
//       D3Q19LatticeView pdfsDst,
//       D3Q19LatticeConstView pdfsSrc);
// ============================================================
class Stream {
public:
  void operator()(CellGrid::CellView cell,
                  D3Q19LatticeView pdfsDst,
                  D3Q19LatticeConstView pdfsSrc) const
  {
    const std::size_t cellIdx{cell};

    // i = 0: stationary population stays in place
    pdfsDst[cellIdx, 0uz] = pdfsSrc[cellIdx, 0uz];

    // i > 0: move population to neighbor cell in direction i
    // Task guarantee: neighbor lists exist for all non-zero velocities and neighbors exist.
    for (std::size_t i = 1uz; i < D3Q19::Q; ++i) {
      const auto& dir = D3Q19::CS.at(i);

      // Build Vec<std::ptrdiff_t,3> without initializer-list (robust with your Vec)
      Vec<std::ptrdiff_t, 3> step{};
      step[0uz] = static_cast<std::ptrdiff_t>(dir[0uz]);
      step[1uz] = static_cast<std::ptrdiff_t>(dir[1uz]);
      step[2uz] = static_cast<std::ptrdiff_t>(dir[2uz]);

      const std::size_t nbrIdx{static_cast<std::size_t>(cell.neighbor(step))};
      pdfsDst[nbrIdx, i] = pdfsSrc[cellIdx, i];
    }
  }
};

} // namespace oktal::lbm
