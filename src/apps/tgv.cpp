// src/apps/tgv.cpp
#include <algorithm>   // std::ranges::for_each
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <ranges>
#include <span>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>     // std::swap

#include "oktal/data/GridVector.hpp"
#include "oktal/io/VtkExport.hpp"
#include "oktal/lbm/D3Q19.hpp"
#include "oktal/lbm/LbmKernels.hpp"
#include "oktal/lbm/TaylorGreen.hpp"
#include "oktal/octree/CellGrid.hpp"
#include "oktal/octree/CellOctree.hpp"

namespace fs = std::filesystem;

namespace {

// ---------------------------
// CLI helpers
// ---------------------------
int usage(std::string_view argv0)
{
  std::cerr << "Usage: " << argv0 << " <refinement-level> <output-directory>\n"
            << "  refinement-level must be at least 5\n";
  return 1;
}

bool parseSizeT(std::string_view s, std::size_t& out)
{
  try {
    const std::string str{s};
    std::size_t pos = 0;
    const unsigned long long v = std::stoull(str, &pos);
    if (pos != str.size()) {
      return false;
    }
    out = static_cast<std::size_t>(v);
    return true;
  } catch (...) {
    return false;
  }
}

// ---------------------------
// Geometry helper: cell center in physical coordinates
// ---------------------------
oktal::Vec3D cellCenterPhysical(oktal::CellGrid::CellView cell)
{
  // CellGrid::CellView provides center() in physical coordinates.
  return cell.center();
}

// ---------------------------
// Taylor-Green adapters
// ---------------------------
double evalAnalyticRho(const oktal::lbm::TaylorGreen& tg,
                       oktal::CellGrid::CellView cell,
                       double t)
{
  const oktal::Vec3D x = cellCenterPhysical(cell);
  return tg.rho(x, t);
}

oktal::Vec3D evalAnalyticU(const oktal::lbm::TaylorGreen& tg,
                           oktal::CellGrid::CellView cell,
                           double t)
{
  const oktal::Vec3D x = cellCenterPhysical(cell);
  return tg.u(x, t);
}

} // namespace

int main(int argc, char** argv)
{
  try {
    using namespace oktal;

    const std::span<char*> args{argv, static_cast<std::size_t>(argc)};
    const std::string_view argv0 =
        (!args.empty() && args[0] != nullptr) ? std::string_view{args[0]} : std::string_view{"tgv"};

    // ---------------------------
    // CLI parsing
    // ---------------------------
    if (args.size() != 3uz) {
      return usage(argv0);
    }

    std::size_t refinement = 0uz;
    if (!parseSizeT(std::string_view{args[1]}, refinement)) {
      return usage(argv0);
    }

    if (refinement < 5uz) {
      std::cerr << "Error: refinement-level must be at least 5.\n";
      return 2;
    }

    const fs::path outDir{std::string{std::string_view{args[2]}}};
    std::error_code ec;
    fs::create_directories(outDir, ec);
    if (ec) {
      std::cerr << "Error: failed to create output directory " << outDir << "\n";
      return 3;
    }

    // ---------------------------
    // Taylor-Green parameters
    // ---------------------------
    const lbm::TaylorGreen tg{refinement};

    const double dt = tg.dt();
    const std::size_t nSteps = tg.numberOfTimesteps();
    const double omega = tg.omega();

    // ---------------------------
    // Grid setup
    // ---------------------------
    auto octree = CellOctree::createUniformGrid(refinement);

    auto cells = CellGrid::create(octree)
                     .levels({refinement})
                     .neighborhood(lbm::D3Q19::CS_NO_CENTER)
                     .periodicityMapper(Torus({true, true, true}))
                     .build();

    // ---------------------------
    // Fields
    // ---------------------------
    lbm::D3Q19Lattice pdfs{cells};
    lbm::D3Q19Lattice pdfsTmp{cells};

    GridVector<double, 1> rho{cells};
    GridVector<double, 3> u{cells};
    GridVector<double, 3> uErr{cells};

    // ---------------------------
    // Kernels
    // ---------------------------
    lbm::ComputeMacroscopicQuantities computeMQ;
    lbm::Collide collide{omega};
    lbm::Stream stream;
    lbm::InitializePdfs init;

    // ---------------------------
    // Initialize rho,u at t=0 and initialize PDFs from equilibrium
    // ---------------------------
    {
      auto rhoV = rho.view();
      auto uV = u.view();

      for (auto cell : cells) {
        const std::size_t idx{cell};

        rhoV[idx] = evalAnalyticRho(tg, cell, 0.0);

        const Vec3D ua = evalAnalyticU(tg, cell, 0.0);
        uV[idx, 0uz] = ua[0uz];
        uV[idx, 1uz] = ua[1uz];
        uV[idx, 2uz] = ua[2uz];
      }

      std::ranges::for_each(cells, [&](auto cellView) {
        init(cellView, pdfs.view(), rho.const_view(), u.const_view());
      });
    }

    // After initializing PDFs, recompute macroscopic quantities from PDFs
    std::ranges::for_each(cells, [&](auto cellView) {
      computeMQ(cellView, pdfs.const_view(), rho.view(), u.view());
    });

    // ---------------------------
    // Export helper (time series)
    // ---------------------------
    const std::size_t exportEvery = 50uz;
    std::size_t exportCounter = 0uz;

    auto exportState = [&]() {
      const fs::path filename =
          outDir / ("step" + std::to_string(exportCounter) + ".vtkhdf");
      ++exportCounter;

      io::vtk::exportCellGrid(cells, filename.string())
          .writeGridVector<double, 1>("rho", rho.const_view())
          .writeGridVector<double, 3>("u", u.const_view())
          .writeGridVector<double, 3>("u_err", uErr.const_view());
    };

    // ---------------------------
    // Integration loop + timing + progress
    // ---------------------------
    const auto tStart = std::chrono::steady_clock::now();

    for (std::size_t step = 0uz; step < nSteps; ++step) {

      // 1) collision (in-place) using current rho,u
      std::ranges::for_each(cells, [&](auto cellView) {
        collide(cellView, pdfs.view(), rho.const_view(), u.const_view());
      });

      // 2) streaming: pdfs -> pdfsTmp
      std::ranges::for_each(cells, [&](auto cellView) {
        stream(cellView, pdfsTmp.view(), pdfs.const_view());
      });

      // 3) swap buffers
      std::swap(pdfs, pdfsTmp);

      // 4) macroscopic quantities from UPDATED pdfs
      std::ranges::for_each(cells, [&](auto cellView) {
        computeMQ(cellView, pdfs.const_view(), rho.view(), u.view());
      });

      // 5) compute error field vs analytic at physical time
      const double time = static_cast<double>(step + 1uz) * dt;

      {
        auto uV = u.view();
        auto eV = uErr.view();

        for (auto cell : cells) {
          const std::size_t idx{cell};
          const Vec3D ua = evalAnalyticU(tg, cell, time);

          eV[idx, 0uz] = std::abs(uV[idx, 0uz] - ua[0uz]);
          eV[idx, 1uz] = std::abs(uV[idx, 1uz] - ua[1uz]);
          eV[idx, 2uz] = std::abs(uV[idx, 2uz] - ua[2uz]);
        }
      }

      // Export + progress report
      if (((step % exportEvery) == 0uz) || ((step + 1uz) == nSteps)) {
        const auto tNow = std::chrono::steady_clock::now();
        const double elapsed = std::chrono::duration<double>(tNow - tStart).count();
        const double frac =
            static_cast<double>(step + 1uz) / static_cast<double>(nSteps);
        const double estTotal = (frac > 0.0) ? (elapsed / frac) : 0.0;
        const double remaining = (estTotal > elapsed) ? (estTotal - elapsed) : 0.0;

        std::cout << "[tgv] step " << (step + 1uz) << "/" << nSteps
                  << " | elapsed " << elapsed << " s"
                  << " | eta " << remaining << " s\n";

        exportState();
      }
    }

    const auto tEnd = std::chrono::steady_clock::now();
    const double runtime = std::chrono::duration<double>(tEnd - tStart).count();
    std::cout << "[tgv] runtime: " << runtime << " s\n";

    // ---------------------------
    // Final relative L2 errors Ex, Ey
    // ---------------------------
    const double tFinal = static_cast<double>(nSteps) * dt;

    double numX = 0.0;
    double denX = 0.0;
    double numY = 0.0;
    double denY = 0.0;

    const auto uC = u.const_view();

    for (auto cell : cells) {
      const std::size_t idx{cell};
      const Vec3D ua = evalAnalyticU(tg, cell, tFinal);

      const double ex = uC[idx, 0uz] - ua[0uz];
      const double ey = uC[idx, 1uz] - ua[1uz];

      numX += ex * ex;
      numY += ey * ey;

      denX += ua[0uz] * ua[0uz];
      denY += ua[1uz] * ua[1uz];
    }

    const double Ex = (denX > 0.0) ? std::sqrt(numX / denX) : 0.0;
    const double Ey = (denY > 0.0) ? std::sqrt(numY / denY) : 0.0;

    // Write errors.txt
    {
      const fs::path errFile = outDir / "errors.txt";
      std::ofstream os(errFile);
      if (!os) {
        std::cerr << "Error: cannot write " << errFile << "\n";
        return 4;
      }
      os << Ex << " " << Ey << "\n";
    }

    // Final export
    {
      const fs::path filename = outDir / "final.vtkhdf";
      io::vtk::exportCellGrid(cells, filename.string())
          .writeGridVector<double, 1>("rho", rho.const_view())
          .writeGridVector<double, 3>("u", u.const_view())
          .writeGridVector<double, 3>("u_err", uErr.const_view());
    }

    return 0;

  } catch (const std::exception& e) {
    std::cerr << "[tgv] fatal error: " << e.what() << "\n";
    return 10;
  } catch (...) {
    std::cerr << "[tgv] fatal error: unknown exception\n";
    return 11;
  }
}
