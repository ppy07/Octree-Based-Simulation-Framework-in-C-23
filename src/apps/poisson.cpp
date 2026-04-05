#include "oktal/octree/CellOctree.hpp"
#include "oktal/octree/CellGrid.hpp"
#include "oktal/io/VtkExport.hpp"
#include "oktal/geometry/Vec.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <span>
#include <string>
#include <vector>

namespace {

bool parseSizeT(const char* s, std::size_t& out) {
    try {
        const std::string str{s};
        std::size_t pos = 0;
        const unsigned long long v = std::stoull(str, &pos);

        if (pos != str.size()) { return false; }
        if (v > std::numeric_limits<std::size_t>::max()) { return false; }

        out = static_cast<std::size_t>(v);
        return true;
    } catch (...) {
        return false;
    }
}

bool parseDouble(const char* s, double& out) {
    try {
        const std::string str{s};
        std::size_t pos = 0;
        out = std::stod(str, &pos);
        return pos == str.size();
    } catch (...) {
        return false;
    }
}

double u_exact(double x, double y, double z) {
    const double pi = std::acos(-1.0);
    return std::sin(pi * x) * std::sin(pi * y) * std::sin(pi * z);
}

// For u = sin(pi x) sin(pi y) sin(pi z), we have -Δu = 3*pi^2*u
double f_rhs(double x, double y, double z) {
    const double pi = std::acos(-1.0);
    return 3.0 * pi * pi * u_exact(x, y, z);
}

} // namespace

int main(int argc, char** argv)
{
    try {
        using namespace oktal;

        const std::span<char*> args(argv, static_cast<std::size_t>(argc));

        // ./poisson <refinementLevel> <max-iterations> <epsilon> <output-file>
        if (args.size() != 5) {
            std::cerr
                << "Usage:\n"
                << "  " << (!args.empty() ? args[0] : "poisson")
                << " <refinementLevel> <max-iterations> <epsilon> <output-file>\n";
            return 1;
        }

        std::size_t level = 0;
        std::size_t maxIter = 0;
        double eps = 0.0;

        if (!parseSizeT(args[1], level) ||
            !parseSizeT(args[2], maxIter) ||
            !parseDouble(args[3], eps) ||
            eps <= 0.0) {
            std::cerr << "Error: invalid arguments. Expected level, max-iterations (>=0), epsilon (>0).\n";
            return 2;
        }

        const std::filesystem::path outFile = args[4];

        // Build uniform octree and finest-level grid
        auto ot = CellOctree::createUniformGrid(level);

        // 6-neighborhood stencil
        const std::initializer_list<Vec<std::ptrdiff_t,3>> stencil = {
            {-1, 0, 0}, {+1, 0, 0},
            { 0,-1, 0}, { 0,+1, 0},
            { 0, 0,-1}, { 0, 0,+1},
        };

        auto cells = CellGrid::create(ot)
            .levels({level})
            .neighborhood(stencil)
            .periodicityMapper(NoPeriodicity{})
            .build();

        const std::size_t n = cells.size();
        if (n == 0) {
            std::cerr << "Error: grid has zero cells.\n";
            return 3;
        }

        // Grid spacing on unit cube
        const double h  = 1.0 / static_cast<double>(std::size_t{1} << level);
        const double h2 = h * h;

        std::vector<double> f(n, 0.0);
        std::vector<double> u(n, 0.0);
        std::vector<double> uTmp(n, 0.0);
        std::vector<double> r(n, 0.0);

        // Neighbor arrays
        const auto left  = cells.neighborIndices({-1,0,0});
        const auto right = cells.neighborIndices({+1,0,0});
        const auto down  = cells.neighborIndices({0,-1,0});
        const auto up    = cells.neighborIndices({0,+1,0});
        const auto back  = cells.neighborIndices({0,0,-1});
        const auto front = cells.neighborIndices({0,0,+1});

        auto isBoundary = [&](std::size_t i) {
            return left[i]  == CellGrid::NO_NEIGHBOR ||
                   right[i] == CellGrid::NO_NEIGHBOR ||
                   down[i]  == CellGrid::NO_NEIGHBOR ||
                   up[i]    == CellGrid::NO_NEIGHBOR ||
                   back[i]  == CellGrid::NO_NEIGHBOR ||
                   front[i] == CellGrid::NO_NEIGHBOR;
        };

        // Initialize f everywhere and enforce boundary values u=phi=0
        for (std::size_t i = 0; i < n; ++i) {
            const auto c = cells[i].center();
            f[i] = f_rhs(c[0], c[1], c[2]);

            if (isBoundary(i)) {
                u[i] = 0.0;      // phi = 0
                uTmp[i] = 0.0;
            }
        }

        auto computeResidualL2 = [&]() {
            double sum2 = 0.0;

            for (std::size_t i = 0; i < n; ++i) {
                if (isBoundary(i)) {
                    r[i] = 0.0; // boundary equation is u=phi; keep residual 0 for stopping
                    continue;
                }

                const double uL = u[left[i]];
                const double uR = u[right[i]];
                const double uD = u[down[i]];
                const double uU = u[up[i]];
                const double uB = u[back[i]];
                const double uF = u[front[i]];

                // -Δ_h u = (6u - sumNeighbors)/h^2
                const double minusLap = (6.0 * u[i] - (uL + uR + uD + uU + uB + uF)) / h2;

                // residual r = f - (-Δ_h u)
                r[i] = f[i] - minusLap;

                sum2 += r[i] * r[i];
            }

            return std::sqrt(sum2 / static_cast<double>(n));
        };

        double resL2 = computeResidualL2();
        std::size_t iter = 0;

        // Jacobi iterations
        for (; iter < maxIter && resL2 > eps; ++iter) {
            for (std::size_t i = 0; i < n; ++i) {
                if (isBoundary(i)) {
                    uTmp[i] = 0.0; // phi
                    continue;
                }

                const double uL = u[left[i]];
                const double uR = u[right[i]];
                const double uD = u[down[i]];
                const double uU = u[up[i]];
                const double uB = u[back[i]];
                const double uF = u[front[i]];

                // Jacobi stencil for -Δu=f
                uTmp[i] = (uL + uR + uD + uU + uB + uF + h2 * f[i]) / 6.0;
            }

            u.swap(uTmp);
            resL2 = computeResidualL2();
        }

        // Print required info
        std::cout << "iterations " << iter << "\n";
        std::cout << "residual_l2 " << resL2 << "\n";

        // Export u, f, residual
        oktal::io::vtk::exportCellGrid(cells, outFile)
            .writeGridVector<double>("u", std::span<const double>(u))
            .writeGridVector<double>("f", std::span<const double>(f))
            .writeGridVector<double>("residual", std::span<const double>(r));

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    } catch (...) {
        std::cerr << "Unknown error\n";
        return 1;
    }
}
