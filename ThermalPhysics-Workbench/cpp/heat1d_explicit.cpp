// Original, self-contained teaching implementation of conservative 1D conduction.
// Generic illustrative material, two fixed-temperature surfaces; no flight hardware.
// Build: g++ -std=c++17 -O2 -Wall -Wextra -pedantic cpp/heat1d_explicit.cpp -o heat1d
// Run:   ./heat1d 80 8000 0.0005 results/cpp_profile.csv
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    try {
        const int cells = argc > 1 ? std::stoi(argv[1]) : 80;
        const int steps = argc > 2 ? std::stoi(argv[2]) : 8000;
        const double dt = argc > 3 ? std::stod(argv[3]) : 0.0005;
        const std::string output = argc > 4 ? argv[4] : "cpp_profile.csv";
        if (cells < 3 || cells > 1000000 || steps < 1 || steps > 10000000 || !std::isfinite(dt) || dt <= 0) {
            throw std::invalid_argument("Expect cells>=3, steps>=1, finite dt>0, within size limits");
        }
        constexpr double length_m = 0.02;
        constexpr double conductivity = 45.0;
        constexpr double density = 7800.0;
        constexpr double heat_capacity = 500.0;
        constexpr double initial_c = 20.0;
        constexpr double wall_c = 100.0;

        const double dx = length_m / cells;
        const double alpha = conductivity / (density * heat_capacity);
        const double r = alpha * dt / (dx * dx);
        // The first and last cell have an exterior half-cell resistance.
        // Their outward conductance is 2k/dx; explicit positivity requires 3r<=1.
        if (r > 1.0/3.0 + 1e-12) {
            throw std::invalid_argument("Explicit Euler stability/positivity bound exceeded: alpha*dt/dx^2 <= 1/3");
        }

        std::vector<double> t(static_cast<std::size_t>(cells), initial_c);
        std::vector<double> next(t.size());
        for (int step = 0; step < steps; ++step) {
            for (int i = 0; i < cells; ++i) {
                const double net = i == 0
                    ? (t[1] - t[0]) + 2.0*(wall_c - t[0])
                    : i == cells - 1
                    ? (t[cells-2] - t[cells-1]) + 2.0*(wall_c - t[cells-1])
                    : t[i-1] - 2.0*t[i] + t[i+1];
                next[i] = t[i] + r * net;
            }
            t.swap(next);
        }

        std::ofstream csv(output);
        if (!csv) throw std::runtime_error("Cannot open output file: " + output);
        csv << "depth_m,temperature_c\n" << std::setprecision(16);
        for (int i = 0; i < cells; ++i) {
            csv << (i + .5) * dx << ',' << t[i] << '\n';
        }
        std::cout << "Completed generic C++ conduction benchmark: "
                  << cells << " cells, " << steps*dt << " s, output=" << output << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
}
