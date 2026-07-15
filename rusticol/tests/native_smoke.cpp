#include <rusticol.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "usage: native_smoke PROCESS_OUTPUT\n";
        return 2;
    }
    try {
        rusticol::Runtime runtime(argv[1]);
        const std::vector<double> momenta = {
            500.0, 0.0, 0.0, 500.0,
            500.0, 0.0, 0.0, -500.0,
            504.15762567199999, -304.10842628649999, 208.76026523528103,
            331.35611794513767,
            495.84237432800001, 304.10842628649999, -208.76026523528103,
            -331.35611794513767,
        };
        const auto total = runtime.evaluate(momenta, 1);
        const auto resolved = runtime.evaluate_resolved(momenta, 1);
        const auto resolved_total = resolved.total();
        if (total.size() != 1 || resolved_total.size() != 1 ||
            std::abs(total[0] - resolved_total[0]) > 1.0e-12 * std::abs(total[0])) {
            std::cerr << "resolved sum does not reproduce the compatibility total\n";
            return 1;
        }
        std::cout << std::setprecision(17) << total[0] << " " << resolved_total[0] << "\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << "\n";
        return 1;
    }
}
