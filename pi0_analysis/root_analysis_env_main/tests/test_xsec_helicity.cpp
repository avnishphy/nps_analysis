#include "../src/xsec_extract/xsec_helicity.h"
#include <cassert>
#include <iostream>

int main() {
    const auto balanced = nps_xsec::helicity_yield_asymmetry(50, 50, 50, 50);
    assert(balanced.valid && balanced.value == 0.0);
    assert(std::abs(balanced.error - 0.1) < 1e-14);
    // Counting limit: Var(A)=(1-A^2)/N, not 1/N away from A=0.
    const auto unequal = nps_xsec::helicity_yield_asymmetry(80, 20, 80, 20);
    assert(unequal.valid && std::abs(unequal.value - 0.6) < 1e-14);
    assert(std::abs(unequal.error - 0.08) < 1e-14);
    const auto weighted = nps_xsec::helicity_yield_asymmetry(3, 2, 5, 2);
    const auto scaled = nps_xsec::helicity_yield_asymmetry(30, 20, 500, 200);
    assert(weighted.valid && scaled.valid);
    assert(std::abs(weighted.error - std::sqrt(152.0 / 625.0)) < 1e-14);
    assert(std::abs(scaled.error - weighted.error) < 1e-14);
    assert(std::abs(scaled.value - weighted.value) < 1e-14);
    assert(!nps_xsec::helicity_yield_asymmetry(10, 0, 10, 0).valid);
    assert(!nps_xsec::helicity_yield_asymmetry(0, 10, 0, 10).valid);
    assert(!nps_xsec::helicity_yield_asymmetry(0, 0, 0, 0).valid);
    assert(!nps_xsec::helicity_yield_asymmetry(-3, 2, 5, 2).valid);
    assert(!nps_xsec::helicity_yield_asymmetry(3, 2, -5, 2).valid);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    assert(!nps_xsec::helicity_yield_asymmetry(nan, 2, 5, 2).valid);
    // Signed background subtraction is allowed when both samples have support.
    const auto signed_yield = nps_xsec::helicity_yield_asymmetry(-1, 3, 1, 5);
    assert(signed_yield.valid && signed_yield.value == -2.0);
    assert(std::abs(signed_yield.error - std::sqrt(3.5)) < 1e-14);
    std::cout << "PASS helicity weighted-yield diagnostic\n";
}
