#include <iostream>
#include <cmath>

int main() {
    const double mass_p = 0.938272;  // Proton mass in GeV
    double inci_energy, scat_energy, w_peak, scat_theta;

    std::cout << "Input the incident beam energy in GeV: ";
    std::cin >> inci_energy;

    std::cout << "Input the scattered beam energy in GeV: ";
    std::cin >> scat_energy;

    std::cout << "Input the W peak in GeV: ";
    std::cin >> w_peak;

    // Compute scattered angle (theta)
    double argument = (1 / (2 * inci_energy * scat_energy)) *
                      (mass_p * (inci_energy - scat_energy) - (std::pow(w_peak, 2) - std::pow(mass_p, 2)) / 2);

    if (argument < -1.0 || argument > 1.0) {
        std::cout << "Error: Invalid input values resulting in out-of-range asin argument." << std::endl;
        return 1;
    }

    scat_theta = 2 * std::asin(std::sqrt(argument));

    // Convert radians to degrees
    scat_theta = scat_theta * (180.0 / M_PI);

    std::cout << "Scattered angle should be: " << scat_theta << " degrees" << std::endl;

    return 0;
}