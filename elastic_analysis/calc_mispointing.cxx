#include <iostream>
#include <cmath>

int main() {
    char prefix;
    double theta_lab;
    double mispointing_y = 999.0, mispointing_x = 999.0;

    // Ask user for inputs
    std::cout << "Enter spectrometer type (h for HMS, p for SHMS): ";
    std::cin >> prefix;
    
    std::cout << "Enter theta_lab (in degrees): ";
    std::cin >> theta_lab;

    // Convert degrees to radians
    // double theta_rad = theta_lab * M_PI / 180.0;
    double theta_rad = theta_lab;

    // Compute y-mispointing
    if (mispointing_y == 999.0) {
        if (prefix == 'h') {
            if (std::abs(theta_lab) < 40) {
                mispointing_y = 0.1 * (0.52 - 0.012 * std::abs(theta_rad) + 0.002 * std::abs(theta_rad) * std::abs(theta_rad));
            } else {
                mispointing_y = 0.1 * (0.52 - 0.012 * 40.0 + 0.002 * 40.0 * 40.0);
            }
        } else if (prefix == 'p') {
            mispointing_y = 0.1 * (-0.6);
        }
    }
    std::cout << prefix << " Mispointing_y = " << mispointing_y << " cm" << std::endl;

    // Compute x-mispointing
    if (mispointing_x == 999.0) {
        if (prefix == 'h') {
            if (std::abs(theta_lab) < 50) {
                mispointing_x = 0.1 * (2.37 - 0.086 * std::abs(theta_rad) + 0.0012 * std::abs(theta_rad) * std::abs(theta_rad));
            } else {
                mispointing_x = 0.1 * (2.37 - 0.086 * 50.0 + 0.0012 * 50.0 * 50.0);
            }
        } else if (prefix == 'p') {
            mispointing_x = 0.1 * (-1.26);
        }
    }
    std::cout << prefix << " Mispointing_x = " << mispointing_x << " cm" << std::endl;

    return 0;
}
