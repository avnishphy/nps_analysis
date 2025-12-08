#ifndef NPS_CHARGE_H
#define NPS_CHARGE_H

#include <iostream>
#include <cmath>
#include "TTree.h"

namespace nps {

/**
 * @brief Compute accumulated charge from a scaler TTree.
 *
 * The function expects the caller to have enabled/connected the tree branches
 * so the provided pointers are updated when scaler_tree->GetEntry(i) is called.
 *
 * Units:
 *  - scalerCurrent: microamps (µA)
 *  - scalerTime: seconds (s)
 *  - returned value: mC (milliCoulombs) because µA * s = µC and we divide by 1000.
 *
 * @param scaler_tree Pointer to the scaler TTree.
 * @param nentries Number of entries in the tree (pass 0 to auto-query from the tree).
 * @param scalerCurrent Pointer to the scaler current variable (e.g. &H_BCM4A_scalerCurrent).
 * @param scalerTime Pointer to the scaler time variable (e.g. &H_1MHz_scalerTime).
 * @param min_current Minimum current (µA) required to include the interval in the sum (default 2.0 µA).
 * @param run Optional run number for verbose output.
 * @param verbose If true prints summary and warnings.
 * @return accumulated charge in mC, or -1 on error / if accumulated charge <= 0.
 */
inline double get_accumulated_charge(
    TTree* scaler_tree,
    Long64_t nentries,
    double* scalerCurrent,
    double* scalerTime,
    double min_current = 2.0,
    int run = 0,
    bool verbose = true
) {
    if (!scaler_tree) {
        std::cerr << "nps::get_accumulated_charge: scaler_tree is null\n";
        return -1;
    }
    if (!scalerCurrent || !scalerTime) {
        std::cerr << "nps::get_accumulated_charge: one or more scaler pointers are null\n";
        return -1;
    }

    if (nentries <= 0) nentries = scaler_tree->GetEntries();
    if (nentries < 2) {
        if (verbose) std::cerr << "nps::get_accumulated_charge: not enough entries in run " << run << "\n";
        return -1;
    }

    double accumulated_charge_uC = 0.0; // µC
    double prev_time = 0.0;
    bool first = true;

    for (Long64_t i = 0; i < nentries; ++i) {
        scaler_tree->GetEntry(i);

        if (first) {
            prev_time = *scalerTime;
            first = false;
            continue;
        }

        double delta_time = (*scalerTime) - prev_time;

        // require positive time increment and current above floor
        if (delta_time > 0.0 && (*scalerCurrent) > min_current) {
            // current in µA * time in s -> µC
            accumulated_charge_uC += (*scalerCurrent) * delta_time;
        }

        prev_time = *scalerTime;
    }

    if (accumulated_charge_uC <= 0.0) {
        if (verbose) std::cerr << "nps::get_accumulated_charge: accumulated_charge <= 0 for run "
                               << run << "\n";
        return -1;
    }

    double accumulated_charge_mC = accumulated_charge_uC / 1000.0; // µC -> mC

    if (verbose) {
        std::cout << "nps::get_accumulated_charge: run " << run
                  << " accumulated charge = " << accumulated_charge_mC << " mC"
                  << " (min_current=" << min_current << " µA)\n";
    }

    return accumulated_charge_mC;
}

} // namespace nps

#endif // NPS_CHARGE_H
