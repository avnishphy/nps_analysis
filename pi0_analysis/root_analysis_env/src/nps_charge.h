#ifndef NPS_CHARGE
#define NPS_CHARGE

#include <iostream>
#include "TTree.h"

namespace nps {

    /**
     * @brief Calculate accumulated charge from scaler variables already connected to a TTree.
     * 
     * @param nentries Number of entries in the TTree.
     * @param BCM2_scalerCharge Pointer to the current BCM2 scaler charge variable.
     * @param BCM2_scalerCurrent Pointer to the current BCM2 scaler current variable.
     * @param H_1MHz_scalerTime Pointer to the 1MHz scaler time variable.
     * @param run Optional run number for printing messages.
     * @param verbose If true, prints accumulated charge info.
     * @return Accumulated charge in mC, or -1 on error.
     */
    inline double get_accumulated_charge(
        Long64_t nentries,
        double* BCM2_scalerCharge,
        double* BCM2_scalerCurrent,
        double* H_1MHz_scalerTime,
        int run = 0,
        bool verbose = true
    ) {

        if (!BCM2_scalerCharge || !BCM2_scalerCurrent || !H_1MHz_scalerTime) {
            std::cerr << "Error: One or more scaler pointers are null." << std::endl;
            return -1;
        }

        if (nentries < 2) {
            std::cerr << "Error: Not enough scaler entries in run " << run << std::endl;
            return -1;
        }

        double accumulated_charge = 0.0;
        double prev_charge = 0.0;
        double prev_time = 0.0;
        bool first_entry = true;

        for (Long64_t i = 0; i < nentries; ++i) {
            // The values are assumed to already be updated by the tree->GetEntry(i) in main
            if (first_entry) {
                prev_charge = *BCM2_scalerCharge;
                prev_time = *H_1MHz_scalerTime;
                first_entry = false;
                continue;
            }

            double delta_charge = *BCM2_scalerCharge - prev_charge;
            double delta_time = *H_1MHz_scalerTime - prev_time;

            if (delta_charge > 0 && delta_time > 0) {
                accumulated_charge += delta_charge;
            }

            prev_charge = *BCM2_scalerCharge;
            prev_time = *H_1MHz_scalerTime;
        }

        // Convert from µC to mC
        double accumulated_charge_mC = accumulated_charge / 1000.0;

        if (verbose) {
            std::cout << "Run " << run << " accumulated charge: "
                      << accumulated_charge_mC << " mC" << std::endl;
        }

        return accumulated_charge_mC;
    }

} // namespace nps

#endif // NPS_CHARGE
