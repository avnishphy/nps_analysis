#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>

using namespace std;

// Function to extract a floating-point value from a line
void extractValue(const std::string &line, double &value) {
    size_t colonPos = line.find(":");
    if (colonPos != std::string::npos) {
        std::string numberPart = line.substr(colonPos + 1); // Extract part after ':'
        std::istringstream iss(numberPart);
        double extractedValue;
        iss >> extractedValue; // Read the first valid number
        value = extractedValue;
    }
}


int main() {
    vector<int> runs = {6828, 6829, 6830, 6831, 6832, 6833, 6834, 6835, 6836, 6837, 6838, 6839, 6840}; // Hardcoded run numbers
    // string base_path = "/w/hallc-scshelf2102/nps/nps-ana/REPORT_OUTPUT_pass1/COIN/";
    string base_path = "/lustre24/expphy/volatile/hallc/nps/singhav/REPORT_OUTPUT/COIN/";
    double total_effective_charge = 0.0;

    for (int run : runs) {
        string filename = base_path + "coin_NPS_HMS_report_" + to_string(run) + "_0_1_-1.report";
        ifstream file(filename);
        if (!file) {
            cerr << "Error opening file: " << filename << endl;
            continue;
        }

        double charge = 0.0, hms_tracking_efficiency = 0.0, hgc_cerenkov_efficiency = 0.0, hod_efficiency = 0.0;
        string line;

        while (getline(file, line)) {
            if (line.find("BCM4A Charge:") != string::npos) {
                extractValue(line, charge);
            } else if (line.find("E SING FID TRACK EFFIC") != string::npos) {
                extractValue(line, hms_tracking_efficiency);
            } else if (line.find("Overall HGC Efficiency:") != string::npos) {
                extractValue(line, hgc_cerenkov_efficiency);
            } else if (line.find("3_of_4 EFF") != string::npos) {
                extractValue(line, hod_efficiency);
            }
        }
        file.close();

        if (charge == 0.0 || hms_tracking_efficiency == 0.0 || hgc_cerenkov_efficiency == 0.0 || hod_efficiency == 0.0) {
            cerr << "Missing values in file: " << filename << endl;
            continue;
        }

        double effective_charge = charge / (1000 * hms_tracking_efficiency * hod_efficiency * hgc_cerenkov_efficiency);
        total_effective_charge += effective_charge;

        cout << "Run: " << run << endl;
        cout << "Charge: " << charge << " uC" << endl;
        cout << "HMS Tracking Efficiency: " << hms_tracking_efficiency << endl;
        cout << "HGC Cerenkov Efficiency: " << hgc_cerenkov_efficiency << endl;
        cout << "Hodoscope Efficiency: " << hod_efficiency << endl;
        cout << "Effective Charge: " << effective_charge << " mC" << endl;
        cout << "------------------------------------" << endl;
    }

    cout << "Total Effective Charge for Combined Data: " << total_effective_charge << " mC" << endl;
    return 0;
}
