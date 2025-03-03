#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>

struct ReportEntry {
    double physics_trigger_rate = 0.0;
    double h1x_rate = 0.0;
    double hms_trig1_live_time = 0.0;
    double hms_trig2_live_time = 0.0;
    double hms_trig3_live_time = 0.0;
    double hms_trig4_live_time = 0.0;
    double hms_trig5_live_time = 0.0;
    double hms_trig6_live_time = 0.0;
    double electronic_live_time = 0.0;
};

// Function to extract the last valid numeric value from a line
std::string extractNumber(const std::string& line) {
    std::istringstream iss(line);
    std::string word, lastValidNumber;

    while (iss >> word) {
        // Remove unwanted characters like `kHz`, `%`, `,` etc.
        word.erase(std::remove_if(word.begin(), word.end(),
                  [](unsigned char c) { return !std::isdigit(c) && c != '.' && c != '-'; }),
                  word.end());

        // If it's a valid number, store it
        if (!word.empty() && (std::isdigit(word[0]) || word[0] == '-')) {
            lastValidNumber = word;
        }
    }
    return lastValidNumber;
}

// Function to parse the report file and extract values
std::vector<ReportEntry> parseReportFiles(const std::string& filename) {
    std::vector<ReportEntry> reports;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return reports;
    }

    ReportEntry entry;
    std::string line;

    while (std::getline(file, line)) {
        if (line.find("Physics Trigger Rate") != std::string::npos) {
            std::string numStr = extractNumber(line);
            if (!numStr.empty()) entry.physics_trigger_rate = std::stod(numStr);
        } 
        else if (line.find("H1X :") != std::string::npos) {
            std::string numStr = extractNumber(line);
            if (!numStr.empty()) entry.h1x_rate = std::stod(numStr);
        } 
        else if (line.find("HMS TRIG1 Computer Live Time") != std::string::npos) {
            std::string numStr = extractNumber(line);
            if (!numStr.empty()) entry.hms_trig1_live_time = std::stod(numStr);
        } 
        else if (line.find("HMS TRIG2 Computer Live Time") != std::string::npos) {
            std::string numStr = extractNumber(line);
            if (!numStr.empty()) entry.hms_trig2_live_time = std::stod(numStr);
        } 
        else if (line.find("HMS TRIG3 Computer Live Time") != std::string::npos) {
            std::string numStr = extractNumber(line);
            if (!numStr.empty()) entry.hms_trig3_live_time = std::stod(numStr);
        } 
        else if (line.find("HMS TRIG4 Computer Live Time") != std::string::npos) {
            std::string numStr = extractNumber(line);
            if (!numStr.empty()) entry.hms_trig4_live_time = std::stod(numStr);
        } 
        else if (line.find("HMS TRIG5 Computer Live Time") != std::string::npos) {
            std::string numStr = extractNumber(line);
            if (!numStr.empty()) entry.hms_trig5_live_time = std::stod(numStr);
        } 
        else if (line.find("HMS TRIG6 Computer Live Time") != std::string::npos) {
            std::string numStr = extractNumber(line);
            if (!numStr.empty()) entry.hms_trig6_live_time = std::stod(numStr);
        } 
        else if (line.find("OG 6 GeV Electronic Live Time") != std::string::npos) {
            std::string numStr = extractNumber(line);
            if (!numStr.empty()) entry.electronic_live_time = std::stod(numStr);
        }
    }

    file.close();
    reports.push_back(entry);
    return reports;
}

// Function to print extracted data for debugging
void printReportEntries(const std::vector<ReportEntry>& reports) {
    for (const auto& entry : reports) {
        std::cout << "Physics Trigger Rate: " << entry.physics_trigger_rate << " kHz\n";
        std::cout << "H1X Rate: " << entry.h1x_rate << " kHz\n";
        std::cout << "HMS TRIG1 Live Time: " << entry.hms_trig1_live_time << " %\n";
        std::cout << "HMS TRIG2 Live Time: " << entry.hms_trig2_live_time << " %\n";
        std::cout << "HMS TRIG3 Live Time: " << entry.hms_trig3_live_time << " %\n";
        std::cout << "HMS TRIG4 Live Time: " << entry.hms_trig4_live_time << " %\n";
        std::cout << "HMS TRIG5 Live Time: " << entry.hms_trig5_live_time << " %\n";
        std::cout << "HMS TRIG6 Live Time: " << entry.hms_trig6_live_time << " %\n";
        std::cout << "Electronic Live Time: " << entry.electronic_live_time << " %\n";
        std::cout << "--------------------------------------\n";
    }
}

int main() {
    std::string filename = "your_report_file.txt"; // Change this to your actual report file path
    std::vector<ReportEntry> reports = parseReportFiles(filename);
    printReportEntries(reports);
    return 0;
}
