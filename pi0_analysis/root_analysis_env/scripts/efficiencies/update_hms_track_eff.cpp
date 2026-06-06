// update_hms_track_eff.cpp
// Fast C++ program to update HMS tracking efficiency in summary_all_runs.csv
// 
// Reads ROOT output files and extracts/calculates HMS tracking efficiency
// Compiles and runs directly without Python overhead
//
// Compile: g++ -O3 update_hms_track_eff.cpp `root-config --cflags --libs` -o update_hms_track_eff
// Usage: ./update_hms_track_eff --root-dir /path/to/root/files --summary-csv summary_all_runs.csv

#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TParameter.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <filesystem>
#include <algorithm>
#include <iomanip>

namespace fs = std::filesystem;

struct EfficiencyResult {
    int run = -1;
    double eff = 0.0;
    double eff_err = 0.0;
    bool found = false;
};

// ============================================================================
// Find ROOT file for a given run
// ============================================================================
std::string findRootFile(const std::string& root_dir, int run_number) {
    std::vector<std::string> patterns = {
        "nps_output_run_" + std::to_string(run_number) + ".root",
        "nps_analysis_run_" + std::to_string(run_number) + ".root",
        "run_" + std::to_string(run_number) + ".root",
        "skim_run" + std::to_string(run_number) + ".root"
    };
    
    for (const auto& pattern : patterns) {
        fs::path search_path = fs::path(root_dir) / pattern;
        if (fs::exists(search_path)) {
            return search_path.string();
        }
    }
    
    // Try wildcard search
    try {
        for (const auto& entry : fs::directory_iterator(root_dir)) {
            if (entry.path().extension() == ".root") {
                std::string filename = entry.path().filename().string();
                if (filename.find(std::to_string(run_number)) != std::string::npos) {
                    return entry.path().string();
                }
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "[WARN] Error searching directory: " << e.what() << std::endl;
    }
    
    return "";
}

// ============================================================================
// Extract efficiency from ROOT file using TParameter
// ============================================================================
EfficiencyResult extractEfficiencyFromROOT(const std::string& root_file, int run_number) {
    EfficiencyResult result;
    result.run = run_number;
    
    TFile* f = TFile::Open(root_file.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "[ERROR] Could not open ROOT file: " << root_file << std::endl;
        return result;
    }
    
    std::string dir_name = "hms_track_eff/run_" + std::to_string(run_number);
    TDirectory* d = f->GetDirectory(dir_name.c_str());
    
    if (!d) {
        f->Close();
        delete f;
        return result;
    }
    
    // Get parameter objects
    std::string eff_param_name = "hms_track_eff_run_" + std::to_string(run_number) + "_eff";
    std::string unc_param_name = "hms_track_eff_run_" + std::to_string(run_number) + "_eff_unc";
    
    d->cd();
    TParameter<double>* eff_param = (TParameter<double>*)gDirectory->Get(eff_param_name.c_str());
    TParameter<double>* unc_param = (TParameter<double>*)gDirectory->Get(unc_param_name.c_str());
    
    if (eff_param && unc_param) {
        result.eff = eff_param->GetVal();
        result.eff_err = unc_param->GetVal();
        result.found = true;
        
        std::cout << "[FOUND] Run " << run_number << ": eff=" << std::fixed << std::setprecision(4) 
                  << result.eff << "±" << result.eff_err << std::endl;
    }
    
    f->Close();
    delete f;
    
    return result;
}

// ============================================================================
// Calculate efficiency from event tree
// ============================================================================
EfficiencyResult calculateEfficiencyFromTree(const std::string& root_file, int run_number, 
                                            const std::string& tree_name = "T") {
    EfficiencyResult result;
    result.run = run_number;
    
    TFile* f = TFile::Open(root_file.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "[ERROR] Could not open ROOT file: " << root_file << std::endl;
        return result;
    }
    
    TTree* tree = (TTree*)f->Get(tree_name.c_str());
    if (!tree) {
        f->Close();
        delete f;
        return result;
    }
    
    // Disable all branches first for speed
    tree->SetBranchStatus("*", 0);
    
    // Enable only the branches we need
    tree->SetBranchStatus("H.hod.goodscinhit", 1);
    tree->SetBranchStatus("H.hod.betanotrack", 1);
    tree->SetBranchStatus("H.cal.etotnorm", 1);
    tree->SetBranchStatus("H.cer.npeSum", 1);
    tree->SetBranchStatus("H.dc.ntrack", 1);
    
    // Branch addresses
    Double_t H_hod_goodscinhit = 0;
    Double_t H_hod_betanotrack = 0;
    Double_t H_cal_etotnorm = 0;
    Double_t H_cer_npeSum = 0;
    Double_t H_dc_ntrack = 0;
    
    // Set branch addresses
    tree->SetBranchAddress("H.hod.goodscinhit", &H_hod_goodscinhit);
    tree->SetBranchAddress("H.hod.betanotrack", &H_hod_betanotrack);
    tree->SetBranchAddress("H.cal.etotnorm", &H_cal_etotnorm);
    tree->SetBranchAddress("H.cer.npeSum", &H_cer_npeSum);
    tree->SetBranchAddress("H.dc.ntrack", &H_dc_ntrack);
    
    Long64_t Nshould = 0;
    Long64_t Ndid = 0;
    
    // Process events
    Long64_t n_entries = tree->GetEntries();
    for (Long64_t i = 0; i < n_entries; ++i) {
        tree->GetEntry(i);
        
        // Check "should" criteria
        if (H_hod_goodscinhit == 1 &&
            H_hod_betanotrack > 0.5 && H_hod_betanotrack < 1.5 &&
            H_cal_etotnorm > 0.6 &&
            H_cer_npeSum > 0.6) {
            ++Nshould;
            
            // Check if has track
            if (H_dc_ntrack > 0) {
                ++Ndid;
            }
        }
    }
    
    f->Close();
    delete f;
    
    if (Nshould > 0) {
        result.eff = static_cast<double>(Ndid) / static_cast<double>(Nshould);
        result.eff_err = std::sqrt(result.eff * (1.0 - result.eff) / static_cast<double>(Nshould));
        result.found = true;
        
        std::cout << "[CALC] Run " << run_number << ": Nshould=" << Nshould << ", Ndid=" << Ndid 
                  << ", eff=" << std::fixed << std::setprecision(4) << result.eff 
                  << "±" << result.eff_err << std::endl;
    }
    
    return result;
}

// ============================================================================
// CSV parsing and updating
// ============================================================================
struct CSVRow {
    std::map<std::string, std::string> data;
    std::vector<std::string> header;
    std::vector<std::string> values;
};

std::vector<std::string> parseCSVLine(const std::string& line) {
    std::vector<std::string> result;
    std::stringstream ss(line);
    std::string item;
    
    while (std::getline(ss, item, ',')) {
        // Trim leading and trailing whitespace
        size_t start = item.find_first_not_of(" \t\r\n");
        size_t end = item.find_last_not_of(" \t\r\n");
        
        if (start != std::string::npos) {
            item = item.substr(start, end - start + 1);
        } else {
            item = "";
        }
        
        result.push_back(item);
    }
    
    return result;
}

std::string escapeCSVField(const std::string& field) {
    if (field.find(',') != std::string::npos || field.find('"') != std::string::npos) {
        std::string escaped = "\"";
        for (char c : field) {
            if (c == '"') escaped += "\"";
            escaped += c;
        }
        escaped += "\"";
        return escaped;
    }
    return field;
}

int main(int argc, char* argv[]) {
    std::string root_dir;
    std::string summary_csv;
    std::string tree_name = "T";
    bool recalculate = false;
    
    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "--root-dir" && i + 1 < argc) {
            root_dir = argv[++i];
        } else if (arg == "--summary-csv" && i + 1 < argc) {
            summary_csv = argv[++i];
        } else if (arg == "--tree-name" && i + 1 < argc) {
            tree_name = argv[++i];
        } else if (arg == "--recalculate") {
            recalculate = true;
        } else if (arg == "-h" || arg == "--help") {
            std::cout << "Usage: " << argv[0] << " [options]\n"
                      << "Options:\n"
                      << "  --root-dir <path>       Directory containing ROOT files\n"
                      << "  --summary-csv <path>    Path to summary_all_runs.csv\n"
                      << "  --tree-name <name>      ROOT tree name (default: T)\n"
                      << "  --recalculate           Recalculate from tree\n"
                      << "  -h, --help              Show this help message\n";
            return 0;
        }
    }
    
    if (root_dir.empty() || summary_csv.empty()) {
        std::cerr << "[ERROR] Missing required arguments. Use -h for help.\n";
        return 1;
    }
    
    // Validate paths
    if (!fs::exists(root_dir)) {
        std::cerr << "[ERROR] ROOT directory not found: " << root_dir << std::endl;
        return 1;
    }
    
    if (!fs::exists(summary_csv)) {
        std::cerr << "[ERROR] Summary CSV not found: " << summary_csv << std::endl;
        return 1;
    }
    
    std::cout << "[INFO] Loading CSV: " << summary_csv << std::endl;
    
    // Convert to absolute path for clarity
    fs::path csv_abs_path = fs::absolute(summary_csv);
    std::cout << "[DEBUG] Absolute CSV path: " << csv_abs_path.string() << std::endl;
    std::ifstream csv_input(summary_csv);
    if (!csv_input) {
        std::cerr << "[ERROR] Could not open CSV file\n";
        return 1;
    }
    
    std::string header_line;
    std::getline(csv_input, header_line);
    
    std::vector<std::string> header = parseCSVLine(header_line);
    
    // Find column indices for new format
    int run_col = -1;
    int eff_col = -1;
    int eff_err_col = -1;
    for (size_t i = 0; i < header.size(); ++i) {
        if (header[i] == "run") run_col = i;
        else if (header[i] == "hms_track_eff") eff_col = i;
        else if (header[i] == "hms_track_eff_err") eff_err_col = i;
    }
    // If efficiency columns are missing, add them
    if (eff_col == -1) {
        header.push_back("hms_track_eff");
        eff_col = header.size() - 1;
    }
    if (eff_err_col == -1) {
        header.push_back("hms_track_eff_err");
        eff_err_col = header.size() - 1;
    }
    if (run_col == -1) {
        std::cerr << "[ERROR] 'run' column not found in CSV\n";
        return 1;
    }
    
    // Read all rows
    std::vector<std::vector<std::string>> rows;
    std::string line;
    while (std::getline(csv_input, line)) {
        rows.push_back(parseCSVLine(line));
    }
    csv_input.close();
    
    std::cout << "[INFO] Read " << rows.size() << " runs from CSV\n";
    std::cout << "[DEBUG] Header has " << header.size() << " columns\n";
    std::cout << "[DEBUG] run_col=" << run_col << ", eff_col=" << eff_col << ", eff_err_col=" << eff_err_col << "\n";
    std::cout << "[DEBUG] Header: ";
    for (const auto& h : header) {
        std::cout << h << " | ";
    }
    std::cout << "\n";
    
    // Show first row for debugging
    if (!rows.empty() && rows[0].size() > static_cast<size_t>(eff_col)) {
        std::cout << "[DEBUG] First row run=" << rows[0][run_col] << ", eff='" << rows[0][eff_col] << "', eff_err='" << rows[0][eff_err_col] << "'\n";
    }
    
    if (rows.empty()) {
        std::cerr << "[ERROR] No rows read from CSV\n";
        return 1;
    }
    
    int updated_count = 0;
    
    // Process each run
    for (auto& row : rows) {
        if (row.empty() || row.size() <= static_cast<size_t>(run_col)) {
            std::cout << "[DEBUG] Skipping row with size " << row.size() << "\n";
            continue;
        }
        // Ensure row has enough columns for new fields
        while (row.size() <= static_cast<size_t>(eff_err_col)) {
            row.push_back("");
        }
        int run = -1;
        try {
            run = std::stoi(row[run_col]);
        } catch (...) {
            std::cout << "[DEBUG] Could not parse run number in row\n";
            continue;
        }
        std::cout << "\n[Processing] Run " << run << std::endl;
        // Find ROOT file
        std::string root_file = findRootFile(root_dir, run);
        if (root_file.empty()) {
            std::cout << "[WARN] Could not find ROOT file for run " << run << std::endl;
            continue;
        }
        std::cout << "[INFO] Found: " << fs::path(root_file).filename().string() << std::endl;
        // Try to extract from ROOT file parameters
        EfficiencyResult result = extractEfficiencyFromROOT(root_file, run);
        // If not found or recalculate, try from tree
        if (!result.found || recalculate) {
            result = calculateEfficiencyFromTree(root_file, run, tree_name);
        }
        // Update CSV if found
        if (result.found) {
            row[eff_col] = std::to_string(result.eff);
            row[eff_err_col] = std::to_string(result.eff_err);
            updated_count++;
        } else {
            std::cout << "[SKIP] Could not calculate efficiency for run " << run << std::endl;
        }
    }
    
    // Write updated CSV
    std::string backup_csv = summary_csv + ".bak";
    // Create output path WITHOUT modifying summary_csv
    std::string output_csv = summary_csv;
    size_t csv_pos = output_csv.rfind(".csv");
    if (csv_pos != std::string::npos) {
        output_csv.replace(csv_pos, 4, "_updated.csv");
    }
    
    // Convert to absolute paths
    fs::path output_csv_abs = fs::absolute(output_csv);
    fs::path backup_csv_abs = fs::absolute(backup_csv);
    fs::path original_csv_abs = fs::absolute(summary_csv);
    
    std::cout << "\n[INFO] Writing updated CSV: " << output_csv_abs.string() << std::endl;
    std::cout << "[DEBUG] Total rows to write: " << rows.size() << "\n";
    
    std::ofstream csv_output(output_csv_abs.string());
    if (!csv_output) {
        std::cerr << "[ERROR] Could not open output CSV file: " << output_csv_abs.string() << "\n";
        return 1;
    }
    
    // Write header
    for (size_t i = 0; i < header.size(); ++i) {
        if (i > 0) csv_output << ",";
        csv_output << escapeCSVField(header[i]);
    }
    csv_output << "\n";
    csv_output.flush();
    
    // Write rows
    int row_count = 0;
    for (const auto& row : rows) {
        if (row.empty()) {
            std::cout << "[DEBUG] Skipping empty row\n";
            continue;
        }
        for (size_t i = 0; i < row.size(); ++i) {
            if (i > 0) csv_output << ",";
            csv_output << escapeCSVField(row[i]);
        }
        csv_output << "\n";
        row_count++;
    }
    csv_output.flush();
    csv_output.close();
    
    std::cout << "[DEBUG] Wrote " << row_count << " data rows to " << output_csv_abs.string() << "\n";
    
    // Verify file was written
    if (!fs::exists(output_csv_abs) || fs::file_size(output_csv_abs) == 0) {
        std::cerr << "[ERROR] Output CSV file is empty or doesn't exist!\n";
        std::cerr << "[ERROR] Expected file: " << output_csv_abs.string() << "\n";
        return 1;
    }
    
    std::cout << "[SAVE] Updated CSV: " << output_csv_abs.string() << " (size: " << fs::file_size(output_csv_abs) << " bytes)\n";
    std::cout << "[INFO] Updated " << updated_count << " runs\n";
    
    // Backup and replace original
    try {
        std::cout << "\n[INFO] Creating backup: " << backup_csv_abs.string() << std::endl;
        if (fs::exists(backup_csv_abs)) {
            fs::remove(backup_csv_abs);  // Remove old backup if exists
        }
        fs::rename(original_csv_abs, backup_csv_abs);
        std::cout << "[INFO] Backup created: " << backup_csv_abs.string() << std::endl;
        
        fs::rename(output_csv_abs, original_csv_abs);
        std::cout << "[INFO] Replaced original CSV with updated version\n";
        std::cout << "[INFO] Final file: " << original_csv_abs.string() << "\n";
    } catch (const fs::filesystem_error& e) {
        std::cerr << "[ERROR] File operation failed: " << e.what() << std::endl;
        std::cerr << "[ERROR] Please manually replace:\n"
                  << "  mv " << output_csv_abs.string() << " " << original_csv_abs.string() << "\n";
        return 1;
    }
    
    return 0;
}
