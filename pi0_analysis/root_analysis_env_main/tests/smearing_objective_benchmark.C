// Deterministic performance/regression fixture; compile with the fitter's flags.
// SMEARING_SOURCE may select a preserved pre-change source for exact comparison.
// Usage: benchmark OUTPUT_PREFIX [EVENTS=12000] [NSMEAR=5] [REPETITIONS=3]
// Compare the four OUTPUT_PREFIX_*.bin files before/after with cmp. Each contains
// histogram bins, sumw2 and objective metrics; timings exclude cache construction.
#ifndef SMEARING_SOURCE
#define SMEARING_SOURCE "../src/simulation_smearing/nps_sim_smearing_new.C"
#endif
#define main smearing_production_main_for_benchmark
#include SMEARING_SOURCE
#undef main

#include <cstdint>
#include <cstring>

template <typename T>
void append(std::vector<unsigned char>& bytes, const T& value) {
    const auto* p = reinterpret_cast<const unsigned char*>(&value);
    bytes.insert(bytes.end(), p, p + sizeof(value));
}

void appendHistogram(std::vector<unsigned char>& bytes, const FastHistogram1D& h) {
    append(bytes, h.nbins);
    append(bytes, h.xmin);
    append(bytes, h.xmax);
    for (double value : h.content) append(bytes, value);
    for (double value : h.error2) append(bytes, value);
}

void appendMetrics(std::vector<unsigned char>& bytes, const HistObjectiveMetrics& m) {
    append(bytes, m.chi2);
    append(bytes, m.informative_bins);
    append(bytes, m.empty_sim_data_positive_bins);
    append(bytes, m.data_integral);
    append(bytes, m.sim_integral);
}

void appendSnapshot(std::vector<unsigned char>& bytes, const FastObjectiveSnapshot& s) {
    appendHistogram(bytes, s.mpi0);
    appendHistogram(bytes, s.mmiss);
    appendHistogram(bytes, s.mpgg2);
    appendMetrics(bytes, s.breakdown.mpi0);
    appendMetrics(bytes, s.breakdown.mmiss);
    appendMetrics(bytes, s.breakdown.mpgg2);
    append(bytes, s.breakdown.total(Config::W_MPI0, Config::W_MMISS, Config::W_MPGG2));
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " OUTPUT_PREFIX [EVENTS=12000] [NSMEAR=5] [REPETITIONS=3]\n";
        return 2;
    }
    const size_t count = argc > 2 ? std::stoul(argv[2]) : 12000;
    const int nsmear = argc > 3 ? std::stoi(argv[3]) : 5;
    const int repetitions = argc > 4 ? std::stoi(argv[4]) : 3;
    TH1::AddDirectory(false);
    Config::BEAM_ENERGY = 10.6;
    Config::NPS_Z_NPS_CM = 407.0;
    Config::NPS_THETA_NPS_DEG = 17.51;

    TH1D data_mpi0("data_mpi0", "", Config::MGGAMMA_NBINS, Config::MGGAMMA_MIN, Config::MGGAMMA_MAX);
    TH1D data_mmiss("data_mmiss", "", Config::MMISS_NBINS, Config::MMISS_MIN, Config::MMISS_MAX);
    TH1D data_mpgg2("data_mpgg2", "", Config::MPGG2_NBINS, Config::MPGG2_MIN, Config::MPGG2_MAX);
    for (TH1D* h : {&data_mpi0, &data_mmiss, &data_mpgg2}) {
        h->Sumw2();
        for (int i = 1; i <= h->GetNbinsX(); ++i) {
            const double weight = 20.0 + ((i * 13) % 31);
            h->SetBinContent(i, weight);
            h->SetBinError(i, std::sqrt(weight * 1.07));
        }
    }
    std::vector<ClusterPair> events;
    events.reserve(count);
    TRandom3 input(87912);
    for (size_t i = 0; i < count; ++i) {
        ClusterPair e;
        e.e1 = input.Uniform(2.1, 3.0);
        e.e2 = input.Uniform(2.0, 2.9);
        e.x1 = input.Uniform(-20.0, -5.0);
        e.x2 = e.x1 + input.Uniform(16.0, 26.0);
        e.y1 = input.Uniform(-24.0, 24.0);
        e.y2 = e.y1 + input.Uniform(-4.0, 4.0);
        e.weight = i % 37 == 0 ? -0.2 : input.Uniform(0.05, 2.0);
        e.px_e = input.Uniform(1.4, 1.8);
        e.py_e = input.Uniform(-0.05, 0.05);
        e.pz_e = input.Uniform(4.7, 5.1);
        e.Ee = std::sqrt(e.px_e * e.px_e + e.py_e * e.py_e + e.pz_e * e.pz_e + 0.000511 * 0.000511);
        e.photon1_in_section = i % 3 != 0;
        e.photon2_in_section = i % 4 != 0;
        e.mu_a1_ext = 0.01;
        e.mu1_ext = 0.995;
        e.mu_c1_ext = -0.002;
        e.sigma1_ext = i % 5 == 0 ? 0.0 : 0.024;
        e.sigma_pos1_ext = i % 7 == 0 ? 0.0 : 0.45;
        e.mu_a2_ext = -0.01;
        e.mu2_ext = 1.008;
        e.mu_c2_ext = 0.003;
        e.sigma2_ext = i % 5 == 0 ? 0.0 : 0.038;
        e.sigma_pos2_ext = i % 7 == 0 ? 0.0 : 0.65;
        if (i % 997 == 0) e.e1 = 0.0;
        if (i % 991 == 0) e.e2 = -0.5;
        events.push_back(e);
    }

    std::vector<std::array<double, 6>> points;
    for (int i = 0; i < 17; ++i) {
        points.push_back({0.003 * (i % 5 - 2), 0.987 + 0.002 * i,
                          0.001 * (i % 7 - 3), i % 4 == 0 ? 0.0 : 0.02 + 0.002 * i,
                          i % 3 == 0 ? 0.0 : 0.3 + 0.04 * i,
                          i % 3 == 0 ? 1.0 : 0.997 + 0.001 * (i % 7)});
    }

    std::vector<unsigned char> cached_scalar;
    for (bool cached : {true, false}) {
        setenv("NPS_FAST_PULL_CACHE_BUDGET_BYTES", cached ? "1073741824" : "0", 1);
        FastObjectiveEvaluator evaluator(events, data_mpi0, data_mmiss, data_mpgg2,
                                         nsmear, Config::RESOLUTION_A_DEFAULT,
                                         Config::RESOLUTION_B_DEFAULT, Config::RESOLUTION_C_DEFAULT);
        if (evaluator.usesPullCache() != cached) return 3;
        std::vector<unsigned char> scalar_bytes, batch_bytes;
        for (const auto& p : points) {
            FastObjectiveSnapshot s;
            evaluator.evaluateBreakdown(p[0], p[1], p[2], p[3], p[4], p[5], &s);
            if (!(s.breakdown.mpi0.sim_integral > 0.0) || !(s.breakdown.mmiss.sim_integral > 0.0)) {
                std::cerr << "Empty synthetic histogram: mpi0=" << s.breakdown.mpi0.sim_integral
                          << " mmiss=" << s.breakdown.mmiss.sim_integral << "\n";
                return 4;
            }
            appendSnapshot(scalar_bytes, s);
        }
        std::vector<FastObjectiveSnapshot> snapshots;
        evaluator.evaluateBreakdownBatch(points, &snapshots);
        for (const auto& s : snapshots) appendSnapshot(batch_bytes, s);
        if (scalar_bytes != batch_bytes) {
            std::cerr << "Scalar/batch bitwise mismatch\n";
            return 5;
        }
        if (cached) cached_scalar = scalar_bytes;
        else if (cached_scalar != scalar_bytes) {
            std::cerr << "Cached/uncached bitwise mismatch\n";
            return 6;
        }
        for (bool batch : {false, true}) {
            const std::string mode = std::string(cached ? "cached" : "uncached") + (batch ? "_batch" : "_scalar");
            const std::string file = std::string(argv[1]) + "_" + mode + ".bin";
            std::ofstream out(file, std::ios::binary);
            const auto& bytes = batch ? batch_bytes : scalar_bytes;
            out.write(reinterpret_cast<const char*>(bytes.data()), bytes.size());
            out.close();
            volatile double sink = 0.0;
            std::vector<double> durations;
            for (int repeat = 0; repeat < repetitions; ++repeat) {
                const auto start = std::chrono::steady_clock::now();
                if (batch) {
                    const auto values = evaluator.evaluateSelectedBatch(points);
                    for (double v : values) sink += v;
                } else {
                    for (const auto& p : points)
                        sink += evaluator.evaluateSelected(p[0], p[1], p[2], p[3], p[4], p[5]);
                }
                durations.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
            }
            std::sort(durations.begin(), durations.end());
            std::cout << std::setprecision(8) << mode << " events=" << count << " nsmear=" << nsmear
                      << " points=" << points.size() << " median_seconds=" << durations[durations.size()/2]
                      << " cache_bytes=" << evaluator.pullCacheBytes() << " bytes=" << bytes.size()
                      << " sink=" << sink << "\n";
        }
    }
    return 0;
}
