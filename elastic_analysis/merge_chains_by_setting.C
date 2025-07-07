void merge_chains_by_setting() {
  const std::string baseDir = "./data/hms_elastics_skimmed/";
  const std::string treeName = "T";
  const std::string outputDir = "./data/merged_elastics/";

  gSystem->mkdir(outputDir.c_str(), true);  // create output dir if it doesn't exist

  std::map<int, std::vector<std::string>> setting_runs = {
    {1, {"1249"}},
    {2, {"1250"}},
    {3, {"1251"}},
    {4, {"1252"}},
    {5, {"1253"}},
    {6, {"1534"}},
    {7, {"1535"}},
    {8, {"1536"}},
    {9, {"1714"}},
    {10, {"1715", "1716"}},
    {11, {"6828"}}
  };

  // Add runs 6828–6840 to setting 10
  for (int run = 6828; run <= 6840; ++run)
    setting_runs[11].push_back(std::to_string(run));

  for (const auto& [setting, runs] : setting_runs) {
    TChain chain(treeName.c_str());
    for (const auto& run : runs) {
      std::string filepath = baseDir + "nps_hms_coin_skimmed_" + run + ".root";
      if (gSystem->AccessPathName(filepath.c_str())) {
        std::cout << "⚠️  Missing file: " << filepath << std::endl;
        continue;
      }
      chain.Add(filepath.c_str());
    }

    if (chain.GetEntries() == 0) {
      std::cout << "⚠️  No entries found for setting " << setting << ", skipping write.\n";
      continue;
    }

    std::string outputFile = outputDir + "merged_setting_" + std::to_string(setting) + ".root";
    chain.Merge(outputFile.c_str());
    std::cout << "✅ Merged setting " << setting << " into " << outputFile << "\n";
  }

  std::cout << "\n🎉 All available chains have been merged.\n";
}
