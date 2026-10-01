#pragma once

// Shared ROOT/STL dependencies. No event state or objects are created here.

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TGraph.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TGaxis.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TDecompLU.h>
#include <TObjString.h>
#include <TParameter.h>
#include <TMath.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#ifdef NPS_ENABLE_PARTONS
#include "partons_pi0_projection.h"
#endif

namespace fs = std::filesystem;
