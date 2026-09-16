#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMethodCall.h"
#include "TClass.h"

#include <fstream>
#include <cstdio>
#include <string>

Int_t Stage2ThresholdBin(Int_t lo, Int_t hi) {
  const Int_t low[11]  = {2,3,4,5,7,9,11,15,20,25,35};
  const Int_t high[11] = {3,4,5,7,9,11,15,20,25,35,-1};

  for (Int_t i = 0; i < 11; ++i) {
    if (lo == low[i] && hi == high[i]) return i;
  }

  return -1;
}

void Stage2_ValidateJetThreshold(
  const char* manifestFile,
  const char* outputReport = "Stage2_jet_threshold_report.txt")
{
  gROOT->Macro("loadMuDst.C");
  gROOT->Macro("LoadLogger.C");

  if (gSystem->Load("StJetSkimEvent") < 0 ||
      gSystem->Load("StJets") < 0 ||
      gSystem->Load("StJetEvent") < 0 ||
      gSystem->Load("StUeEvent") < 0) {
    Error("Stage2_ValidateJetThreshold",
          "Required STAR dictionaries are unavailable.");
    return;
  }

  TClass* jetEventClass = TClass::GetClass("StJetEvent");
  TClass* jetClass = TClass::GetClass("StJetCandidate");

  if (!jetEventClass || !jetClass) {
    Error("Stage2_ValidateJetThreshold",
          "Required STAR class dictionaries are unavailable.");
    return;
  }

  TMethodCall numberOfJets;
  TMethodCall jetAt;
  TMethodCall jetPt;

  numberOfJets.InitWithPrototype(
    jetEventClass, "numberOfJets", "");
  jetAt.InitWithPrototype(
    jetEventClass, "jet", "Int_t");
  jetPt.InitWithPrototype(
    jetClass, "pt", "");

  if (!numberOfJets.IsValid() ||
      !jetAt.IsValid() ||
      !jetPt.IsValid()) {
    Error("Stage2_ValidateJetThreshold",
          "Cannot resolve required STAR methods.");
    return;
  }

  std::ifstream manifest(manifestFile);
  std::ofstream report(outputReport);

  if (!manifest || !report) {
    Error("Stage2_ValidateJetThreshold",
          "Cannot open manifest or output report.");
    return;
  }

  const Int_t low[11]  = {2,3,4,5,7,9,11,15,20,25,35};
  const Int_t high[11] = {3,4,5,7,9,11,15,20,25,35,-1};

  Bool_t tested[11];
  Double_t minimum[11];
  Long64_t entries[11];
  Long64_t jets[11];

  for (Int_t i = 0; i < 11; ++i) {
    tested[i] = kFALSE;
    minimum[i] = 1.e30;
    entries[i] = 0;
    jets[i] = 0;
  }

  std::string line;

  while (std::getline(manifest, line)) {
    if (line.empty() || line[0] == '#') continue;

    Int_t lo = 0;
    Int_t hi = 0;
    char jetName[4096];
    char ueName[4096];

    if (sscanf(line.c_str(), "%d %d %4095s %4095s",
               &lo, &hi, jetName, ueName) != 4) {
      continue;
    }

    const Int_t bin = Stage2ThresholdBin(lo, hi);
    if (bin < 0 || tested[bin]) continue;

    TFile input(jetName, "READ");
    TTree* tree = dynamic_cast<TTree*>(input.Get("jet"));

    if (!tree || !tree->GetBranch("AntiKtR060Particle")) {
      report << lo << " " << hi << " ERROR\n";
      tested[bin] = kTRUE;
      continue;
    }

    void* event = 0;
    TBranch* branch = tree->GetBranch("AntiKtR060Particle");
    branch->SetAddress(&event);

    entries[bin] = tree->GetEntries();

    for (Long64_t ie = 0; ie < tree->GetEntries(); ++ie) {
      branch->GetEntry(ie);

      Long_t nJets = 0;
      numberOfJets.Execute(event, nJets);

      for (Long_t ij = 0; ij < nJets; ++ij) {
        jetAt.ResetParam();
        jetAt.SetParam(ij);

        Long_t jetAddress = 0;
        jetAt.Execute(event, jetAddress);

        if (!jetAddress) continue;

        Double_t pt = 0.;
        jetPt.Execute((void*)jetAddress, pt);

        ++jets[bin];
        if (pt < minimum[bin]) minimum[bin] = pt;
      }
    }

    tested[bin] = kTRUE;
  }

  report << "# pthat_min pthat_max entries jets minimum_stored_jet_pt\n";

  for (Int_t i = 0; i < 11; ++i) {
    report << low[i] << " " << high[i] << " "
           << entries[i] << " " << jets[i] << " ";

    if (!tested[i] || minimum[i] > 1.e20) report << "NA";
    else report << minimum[i];

    report << "\n";
  }

  report.close();
  printf("Wrote %s\n", outputReport);
}
