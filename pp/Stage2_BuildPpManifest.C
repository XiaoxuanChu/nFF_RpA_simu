// Build a read-only manifest of all matched pp particle-jet/UE file pairs.
// The expected storage layout is
//   <base>/<run>/pt<min>_<max>_<run>_<chunk>.jets.root
//   <base>/<run>/pt<min>_<max>_<run>_<chunk>.ueoc.root
// with one run number per line in filelist.list.
//
// Example (run from the directory holding filelist.list):
// root4star -b -q 'Stage2_BuildPpManifest.C("/star/u/tinglin/Run15/2015pp200Embedding/pythia","filelist.list","Stage2_pp_inputs.manifest")'
//
// Output columns: pthat_min pthat_max jet_file ue_file

#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"
#include "TString.h"
#include <fstream>
#include <cstdio>
#include <string>

void Stage2_BuildPpManifest(const char* baseDirectory,
                            const char* runListFile = "filelist.list",
                            const char* outputManifest = "Stage2_pp_inputs.manifest")
{
  std::ifstream runList(runListFile);
  if (!runList) { Error("Stage2_BuildPpManifest","Cannot open %s",runListFile); return; }
  std::ofstream output(outputManifest);
  if (!output) { Error("Stage2_BuildPpManifest","Cannot create %s",outputManifest); return; }
  output << "# pthat_min pthat_max jet_file ue_file\n";
  const Int_t nHardBins = 11;
  const Int_t hardMin[11] = {2,3,4,5,7,9,11,15,20,25,35};
  const Int_t hardMax[11] = {3,4,5,7,9,11,15,20,25,35,-1};
  Long64_t nPairs[11];
  for (Int_t i=0; i<nHardBins; ++i) nPairs[i] = 0;
  Long64_t totalPairs = 0, missingUe = 0;
  std::string run;
  while (std::getline(runList,run)) {
    if (run.empty() || run[0] == '#') continue;
    TString runDirectory = TString(baseDirectory) + "/" + run.c_str();
    TSystemDirectory directory(runDirectory,runDirectory);
    TList* files = directory.GetListOfFiles();
    if (!files) { Warning("Stage2_BuildPpManifest","Cannot list %s",runDirectory.Data()); continue; }
    TIter next(files);
    TSystemFile* file = 0;
    while ((file = static_cast<TSystemFile*>(next()))) {
      if (file->IsDirectory()) continue;
      TString name = file->GetName();
      if (!name.BeginsWith("pt") || !name.EndsWith(".jets.root")) continue;
      Int_t ptMin = 0, ptMax = 0;
      if (sscanf(name.Data(),"pt%d_%d_",&ptMin,&ptMax) != 2) {
        Warning("Stage2_BuildPpManifest","Cannot parse pThat bin from %s",name.Data());
        continue;
      }
      TString jetFile = runDirectory + "/" + name;
      TString ueName = name;
      ueName.ReplaceAll(".jets.root",".ueoc.root");
      TString ueFile = runDirectory + "/" + ueName;
      if (gSystem->AccessPathName(ueFile)) {
        Warning("Stage2_BuildPpManifest","Missing UE partner for %s",jetFile.Data());
        ++missingUe;
        continue;
      }
      Int_t hardBin = -1;
      for (Int_t i=0; i<nHardBins; ++i)
        if (ptMin == hardMin[i] && ptMax == hardMax[i]) hardBin = i;
      if (hardBin < 0) {
        Warning("Stage2_BuildPpManifest","pThat bin %d_%d is absent from locked pp luminosity inputs",ptMin,ptMax);
        continue;
      }
      ++nPairs[hardBin];
      output << ptMin << " " << ptMax << " " << jetFile.Data() << " " << ueFile.Data() << "\n";
      ++totalPairs;
    }
    delete files;
  }
  output.close();
  printf("Wrote %s with %lld matched jet/UE pairs; missing UE partners=%lld\n",outputManifest,totalPairs,missingUe);
  for (Int_t i=0; i<nHardBins; ++i)
    printf("  pThat %d_%d : %lld pairs\n",hardMin[i],hardMax[i],nPairs[i]);
}

