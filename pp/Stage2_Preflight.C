// Read-only STAR-runtime preflight for Stage 2.  This macro discovers the
// actual StJetEvent and UE class interfaces in the user's STAR release before
// the baseline reader is compiled.  It prevents us from guessing the APIs
// that associate a particle with a jet and an off-axis cone with that jet.
//
// Run from this directory with, for example:
// root4star -b -q 'Stage2_Preflight.C("../pp/Input/pt5_7_16092051_1.jets.root","../pp/Input/pt5_7_16092051_1.ueoc.root")'
// It creates Stage2_StarApiReport.txt in the current working directory.

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClass.h"
#include "TList.h"
#include "TSystem.h"
#include "TROOT.h"
#include <fstream>

namespace Stage2Preflight {
void Load(const char* library, std::ostream& out) {
  const Int_t rc = gSystem->Load(library);
  out << "LOAD " << library << " : " << rc << "\n";
}

void PrintClassMethods(const char* className, std::ostream& out) {
  TClass* cl = TClass::GetClass(className, kTRUE, kTRUE);
  out << "\nCLASS " << className << "\n";
  if (!cl) { out << "  DICTIONARY NOT AVAILABLE\n"; return; }
  const char* libraries = cl->GetSharedLibs();
  out << "  library=" << (libraries ? libraries : "<not reported>") << "\n";
  TList* methods = cl->GetListOfMethods();
  if (!methods) { out << "  no methods reported\n"; return; }
  // ROOT 5/CINT can enumerate these as TObject even when TMethod's CINT
  // dictionary is unavailable.  The title carries the signature.
  TIter next(methods);
  TObject* method = 0;
  while ((method = next()))
    out << "  " << method->GetName() << " : " << method->GetTitle() << "\n";
}

void PrintTree(const char* label, TFile& file, const char* treeName, std::ostream& out) {
  TTree* tree = dynamic_cast<TTree*>(file.Get(treeName));
  out << "\n" << label << " tree=" << treeName << "\n";
  if (!tree) { out << "  MISSING\n"; return; }
  out << "  entries=" << tree->GetEntries() << "\n";
  TObjArray* branches = tree->GetListOfBranches();
  for (Int_t i = 0; i < branches->GetEntriesFast(); ++i) {
    TBranch* branch = static_cast<TBranch*>(branches->At(i));
    out << "  BRANCH name=" << branch->GetName()
        << " class=" << branch->GetClassName()
        << " title=" << branch->GetTitle() << "\n";
  }
}
}

void Stage2_Preflight(const char* jetsFile,
                      const char* ueFile,
                      const char* reportFile = "Stage2_StarApiReport.txt")
{
  std::ofstream out(reportFile);
  if (!out) { Error("Stage2_Preflight", "Cannot create %s", reportFile); return; }
  out << "Stage 2 STAR API preflight, revision 4\n";
  out << "ROOT=" << gROOT->GetVersion() << "\n";
  out << "STAR=" << (gSystem->Getenv("STAR") ? gSystem->Getenv("STAR") : "<unset>") << "\n";
  out << "STAR_VERSION=" << (gSystem->Getenv("STAR_VERSION") ? gSystem->Getenv("STAR_VERSION") : "<unset>") << "\n";

  // Match the supplied jet-maker's initialization.  In legacy STAR ROOT5
  // releases loadMuDst.C also extends the dynamic-library search path.
  gROOT->Macro("loadMuDst.C");
  gROOT->Macro("LoadLogger.C");
  out << "DYNAMIC_PATH=" << gSystem->GetDynamicPath() << "\n";

  // These four libraries are sufficient to read the persisted objects.  Do
  // not load StJetFinder/StJetMaker here: their old FastJet/DetectorDb ABI is
  // incompatible with SL24c, although the event dictionaries remain usable.
  Stage2Preflight::Load("StJetSkimEvent", out);
  Stage2Preflight::Load("StJets", out);
  Stage2Preflight::Load("StJetEvent", out);
  Stage2Preflight::Load("StUeEvent", out);

  TFile jetInput(jetsFile, "READ");
  TFile ueInput(ueFile, "READ");
  Stage2Preflight::PrintTree("JET", jetInput, "jet", out);
  Stage2Preflight::PrintTree("UE", ueInput, "ue", out);

  const char* classes[] = {
    "StJetEvent", "StJetCandidate", "StJetParticle",
    "StUeEvent", "StUeOffAxisConesEvent", "StUeOffAxisConesJet"
    "StUeEvent", "StUeVertex", "StUeJet", "StUeParticle",
    "StUeOffAxisConesEvent", "StUeOffAxisConesJet", "StUeOffAxisCones"
  };
  for (UInt_t i = 0; i < sizeof(classes)/sizeof(classes[0]); ++i)
    Stage2Preflight::PrintClassMethods(classes[i], out);
  out.close();
  printf("Wrote %s\n", reportFile);
}
