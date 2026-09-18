#include "TSystem.h"
#include "TString.h"
#include "Stage2_pp_Baseline.C"
#include <iostream>
void Stage2_BatchShard() {
  const char* manifest = gSystem->Getenv("STAGE2_MANIFEST");
  const char* tag = gSystem->Getenv("STAGE2_TAG");
  if (!manifest || !tag) { std::cout << "ERROR: missing STAGE2_MANIFEST or STAGE2_TAG" << std::endl; return; }
  Stage2_pp_Baseline(manifest, Form("partial_%s.root",tag), Form("partial_%s.txt",tag), 0);
}
