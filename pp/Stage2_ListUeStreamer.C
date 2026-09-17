#include "TFile.h"
#include "TList.h"
#include "TCollection.h"
#include "TObjArray.h"
#include "TStreamerInfo.h"
#include "TStreamerElement.h"
#include <iostream>
#include <string.h>

void Stage2_ListUeStreamer(const char *ueFileName)
{
  TFile *file = TFile::Open(ueFileName, "READ");

  if (!file || file->IsZombie()) {
    std::cout << "ERROR: cannot open file." << std::endl;
    return;
  }

  TList *infos = file->GetStreamerInfoList();
  TIter nextInfo(infos);
  TStreamerInfo *info = 0;

  while ((info = (TStreamerInfo *)nextInfo())) {
    if (strcmp(info->GetName(), "StUeParticle") != 0) continue;

    std::cout << "Class: " << info->GetName()
              << " version=" << info->GetClassVersion()
              << std::endl;

    TObjArray *elements = info->GetElements();

    for (Int_t i = 0; i < elements->GetEntriesFast(); ++i) {
      TStreamerElement *element =
          (TStreamerElement *)elements->At(i);

      std::cout << "field=" << element->GetName()
                << " type=" << element->GetTypeName()
                << " offset=" << element->GetOffset()
                << std::endl;
    }

    file->Close();
    return;
  }

  std::cout << "ERROR: StUeParticle streamer info not found."
            << std::endl;
  file->Close();
}
