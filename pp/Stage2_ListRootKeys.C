void Stage2_ListRootKeys(char *jetFileName, char *ueFileName)
{
  std::cout << std::endl;
  std::cout << "===== JET FILE =====" << std::endl;

  TFile *jetFile = TFile::Open(jetFileName, "READ");

  if (!jetFile || jetFile->IsZombie()) {
    std::cout << "ERROR: cannot open jet file" << std::endl;
  } else {
    jetFile->ls();
    jetFile->Close();
  }

  std::cout << std::endl;
  std::cout << "===== UE FILE =====" << std::endl;

  TFile *ueFile = TFile::Open(ueFileName, "READ");

  if (!ueFile || ueFile->IsZombie()) {
    std::cout << "ERROR: cannot open UE file" << std::endl;
  } else {
    ueFile->ls();
    ueFile->Close();
  }
}
