// Split the four-column Stage 2 manifest into independent batch-job shards.
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
#include "TSystem.h"

void Stage2_SplitManifest(const char* input="Stage2_pp_inputs_absolute.manifest",
                          const char* outDir="Stage2_shards",
                          Int_t pairsPerShard=100)
{
  if (pairsPerShard <= 0) { std::cout << "ERROR: pairsPerShard must be positive" << std::endl; return; }
  std::ifstream in(input); if (!in) { std::cout << "ERROR: cannot open " << input << std::endl; return; }
  gSystem->mkdir(outDir, kTRUE);
  std::ofstream list((std::string(outDir)+"/shard_list.list").c_str());
  std::string line; Int_t n=0, shard=0; std::ofstream out;
  while (std::getline(in,line)) {
    if (line.size()==0 || line[0]=='#') continue;
    if (n % pairsPerShard == 0) {
      if (out.is_open()) out.close();
      std::ostringstream name; name << outDir << "/Stage2_pp_shard_" << std::setfill('0') << std::setw(4) << shard << ".manifest";
      out.open(name.str().c_str()); out << "# pthat_min pthat_max jet_file ue_file\n";
      list << name.str() << "\n"; ++shard;
    }
    out << line << "\n"; ++n;
  }
  if (out.is_open()) out.close(); list.close(); in.close();
  std::cout << "Wrote " << shard << " shards containing " << n << " file pairs in " << outDir << std::endl;
}

