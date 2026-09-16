//
// Pibero Djawotho <pibero@tamu.edu>
// Texas A&M University
// 11 September 2009
//

void RunPythia_Perugia2012(const int nevents = 100, const float ptmin = 9, const float ptmax = 11, const int seed = 1, const char* outfile = "pythia.root")
{
  gSystem->Load("libStJetSkimEvent.so");
  assert(gSystem->Load("/star/u/tinglin/lib/lib/libLHAPDF.so") == 0);
  assert(gSystem->Load("/star/u/tinglin/lib/pythia6/libPythia6-6.4.28.so") == 0);

  // Create an instance of Pythia
  TPythia6* pythia = new TPythia6;

  printf("Seed = %d\n",seed);

  // Physics processes
  pythia->SetMSEL(1);		// QCD jets
  //  printf("Seed is %d\n",(seed+1));
  pythia->SetMRPY(1,seed);      // Seed random number generator
  pythia->SetCKIN(3,ptmin);	// Lower partonic pT bound in GeV
  pythia->SetCKIN(4,ptmax);	// Higher partonic pT bound in GeV
  pythia->SetMSTP(5,370);	// Perugia 2012 tune

  //pythia->SetMSTJ(1,0); // Turn Off Fragmentation

  printf("CKIN[3] = %f\n",pythia->GetCKIN(3));
  printf("CKIN[4] = %f\n",pythia->GetCKIN(4));

  // Make the following stable
  pythia->SetMDCY(102,1,0);  // PI0 111
  pythia->SetMDCY(106,1,0);  // PI+ 211
  pythia->SetMDCY(109,1,0);  // ETA 221
  pythia->SetMDCY(116,1,0);  // K+ 321
  pythia->SetMDCY(112,1,0);  // K_SHORT 310
  pythia->SetMDCY(105,1,0);  // K_LONG 130
  pythia->SetMDCY(164,1,0);  // LAMBDA0 3122
  pythia->SetMDCY(167,1,0);  // SIGMA0 3212
  pythia->SetMDCY(162,1,0);  // SIGMA- 3112
  pythia->SetMDCY(169,1,0);  // SIGMA+ 3222
  pythia->SetMDCY(172,1,0);  // Xi- 3312
  pythia->SetMDCY(174,1,0);  // Xi0 3322
  pythia->SetMDCY(176,1,0);  // OMEGA- 3334

  pythia->Initialize("cms","p","p",200); // p+p collisions at sqrt(s)=200 GeV
  //  pythia->Initialize("cms","p","p",500); // p+p collisions at sqrt(s)=500 GeV

  // Choice of PDF
//  pythia->SetMSTP(51,7);
//  pythia->SetMSTP(52,1);
  pythia->SetPARP(90,0.213); // 2012 pp500 official usage

  //printf("PARP[91] = %f\n",pythia->GetPARP(91));

  TFile* ofile = TFile::Open(outfile,"recreate");
  assert(ofile);

  StPythiaEvent* event = new StPythiaEvent;
  TTree* tree = new TTree("PythiaTree","Pythia Record");
  tree->Branch("PythiaBranch","StPythiaEvent",&event);

  TClonesArray* particles = new TClonesArray("TParticle");

  // Event loop
  for (int iEvent = 1; iEvent <= nevents; ++iEvent) {
    pythia->GenerateEvent();
    if (iEvent % 10000 == 0) pythia->Pylist(1);
    pythia->ImportParticles(particles,"All"); // Store all particles

    event->setRunId(1);
    event->setEventId(iEvent);
    event->setProcessId(pythia->GetMSTI(1));
    event->setS(pythia->GetPARI(14));
    event->setT(pythia->GetPARI(15));
    event->setU(pythia->GetPARI(16));
    event->setPt(pythia->GetPARI(17));
    event->setCosTheta(pythia->GetPARI(41));
    event->setX1(pythia->GetPARI(33));
    event->setX2(pythia->GetPARI(34));
    event->setVertex(TVector3(0,0,0));
    event->setMstu72(pythia->GetMSTU(72));
    event->setMstu73(pythia->GetMSTU(73));
    event->setMstp111(pythia->GetMSTP(111));

    for (int i = 0; i < particles->GetEntriesFast(); ++i)
      event->addParticle(*(TParticle*)particles->At(i));

    tree->Fill();

    event->Clear();
  } // End event loop

  pythia->Pystat(1);		// Print PYTHIA statistics

  printf("PARI(1) = %e\n",pythia->GetPARI(1));

  pythia->Pylist(1);

  ofile->Write();
  ofile->Close();
}
