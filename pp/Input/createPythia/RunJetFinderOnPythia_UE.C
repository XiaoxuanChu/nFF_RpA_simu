//
// Pibero Djawotho <pibero@tamu.edu>
// Texas A&M
// 25 Nov 2012
//

void RunJetFinderOnPythia_UE(int nevents = 100000,
			  const char* pythiafile = "pythia.root",
                          const char* uefile = "pythia.ueoc.root",
			  const char* jetfile = "pythia.jets.root")
{
  // Load shared libraries
  gROOT->Macro("loadMuDst.C");
  gROOT->Macro("LoadLogger.C");

  gSystem->Load("StDetectorDbMaker");
  gSystem->Load("StTpcDb");
  gSystem->Load("StDbUtilities");
  gSystem->Load("StMcEvent");
  gSystem->Load("StMcEventMaker");
  gSystem->Load("StDaqLib");
  gSystem->Load("StEmcRawMaker");
  gSystem->Load("StEmcADCtoEMaker");
  gSystem->Load("StEpcMaker");
  gSystem->Load("StEmcSimulatorMaker");
  gSystem->Load("StDbBroker");
  gSystem->Load("St_db_Maker");
  gSystem->Load("StEEmcUtil");
  gSystem->Load("StEEmcDbMaker");
  gSystem->Load("StSpinDbMaker");
  gSystem->Load("StEmcTriggerMaker");
  gSystem->Load("StTriggerUtilities");
  gSystem->Load("StMCAsymMaker");
  gSystem->Load("StRandomSelector");

  gSystem->Load("libfastjet.so");
  gSystem->Load("libsiscone.so");
  gSystem->Load("libsiscone_spherical.so");
  gSystem->Load("libfastjetplugins.so");

  gSystem->Load("StJetFinder");
  gSystem->Load("StJetSkimEvent");
  gSystem->Load("StJets");
  gSystem->Load("StJetEvent");
  gSystem->Load("StUeEvent");
  gSystem->Load("StJetMaker");
  gSystem->Load("StTriggerFilterMaker");

  cout << "Loading shared libraries done" << endl;

  // Create chain
  StChain* chain = new StChain;

  // STAR database: Dummy. Not used for PYTHIA jets, but required by StJetMaker.
  St_db_Maker* starDb = new St_db_Maker("StarDb","MySQL:StarDb");
  starDb->SetDateTime(20090628,53220); // Run 10179006

  // Endcap database
  StEEmcDbMaker* eemcDb = new StEEmcDbMaker;

  // Read Pythia records
  St_pythia_Maker* pythia = new St_pythia_Maker;
  pythia->SetFile(pythiafile);

  // Pythia4pMaker
  StPythiaFourPMaker* pythia4pMaker = new StPythiaFourPMaker;

  // Instantiate the JetMaker
  StJetMaker2012* jetmaker = new StJetMaker2012;
  jetmaker->setJetFile(jetfile);
  jetmaker->setJetFileUe(uefile);
  //------------------------------------------------------------------------------------

  // Set analysis cuts for particle jets branch
  StAnaPars* anaparsParticle = new StAnaPars;
  anaparsParticle->useMonteCarlo = true;

  // MC cuts
  anaparsParticle->addMcCut(new StjMCParticleCutStatus(1)); // final state particles

  // Jet cuts
  anaparsParticle->addJetCut(new StProtoJetCutPt(6,100));
  anaparsParticle->addJetCut(new StProtoJetCutEta(-0.8,1.8));

  // Set analysis cuts for parton jets branch
  StAnaPars* anaparsParton = new StAnaPars;
  anaparsParton->useMonteCarlo = true;

  // MC cuts
  anaparsParton->addMcCut(new StjMCParticleCutParton);

  // Jet cuts
  anaparsParton->addJetCut(new StProtoJetCutPt(6,100));
  anaparsParton->addJetCut(new StProtoJetCutEta(-0.8,1.8));

  StFastJetAreaPars *JetAreaPars = new StFastJetAreaPars;
  // Set anti-kt R=0.6 parameters
  StFastJetPars* AntiKtR060Pars = new StFastJetPars;
  AntiKtR060Pars->setJetAlgorithm(StFastJetPars::antikt_algorithm);
  AntiKtR060Pars->setRparam(0.6);
  AntiKtR060Pars->setRecombinationScheme(StFastJetPars::E_scheme);
  AntiKtR060Pars->setStrategy(StFastJetPars::Best);
  AntiKtR060Pars->setPtMin(6.0);
  AntiKtR060Pars->setJetArea(JetAreaPars);

  jetmaker->addBranch("AntiKtR060Particle",anaparsParticle,AntiKtR060Pars);
  jetmaker->addBranch("AntiKtR060Parton",anaparsParton,AntiKtR060Pars);

  StOffAxisConesPars *off060 = new StOffAxisConesPars(0.6);
  jetmaker->addUeBranch("OffAxisConesR060", off060);
  //------------------------------------------------------------------------------------

  chain->Init();
  chain->EventLoop(nevents);
}
