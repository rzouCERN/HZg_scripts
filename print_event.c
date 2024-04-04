// Load the library at macro parsing time: we need this to use its content in the code
//R__LOAD_LIBRARY($ROOTSYS/test/libEvent.so)
#include <string>

void print_event(const char* in_filename)
{
 
  cout << "L30" << endl;
  //  TString dir = "/net/cms17/cms17r0/pico/NanoAODv2/nano/2017/signal/";
  //  gSystem->ExpandPathName(dir);
  //const auto filename = gSystem->AccessPathName(dir) ? "./Event.root" : "$ROOTSYS/test/Event.root";
  TString dir = "/net/cms37/data2/rui/";
  TString filename=in_filename;// = "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v2__2430000__FDE1C172-542B-3244-8611-978D53C5542D.root";

  TFile oldfile(dir+filename);
  TTree *oldtree = (TTree *) oldfile.Get("Events");
  //  oldfile.GetObject("Events", oldtree);

  // Deactivate all branches
  oldtree->SetBranchStatus("*", 0);

  // Activate only four of them
  //  for (auto activeBranchName : {"nElectron","nMuon","nPhoton","Generator_weight","HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL","Electron_dz","Electron_dxy","Muon_dz","Muon_dxy","Photon_pt","Photon_eta","Photon_phi","Photon_mvaID","Photon_energyErr","Photon_pfRelIso03_all","Photon_isScEtaEB","Photon_isScEtaEE","Photon_electronVeto","Electron_charge","Electron_pt","Electron_eta","Electron_phi","Electron_deltaEtaSC","Electron_mvaFall17V2Iso_WPL","Electron_energyErr","Muon_charge","Muon_pt","Muon_eta","Muon_phi","Muon_looseId","Muon_highPtId","Muon_pfRelIso03_all","Muon_sip3d"})
  ULong64_t event;
  cout << "L54" << endl;
  oldtree->SetBranchStatus("event",1);
  oldtree->SetBranchAddress("event",&event);
  cout << "L61" << endl;
 
  Int_t nevents = (Int_t)oldtree->GetEntries();
  //  Int_t nevents = oldtree->GetEntries(); 
  for(int ievent=0; ievent<nevents; ievent++)
    {
      oldtree->GetEntry(ievent);
      cout << event << endl;
    }
}

