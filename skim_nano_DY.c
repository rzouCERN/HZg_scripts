// Load the library at macro parsing time: we need this to use its content in the code
//R__LOAD_LIBRARY($ROOTSYS/test/libEvent.so)
#include <string>

void skim_nano_DY(const char* in_filename)
{
 
  vector<Long64_t> match_list;
  vector<Long64_t> nomatch_list;
  string line;
  ifstream ReadFile("/net/cms37/data1/xf82/jet_ph_match/el_match.txt");
  while (getline(ReadFile, line)){
    if (line[0] == 'e') {
      line = line.substr(11,  line.size()-11);
      Long64_t numList = static_cast<Long64_t>(stoi(line));
      match_list.push_back(numList);
    }
  }
  ReadFile.close();

  ifstream ReadFileNo("/net/cms37/data1/xf82/jet_ph_match/el_nomatch.txt");
  while (getline(ReadFileNo, line)){
    if (line[0] == 'e') {
      line = line.substr(11,  line.size()-11);
      Long64_t numList = static_cast<Long64_t>(stoi(line));
      nomatch_list.push_back(numList);
    }
  }
  ReadFileNo.close();
  vector<Long64_t> match_list_mu;
  vector<Long64_t> nomatch_list_mu;

  ifstream ReadFile_mu("/net/cms37/data1/xf82/jet_ph_match/mu_match.txt");
  while (getline(ReadFile_mu, line)){
    if (line[0] == 'e') {
      line = line.substr(11,  line.size()-11);
      Long64_t numList = static_cast<Long64_t>(stoi(line));
      match_list_mu.push_back(numList);
    }
  }
  ReadFile_mu.close();

  ifstream ReadFileNo_mu("/net/cms37/data1/xf82/jet_ph_match/mu_nomatch.txt");
  while (getline(ReadFileNo_mu, line)){
    if (line[0] == 'e') {
      line = line.substr(11,  line.size()-11);
      Long64_t numList = static_cast<Long64_t>(stoi(line));
      nomatch_list_mu.push_back(numList);
    }
  }
  ReadFileNo_mu.close();

  cout << "L30" << endl;
  TString dir = "/net/cms17/cms17r0/pico/NanoAODv9/nano/2017/mc/";
  //  TString dir = "/net/cms17/cms17r0/pico/NanoAODv2/nano/2017/signal/";
  //  gSystem->ExpandPathName(dir);
  //const auto filename = gSystem->AccessPathName(dir) ? "./Event.root" : "$ROOTSYS/test/Event.root";
  TString filename=in_filename;// = "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v2__2430000__FDE1C172-542B-3244-8611-978D53C5542D.root";

  TFile oldfile(dir+filename);
  TTree *oldtree = (TTree *) oldfile.Get("Events");
  //  oldfile.GetObject("Events", oldtree);

  float Photon_pt[20];
  UInt_t nPhoton;
  UInt_t nElectron;
  UInt_t  nMuon;
  //  oldtree->SetBranchAddress("Photon_pt",&Photon_pt);
  //  oldtree->SetBranchAddress("nPhoton",&nPhoton);
 
  // Deactivate all branches
  oldtree->SetBranchStatus("*", 0);

  // Activate only four of them
  //  for (auto activeBranchName : {"nElectron","nMuon","nPhoton","Generator_weight","HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL","Electron_dz","Electron_dxy","Muon_dz","Muon_dxy","Photon_pt","Photon_eta","Photon_phi","Photon_mvaID","Photon_energyErr","Photon_pfRelIso03_all","Photon_isScEtaEB","Photon_isScEtaEE","Photon_electronVeto","Electron_charge","Electron_pt","Electron_eta","Electron_phi","Electron_deltaEtaSC","Electron_mvaFall17V2Iso_WPL","Electron_energyErr","Muon_charge","Muon_pt","Muon_eta","Muon_phi","Muon_looseId","Muon_highPtId","Muon_pfRelIso03_all","Muon_sip3d"})
  ULong64_t event;
  cout << "L54" << endl;
  for (auto activeBranchName : {"event","nElectron","nMuon","nPhoton","Photon_pt","nJet","Jet_area","Jet_eta","Jet_mass","Jet_phi","Jet_pt","Jet_qgl","Jet_electronIdx1","Jet_electronIdx2","Jet_jetId","Jet_muonIdx1","Jet_muonIdx2","Jet_nConstituents","Jet_nElectrons","Jet_nMuons","Jet_puId","Jet_puIdDisc","Jet_genJetIdx","Jet_hadronFlavour","Jet_partonFlavour","Jet_cleanmask"})
      oldtree->SetBranchStatus(activeBranchName, 1);
  oldtree->SetBranchAddress("nPhoton",&nPhoton);
  oldtree->SetBranchAddress("Photon_pt",&Photon_pt);
  oldtree->SetBranchAddress("nElectron",&nElectron);
  oldtree->SetBranchAddress("nMuon",&nMuon);
  oldtree->SetBranchAddress("event",&event);
  cout << "L61" << endl;
 
  // Create a new file + a clone of old tree in new file
  TFile newfile("/homes/rui/workspace2/nano_skim/"+filename.ReplaceAll(".root","_skim_DY.root"), "recreate");
  auto newtree = oldtree->CloneTree(0);
  bool match_event;
  bool is_el;
  TBranch *match_branch = newtree->Branch("match",&match_event,"match/B");
  TBranch *is_el_branch = newtree->Branch("is_el",&is_el,"is_el/B");
  cout << "L68" << endl;

  Int_t nevents = (Int_t)oldtree->GetEntries();
  //  Int_t nevents = oldtree->GetEntries(); 
  for(int ievent=0; ievent<nevents; ievent++)
    {
      oldtree->GetEntry(ievent);
      if (nPhoton < 1 || (nElectron < 2 && nMuon < 2)) continue;
      //      cout << "L74" << endl;
      float maxpt = 0;
      //      cout << nPhoton << " "<<  nElectron << " " << nMuon << endl;
      /* for (int i(0); i<nPhoton; i++){ */
      /* 	if (Photon_pt[i] > maxpt) maxpt = Photon_pt[i];  */
      /* } */
      /* if( maxpt <= 15)	continue; */
      //      cout << "L78" << endl;
      int mat =0;
      int nomat =0;
      int mat_mu =0;
      int nomat_mu =0;
      //      cout << " ievent " << ievent << " event: " << unsigned(event) << endl;

      for (auto entry: nomatch_list){
	//	cout << entry << " event: " << event << endl;
	if (entry==event) nomat =1;
      }
      //      cout << "L85" << endl;
      if (nomat !=1){
	for (auto entry: nomatch_list_mu){
	  if (entry==event) nomat_mu =1;
	}
	if (nomat_mu !=1){
	  for (auto entry: match_list){ 
	    if (entry==event) mat =1;
	  }
	  if (mat != 1) {
	    for (auto entry: match_list_mu){
	      if (entry==event) mat_mu =1;
	    }
	    if (mat_mu !=1){
	      continue;
	    } else {
	      match_event = true;
	      is_el = false;
	    }
	  } else {
	    match_event = true;
	    is_el = true;
	  }
	} else {
	  match_event = false;
	  is_el = false;
	}
      } else {
	match_event = false;
	is_el = true;
      }
      
      //      match_branch->Fill();
      newtree->Fill();
    }
 
  newtree->Print();
  newfile.Write();
}

