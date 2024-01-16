#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
 
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
 
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TVector3.h"

TLorentzVector AssignL1(vector<int> ll_lepid, vector<int> ll_i1, vector<int> ll_i2, vector<int> el_charge, vector<int> mu_charge, vector<double> el_pt, vector<double> el_eta, vector<double> el_phi,vector<double> mu_pt, vector<double> mu_eta, vector<double> mu_phi){
    TLorentzVector l1;
    int i1 = -1;
    if (ll_lepid[0] == 11){
        if(el_charge[ll_i1[0]] < 0) i1 = ll_i1[0];
        else i1 = ll_i2[0];
        l1.SetPtEtaPhiM(el_pt[i1], el_eta[i1], el_phi[i1], 0.00511);
    
        /*if(el_charge[0] < 0) i1 = 0;
        else i1 = 1;
        l1.SetPtEtaPhiM(el_pt[i1], el_eta[i1], el_phi[i1], 0.00511);*/
    }
    else if (ll_lepid[0] == 13){
        if(mu_charge[ll_i1[0]] < 0) i1 = ll_i1[0]; 
        else i1 = ll_i2[0];
        l1.SetPtEtaPhiM(mu_pt[i1], mu_eta[i1], mu_phi[i1], 0.10565);
	  /*if(mu_charge[0] < 0) i1 = 0; 
        else i1 = 1;
        l1.SetPtEtaPhiM(mu_pt[i1], mu_eta[i1], mu_phi[i1], 0.10565);*/

    }
    //cout<< "l1_id = "<<i1<< endl;
    return l1;
    }

TLorentzVector AssignL2(vector<int> ll_lepid, vector<int> ll_i1, vector<int> ll_i2, vector<int> el_charge, vector<int> mu_charge, vector<double> el_pt, vector<double> el_eta, vector<double> el_phi,vector<double> mu_pt, vector<double> mu_eta, vector<double> mu_phi){
    TLorentzVector l2;
    int i1 = -1;
    if (ll_lepid[0] == 11){
        if(el_charge[ll_i1[0]] < 0) i1 = ll_i2[0];
        else i1 = ll_i1[0];
        l2.SetPtEtaPhiM(el_pt[i1], el_eta[i1], el_phi[i1], 0.00511);
        /*if(el_charge[0] < 0) i1 = 1;
        else i1 = 0;
        l2.SetPtEtaPhiM(el_pt[i1], el_eta[i1], el_phi[i1], 0.00511);*/
    }
    else if (ll_lepid[0] == 13){
        if(mu_charge[ll_i1[0]] < 0) i1 = ll_i2[0];
        else i1 = ll_i1[0];
        l2.SetPtEtaPhiM(mu_pt[i1], mu_eta[i1], mu_phi[i1], 0.10565);
        /*if(mu_charge[0] < 0) i1 = 1;
        else i1 = 0;
        l2.SetPtEtaPhiM(mu_pt[i1], mu_eta[i1], mu_phi[i1], 0.10565);*/
    }
    //cout<< "l2_id = "<<i1<<endl;
    return l2;
    }

TLorentzVector AssignZ(vector<double> ll_pt, vector<double> ll_eta, vector<double> ll_phi, vector<double> ll_m){
    TLorentzVector Z;
    Z.SetPtEtaPhiM(ll_pt[0], ll_eta[0], ll_phi[0], ll_m[0]);
    return Z;
}

TLorentzVector AssignGamma(vector<double> photon_pt, vector<double> photon_eta, vector<double> photon_phi, vector<int>llphoton_iph){
    TLorentzVector g;
    g.SetPtEtaPhiM(photon_pt[llphoton_iph[0]],photon_eta[llphoton_iph[0]], photon_phi[llphoton_iph[0]],0);
    return g;
}

TLorentzVector AssignH(vector<double> llphoton_pt, vector<double> llphoton_eta, vector<double> llphoton_phi, vector<double> llphoton_m){
    TLorentzVector h;
    h.SetPtEtaPhiM(llphoton_pt[0],llphoton_eta[0], llphoton_phi[0],llphoton_m[0]);
    return h;
}

TLorentzVector AssignQ1(TLorentzVector h){
    TVector3 htran = h.BoostVector();
    htran.SetZ(0);
    h.Boost(-1*htran);
    TLorentzVector k1;
    double pz, E;
    pz = h.Pz() + h.E();
    E  = h.E()  + h.Pz();
    k1.SetPxPyPzE(0,0,pz/2,E/2);
    k1.Boost(htran);
    return k1;
}

TLorentzVector AssignQ2(TLorentzVector h) {
    TLorentzVector k2;
    TVector3 htran = h.BoostVector();
    htran.SetZ(0);
    h.Boost(-1*htran);
    double pz, E;
    pz = h.Pz() - h.E();
    E  = h.E()  - h.Pz();
    k2.SetPxPyPzE(0,0,pz/2,E/2);
    k2.Boost(htran);
    return k2;
  }

double lambdaZ(TLorentzVector h, TLorentzVector z){
    double M = h.M();
    double mll = z.M();
    return sqrt(pow(h.Dot(z)/M, 2) - pow(mll,2));
}

double GetmindR (TLorentzVector l1, TLorentzVector l2, TLorentzVector g){
    double dR1 = l1.DeltaR(g);
    double dR2 = l2.DeltaR(g);
    if (dR1>dR2) return dR2;
    else return dR1;
}

double GetmaxdR (TLorentzVector l1, TLorentzVector l2, TLorentzVector g){
    double dR1 = l1.DeltaR(g);
    double dR2 = l2.DeltaR(g);
    if (dR1>dR2) return dR1;
    else return dR2;
}


void producer_Input(const char* in_name, const char* output_dir, const char* input_dir, const char* in_year) {
//Open files
  TString filename_pico = "";
  TString filename_out = output_dir;
  int stitch = 0;
  if (strcmp( in_name, "DY")==0) {
    stitch = 1; //DY
    if (strcmp( in_year, "2017")==0)
      filename_pico = "/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/2017/mc/merged_zgmc_llg/merged_pico_llg_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_zgmc_llg_nfiles_153.root";
    else if (strcmp( in_year, "2016")==0)
      filename_pico = "/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/2016/mc/merged_zgmc_llg/merged_pico_llg_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_zgmc_llg_nfiles_41.root";
    else if (strcmp( in_year, "2018")==0)
      filename_pico = "/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/2018/mc/merged_zgmc_llg/merged_pico_llg_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_zgmc_llg_nfiles_204.root";
    else if (strcmp( in_year, "2016APV")==0)
      filename_pico = "/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/2016APV/mc/merged_zgmc_llg/merged_pico_llg_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_zgmc_llg_nfiles_62.root";

    //    filename_pico = "/net/cms37/data1/xf82/BDT/Output/DY_rui_output.root";
  }  else if (strcmp( in_name, "SMZg")==0) {
    stitch = 2;
    if (strcmp( in_year, "2017")==0)
      filename_pico = "/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/2017/mc/merged_zgmc_llg/merged_pico_llg_ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_zgmc_llg_nfiles_45.root";
    else if (strcmp( in_year, "2016")==0)
      filename_pico = "/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/2016/mc/merged_zgmc_llg/merged_pico_llg_ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_zgmc_llg_nfiles_26.root";
    else if (strcmp( in_year, "2018")==0)
      filename_pico = "/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/2018/mc/merged_zgmc_llg/merged_pico_llg_ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_zgmc_llg_nfiles_55.root";
    else if (strcmp( in_year, "2016APV")==0)
      filename_pico = "/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/2016APV/mc/merged_zgmc_llg/merged_pico_llg_ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_zgmc_llg_nfiles_25.root";
    //    filename_pico = "/net/cms37/data1/xf82/BDT/Output/SM_rui_output.root";
  } else if (strcmp( in_name, "GGF")==0) {
    if (strcmp( in_year, "2017")==0)
      filename_pico = "/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/2017/mc/merged_zgmc_llg/merged_pico_llg_GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8_zgmc_llg_nfiles_4.root";
    else if (strcmp( in_year, "2016")==0)
      filename_pico = "/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/2016/mc/merged_zgmc_llg/merged_pico_llg_GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8_zgmc_llg_nfiles_19.root";
    else if (strcmp( in_year, "2018")==0)
      filename_pico = "/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/2018/mc/merged_zgmc_llg/merged_pico_llg_GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8_zgmc_llg_nfiles_2.root";
    else if (strcmp( in_year, "2016APV")==0)
      filename_pico = "/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/2016APV/mc/merged_zgmc_llg/merged_pico_llg_GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8_zgmc_llg_nfiles_20.root";
  } else if (strcmp( in_name, "VBF")==0) {
    if (strcmp( in_year, "2017")==0)
      filename_pico = "/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/2017/mc/merged_zgmc_llg/merged_pico_llg_VBFHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8_zgmc_llg_nfiles_13.root";
    else if (strcmp( in_year, "2016")==0)
      filename_pico = "/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/2016/mc/merged_zgmc_llg/merged_pico_llg_VBFHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8_zgmc_llg_nfiles_4.root";
    else if (strcmp( in_year, "2018")==0)
      filename_pico = "/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/2018/mc/merged_zgmc_llg/merged_pico_llg_VBFHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8_zgmc_llg_nfiles_5.root";
    else if (strcmp( in_year, "2016APV")==0)
      filename_pico = "/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/2016APV/mc/merged_zgmc_llg/merged_pico_llg_VBFHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8_zgmc_llg_nfiles_8.root";

    //    filename_pico = "/net/cms37/data1/xf82/BDT/Output/GGF_rui_output.root";
  } else {
    cout << "ERROR: not an available option" << endl;
    return;
  }

  float w_year = 1.;
    if (strcmp( in_year, "2016")==0)
      w_year = 16.80;//19.5 
    else if(strcmp( in_year, "2016APV")==0)
      w_year = 19.51;
    else if(strcmp( in_year, "2017")==0)
      w_year = 41.48;
    else if (strcmp( in_year, "2018")==0)
      w_year = 59.83;
  //    w_year = 137.52397

  TFile *fpico = new TFile(filename_pico);
  //  TChain *chain = new TChain("tree");
  //  chain->AddFile()
  filename_out += "/";
  filename_out += in_name;
  filename_out += "_";
  filename_out += in_year;
  filename_out += "_output.root";
  TFile *fout = new TFile(filename_out, "RECREATE");//"Output_35_low/DY_output.root", "RECREATE");

    TTree *picotree;
    fpico->GetObject("tree", picotree); //TTree *picotree = (TTree *)fpico->Get("tree");
    TTree *outtree = new TTree("outtree", "Just a little tree");

    //float llphoton_m, ll_m, photon_mva, pt_mass, photon_res, photon_rapidity, l1_rapidity, l2_rapidity, cosTheta,\
    //costheta, phi, min_dR, max_dR, weight;

//Set up a MVA reader    
    TMVA::Tools::Instance();
    TMVA::Reader *reader = new TMVA::Reader( "V:Color:!Silent" ); 
    Float_t var[11];
    reader->AddVariable("costheta", &var[0]);
    reader->AddVariable("cosTheta", &var[1]);
    reader->AddVariable("pt_mass", &var[2]);
    reader->AddVariable("l1_rapidity", &var[3]);
    reader->AddVariable("l2_rapidity", &var[4]);
    reader->AddVariable("photon_rapidity", &var[5]);
    reader->AddVariable("phi", &var[6]);
    reader->AddVariable("photon_mva", &var[7]);
    reader->AddVariable("photon_res", &var[8]);                                                                                                                                                                                                  
    reader->AddVariable("min_dR", &var[9]);                                                                                                                                                                                                      
    reader->AddVariable("max_dR", &var[10]);                     
    //    reader->AddVariable("photon_pt_ratio", &var[11]);

    TString BDTdir = "/net/cms37/data1/rui/Training/";
    BDTdir = BDTdir+input_dir;
    BDTdir = BDTdir+"/weights/TMVAClassification_BDTG.weights.xml";
    reader->BookMVA("BDT method", BDTdir);
    

//Read old branches from input file
    int nll, nphoton, nel, nmu, njet;
    float weight,w_lumi;
    bool el_trig = false;
    bool mu_trig = false;
    //    bool mu_trig2 = false;
    bool stitch_dy = false;

    std::vector<double> *mu_dz = new std::vector<double>;
    std::vector<double> *mu_dxy = new std::vector<double>;
    std::vector<double> *el_dz = new std::vector<double>;
    std::vector<double> *el_dxy = new std::vector<double>;
    std::vector<double> *photon_drmin = new std::vector<double>;

    //std::vector<int> *ll_lepid, *ll_i1, *ll_i2, *el_charge, *mu_charge;
    //std::vector<double> *llphoton_eta, *llphoton_phi, *llphoton_pt, *ll_pt, *ll_eta, *ll_phi, *photon_eta, *photon_phi, *el_pt, *el_eta, *el_phi, *mu_pt, *mu_eta, *mu_phi;
    //std::vector<double> *ll_m, *llphoton_m, *photon_idmva, *photon_pt, *photon_pterr;
    std::vector<double> *dijet_pt = new std::vector<double>;
    std::vector<double> *dijet_eta = new std::vector<double>;
    std::vector<double> *dijet_phi = new std::vector<double>;
    std::vector<double> *dijet_m = new std::vector<double>;
    std::vector<double> *dijet_dr = new std::vector<double>;
    std::vector<double> *dijet_dphi = new std::vector<double>;
    std::vector<double> *dijet_deta = new std::vector<double>;
    std::vector<double> *jet_pt = new std::vector<double>;
    std::vector<double> *jet_eta = new std::vector<double>;
    std::vector<double> *jet_phi = new std::vector<double>;

    std::vector<int> *ll_lepid = new std::vector<int>;
    std::vector<int> *ll_i1 = new std::vector<int>;
    std::vector<int> *ll_i2 = new std::vector<int>;
    std::vector<int> *el_charge = new std::vector<int>;
    std::vector<int> *mu_charge = new std::vector<int>;

    std::vector<double> *llphoton_eta = new std::vector<double>;
    std::vector<double> *llphoton_phi = new std::vector<double>;
    std::vector<double> *llphoton_pt = new std::vector<double>;
    std::vector<double> *llphoton_psi = new std::vector<double>;
    std::vector<int> *llphoton_iph = new std::vector<int>;
    std::vector<int> *llphoton_ill = new std::vector<int>;
    std::vector<double> *llphoton_cosTheta = new std::vector<double>;
    std::vector<double> *llphoton_costheta = new std::vector<double>;
    std::vector<double> *ll_pt = new std::vector<double>;
    std::vector<double> *ll_eta = new std::vector<double>;
    std::vector<double> *ll_phi = new std::vector<double>;
    std::vector<double> *photon_eta = new std::vector<double>;
    std::vector<double> *photon_phi = new std::vector<double>;
    std::vector<double> *el_pt = new std::vector<double>;
    std::vector<double> *el_eta = new std::vector<double>;
    std::vector<double> *el_phi = new std::vector<double>;
    std::vector<double> *mu_pt = new std::vector<double>;
    std::vector<double> *mu_eta = new std::vector<double>;
    std::vector<double> *mu_phi = new std::vector<double>;
    std::vector<double> *ll_m = new std::vector<double>;
    std::vector<double> *llphoton_m = new std::vector<double>;
    std::vector<double> *photon_idmva = new std::vector<double>;
    std::vector<double> *photon_pt = new std::vector<double>;
    std::vector<double> *photon_pterr = new std::vector<double>;

    bool HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL = false;
    bool HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = false;
    bool HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL = false;
    bool HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ = false;
    bool HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 = false;
    bool HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 = false;

//Branches for baseline selection
    picotree->SetBranchAddress("nll", &nll);
    picotree->SetBranchAddress("nphoton", &nphoton);
    picotree->SetBranchAddress("weight", &weight);
    picotree->SetBranchAddress("w_lumi", &w_lumi);
    picotree->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL);
    picotree->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
    picotree->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL);
    picotree->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);
    //    picotree->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24);
    //    picotree->SetBranchAddress("HLT_IsoTkMu24", &HLT_IsoTkMu24);
    //    picotree->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27);
    // picotree->SetBranchAddress("HLT_Mu50", &HLT_Mu50);
    // picotree->SetBranchAddress("HLT_OldMu100", &HLT_OldMu100);
    // picotree->SetBranchAddress("HLT_TkMu100", &HLT_TkMu100);
    picotree->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8);
    picotree->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8);
    picotree->SetBranchAddress("stitch_dy", &stitch_dy);
    picotree->SetBranchAddress("nel", &nel);
    picotree->SetBranchAddress("nmu", &nmu);
    picotree->SetBranchAddress("njet", &njet);
    picotree->SetBranchAddress("el_dz", &el_dz);
    picotree->SetBranchAddress("el_dxy", &el_dxy);
    picotree->SetBranchAddress("mu_dz", &mu_dz);
    picotree->SetBranchAddress("mu_dxy", &mu_dxy);
    picotree->SetBranchAddress("photon_drmin", &photon_drmin);

    picotree->SetBranchAddress("dijet_pt", &dijet_pt);
    picotree->SetBranchAddress("dijet_pt", &dijet_pt);
    picotree->SetBranchAddress("dijet_pt", &dijet_pt);
//Branches for MVA variables
    picotree->SetBranchAddress("ll_lepid", &ll_lepid);
    picotree->SetBranchAddress("ll_m", &ll_m);
    picotree->SetBranchAddress("ll_pt", &ll_pt);
    picotree->SetBranchAddress("ll_eta", &ll_eta);
    picotree->SetBranchAddress("ll_phi", &ll_phi);
    picotree->SetBranchAddress("ll_i1", &ll_i1);
    picotree->SetBranchAddress("ll_i2", &ll_i2);
    picotree->SetBranchAddress("llphoton_iph", &llphoton_iph);
    picotree->SetBranchAddress("llphoton_ill", &llphoton_ill);
    picotree->SetBranchAddress("llphoton_m", &llphoton_m);
    picotree->SetBranchAddress("llphoton_pt", &llphoton_pt);
    picotree->SetBranchAddress("llphoton_eta", &llphoton_eta);
    picotree->SetBranchAddress("llphoton_phi", &llphoton_phi);
    picotree->SetBranchAddress("llphoton_psi", &llphoton_psi);
    picotree->SetBranchAddress("photon_pt", &photon_pt);
    picotree->SetBranchAddress("photon_eta", &photon_eta);
    picotree->SetBranchAddress("photon_phi", &photon_phi);
    picotree->SetBranchAddress("photon_idmva", &photon_idmva);
    picotree->SetBranchAddress("photon_pterr", &photon_pterr);
    picotree->SetBranchAddress("llphoton_cosTheta", &llphoton_cosTheta);
    picotree->SetBranchAddress("llphoton_costheta", &llphoton_costheta);

    picotree->SetBranchAddress("el_charge", &el_charge);
    picotree->SetBranchAddress("el_pt", &el_pt);
    picotree->SetBranchAddress("el_eta", &el_eta);
    picotree->SetBranchAddress("el_phi", &el_phi);
    picotree->SetBranchAddress("mu_charge", &mu_charge);
    picotree->SetBranchAddress("mu_pt", &mu_pt);
    picotree->SetBranchAddress("mu_eta", &mu_eta);
    picotree->SetBranchAddress("mu_phi", &mu_phi);

//Create new branches for output file
    outtree->Branch("photon_pt", "std::vector<double>", &photon_pt);
    outtree->Branch("photon_eta", "std::vector<double>", &photon_eta);
    outtree->Branch("photon_phi", "std::vector<double>", &photon_phi);
    outtree->Branch("photon_idmva", "std::vector<double>", &photon_idmva);
    outtree->Branch("el_pt", "std::vector<double>", &el_pt);
    outtree->Branch("el_eta", "std::vector<double>", &el_pt);
    outtree->Branch("el_phi", "std::vector<double>", &el_phi);
    outtree->Branch("mu_pt", "std::vector<double>", &mu_pt);
    outtree->Branch("mu_eta", "std::vector<double>", &mu_eta);
    outtree->Branch("mu_phi", "std::vector<double>", &mu_phi);
    outtree->Branch("ll_m", "std::vector<double>", &ll_m);
    outtree->Branch("llphoton_m", "std::vector<double>", &llphoton_m);
    outtree->Branch("ll_lepid", "std::vector<int>", &ll_lepid);
    outtree->Branch("photon_drmin", "std::vector<double>", &photon_drmin);
    outtree->Branch("ll_i1", &ll_i1);
    outtree->Branch("ll_i2", &ll_i2);
    outtree->Branch("el_trig", &el_trig);
    outtree->Branch("mu_trig", &mu_trig);
    outtree->Branch("llphoton_iph", &llphoton_iph);
    outtree->Branch("llphoton_ill", &llphoton_ill);
		    
    //    outtree->Branch("el_dz", "std::vector<double>", &el_dz);
    //    outtree->Branch("el_dxy", "std::vector<double>", &el_dxy);
    //    outtree->Branch("mu_dz", "std::vector<double>", &mu_dz);
    //    outtree->Branch("mu_dxy", "std::vector<double>", &mu_dxy);
    float BDT_score;
    outtree->Branch("BDT_score", &BDT_score, "BDT_score/F");
    outtree->Branch("weight", &weight, "weight/F");
    outtree->Branch("w_lumi", &w_lumi, "w_lumi/F");
    outtree->Branch("w_year", &w_year, "w_year/F");
    outtree->Branch("nphoton", &nphoton, "nphoton/I");
    outtree->Branch("nel", &nel, "nel/I");
    outtree->Branch("nmu", &nmu, "nmu/I");
    outtree->Branch("njet", &njet, "njet/I");
    outtree->Branch("nll", &nll, "nll/I");
    outtree->Branch("stitch_dy", &stitch_dy, "stitch_dy/O");
    bool is_test = false;
    outtree->Branch("is_test", &is_test, "is_test/B");



//Run over events
    const Long64_t nevents = picotree->GetEntries();
    for (Long64_t i = 0; i < nevents; i++){
      if ( i % 2 == 0)
	is_test = true;
      else
	is_test = false;
      if (HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ)
	el_trig = true;
      else
	el_trig = false;
      if (HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL|| HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8)
	mu_trig = true;
      else 
	mu_trig = false;
      if (i%100000 == 0) std::cout << "--- ... Processing event: " << i << std::endl;
        picotree->GetEntry(i);

	if (stitch == 1 && !stitch_dy)
	  continue;
	else if (stitch == 2 && stitch_dy)
	  continue;
        if (nphoton>0 && nll>0 && (el_trig || mu_trig) ) {
            TLorentzVector l1 = AssignL1(*ll_lepid, *ll_i1, *ll_i2, *el_charge, *mu_charge, *el_pt, *el_eta, *el_phi, *mu_pt, *mu_eta, *mu_phi);
            TLorentzVector l2 = AssignL2(*ll_lepid, *ll_i1, *ll_i2, *el_charge, *mu_charge, *el_pt, *el_eta, *el_phi, *mu_pt, *mu_eta, *mu_phi);
            TLorentzVector g = AssignGamma(*photon_pt, *photon_eta, *photon_phi, *llphoton_iph);
            TLorentzVector z = AssignZ(*ll_pt, *ll_eta, *ll_phi, *ll_m);
            TLorentzVector h = AssignH(*llphoton_pt, *llphoton_eta, *llphoton_phi, *llphoton_m);
            TLorentzVector q1 = AssignQ1(h);
            TLorentzVector q2 = AssignQ2(h);

            //cout << "ll_pt = " << (*ll_pt)[0] << endl;
            //cout << "ll_i1[0] = " << (*ll_i1)[0] << " ll_i1[1] = " << (*ll_i1)[1] << " ll_i1[2] = " << (*ll_i1)[2] << endl << "ll_i2[0] = " << (*ll_i2)[0] << " ll_i2[1] = " << (*ll_i1)[1] << " ll_i2[2] = " << (*ll_i1)[2] <<endl;
            //if ((*ll_i1)[0] == (*ll_i2)[0]) strange++;

	    //Feed MVA with numerical values of training variables
	    int iph = (*llphoton_iph)[0];
	    var[0] = (Float_t)(*llphoton_costheta)[0];//cos_theta(h, z, l1, l2);
	    var[1] = (Float_t)(*llphoton_cosTheta)[0];//cos_Theta(h, z, q1, q2);
	    var[2] = (Float_t)(*llphoton_pt)[0]/ (*llphoton_m)[0];
	    var[3] = (Float_t)l1.Eta();
	    var[4] = (Float_t)l2.Eta();
	    var[5] = (Float_t)(*photon_eta)[iph];
	    var[6] = (Float_t)(*llphoton_psi)[0];//Getphi(q1, z, l1, l2);
	    var[7] = (Float_t)(*photon_idmva)[iph];
	    var[8] = (Float_t)(*photon_pterr)[iph]/std::cosh((*photon_eta)[iph])/(*photon_pt)[iph];
	    var[9] = (Float_t)GetmindR(l1, l2, g);
	    var[10] =(Float_t)GetmaxdR(l1, l2, g);
	    //	    var[11] = (Float_t)(*photon_pt)[iph]/(*llphoton_m)[0]; 
	    //      var[11] = (Float_t)(*photon_pt)[iph];  
//Compute BDT score
            BDT_score = (float)reader -> EvaluateMVA("BDT method");
            outtree->Fill();
        }
        
    }
    fout->Write();
    fout->Close();
    fpico->Close();

}
