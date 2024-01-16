#include <string>
#include <chrono>
#include <thread>
using namespace std::chrono;
using namespace std;
using namespace ROOT; 
//TH1::SetDefaultSumw2();

void optbins_dataset(int sample = 0, bool test=true, string tag=""){
  auto start = high_resolution_clock::now();
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  //  string dir = "Output_Input_G_narrow_strict_ratio_stitch";
  string dir = "Output_Input_Run2_stitch";
  //  string dir = "Output_Input_G_narrow_strict_stitch";
  //  string dir = "Output_Input_G_strict_stitch"; 
  //  string dir = "Output_Input_G_strict_narrow_stitch_pt35";
  //  string dir = "Output_Input_G_strict_stitch_pt35";
  //  string dir = "Output_Input_Run2_stitch";
  //  string dir = "Output_Run2_fixed_photon_pt_ratio";
  //  string dir = "Output_Run2_fixed_G_narrow_strict";
  //  string dir = "Output_Run2_fixed";
  //  string dir = "Output_Run2_fixed_G_narrow_strict_photon_pt_ratio";
  std::string years[4] = {"2016APV","2016","2017","2018"};
  //  std::string bkg_samples[2] = {"DY","SMZg"};
  std::string sig_samples[2] {"GGF","VBF"};
  std::vector<std::string> files;
  std::string sample_name;

  for (auto year : years) {
    if (sample == 0) {
      for (auto sig_sample : sig_samples){
	files.push_back(((string)(dir + "/"+sig_sample+"_"+year+"_output.root")).c_str());
	sample_name = "Signal";
      }
    } else if (sample == 1) {
      files.push_back(((string)(dir+"/SMZg_"+year+"_output.root")).c_str());
      sample_name = "SMZg";
    } else if (sample == 2) {
      files.push_back(((string)(dir+"/DY_"+year+"_output.root")).c_str());
      sample_name = "DY";      
    }
  }
  //  std::vector<std::string> files_bkg = {((string)(dir + "/DY_2017_output.root")).c_str(),((string)(dir + "/SMZg_2017_output.root")).c_str()};
  //  std::vector<std::string> files_sig = {((string)(dir + "/GGF_2017_output.root")).c_str()};
  //  ROOT::EnableImplicitMT(); 

  string is_test = "";
  if (test)
    is_test = "is_test==1";
  else
    is_test = "is_test==0";   

  string el_cuts = "ll_lepid[0] == 11 && nel >= 2 && el_pt[0] > 25. && el_pt[1] > 15.";// && (nel==2 || el_pt[3] <7)";                                                                                                                                                          
  string mu_cuts = "ll_lepid[0] == 13 && nmu >= 2 && mu_pt[0] > 20. && mu_pt[1]> 10.";// && (nmu==2 || mu_pt[3] <7)";                                                                                                                                                             
  string baseline = "nphoton > 0 && ll_m[0] > 50 && (nel+nmu) ==2 && llphoton_m[0]+ll_m[0]>185. && photon_drmin[llphoton_iph[0]] > 0.4 && photon_pt[llphoton_iph[0]]> 15. && ll_m[0] > 80. && ll_m[0] < 100. && "+is_test+" && ((abs(photon_eta[llphoton_iph[0]]) < 1.4442 &&  photon_idmva[llphoton_iph[0]] > 0.42) || (abs(photon_eta[llphoton_iph[0]]) > 1.566 && abs(photon_eta[llphoton_iph[0]]) < 2.5 && photon_idmva[llphoton_iph[0]]>0.14))";
  //  string baseline = "nphoton > 0";
  //  string baseline = "nphoton > 0 && ll_m[0] > 50 && llphoton_m[0]+ll_m[0]>185 && photon_drmin[llphoton_iph[0]] > 0.4 && photon_pt[llphoton_iph[0]]> 15 && "+is_test+" && njet <=1 && (nel+nmu) == 2 && ((abs(photon_eta[llphoton_iph[0]]) < 1.4442 &&  photon_idmva[llphoton_iph[0]] > -0.4) || (abs(photon_eta[llphoton_iph[0]]) > 1.566 && abs(photon_eta[llphoton_iph[0]]) < 2.5 && photon_idmva[llphoton_iph[0]]>-0.59))";
  //  string allcuts =  "((el_trigger = 1 &&"+el_cuts+") || (mu_trigger ==1 && "+mu_cuts+")) && ("+baseline+")";
  string narrowrange = "llphoton_m[0] > 120 && llphoton_m[0] < 130";
  //  string narrowrange = "llphoton_m[0] > 120 && llphoton_m[0] < 130"; 
  //  string allcuts = "(("+el_cuts+") || ("+mu_cuts+")) && ("+baseline+ "&&"+narrowrange+ ")";
  //  string allcuts_l = "(("+el_cuts+") || ("+mu_cuts+")) && ("+baseline+ "&&"+narrowrange+ "&& photon_pt[llphoton_iph[0]] <= 35)";
  //  string allcuts_h = "(("+el_cuts+") || ("+mu_cuts+")) && ("+baseline+ "&&"+narrowrange+ "&& photon_pt[llphoton_iph[0]] > 35)";
  string allcuts = "("+baseline+ ")";// && llphoton_m[0] > 100 && llphoton_m[0] < 180)";//&&"+fullrange+ ")";
  //  double bestBDT[4] = {-0.46, 0.18, -0.58, 0.06};// {-0.62, -0.1, 0.3, 0.66};  pt
  double bestBDT[4] = {-0.3, -0.06, 0.1, 0.18}; // Run 2
  //  double bestBDT[4] = {-0.62, -0.1, 0.3, 0.66}; // strict
  //  double bestBDT[4] = {-0.58, 0.02, 0.38, 0.66}; //ratio
  //  string categories[4] = {"photon_pt[llphoton_iph[0]] > 35 && BDT_score>="+to_string(bestBDT[3]),"photon_pt[llphoton_iph[0]] > 35 && BDT_score>="+to_string(bestBDT[2])+" && BDT_score<"+to_string(bestBDT[3]), "photon_pt[llphoton_iph[0]] <= 35 && BDT_score>="+to_string(bestBDT[1]), "photon_pt[llphoton_iph[0]] <= 35 && BDT_score>="+to_string(bestBDT[0])+" && BDT_score<"+to_string(bestBDT[1])};
  string categories[4] = {"BDT_score>="+to_string(bestBDT[3]),"BDT_score>="+to_string(bestBDT[2])+" && BDT_score<"+to_string(bestBDT[3]), "BDT_score>="+to_string(bestBDT[1])+"&& BDT_score<"+to_string(bestBDT[2]), "BDT_score>="+to_string(bestBDT[0])+" && BDT_score<"+to_string(bestBDT[1])}; 
  //  const int nbins = 50;
  //  auto h_bkg = d_bkg.Filter(allcuts).Histo1D(("BDT_bkg", "BDT_bkg", nbins, -1., 1.), "BDT_score", "w_lumi");
  //  auto h_bkg = d_bkg.Filter(allcuts).Define("weightY", "w_lumi*w_year*2").Histo1D({"h_bkg","BDT score",nbins,0.,0.},"BDT_score", "weightY");
  //  auto h_sig = d_sig.Filter(allcuts).Define("weightY", "w_lumi*w_year*2").Histo1D({"h_sig","BDT score",nbins,0.,0.},"BDT_score", "weightY");
  //  */
  //    double bestBDT[4] = {-0.62, -0.1, 0.3, 0.66};
  //  string bestBDTStrings[4] = {"BDT_score>="+to_string(bestBDT[0])+" && BDT_score<"+to_string(bestBDT[1])+"&& photon_pt[llphoton_iph[0]] <= 35", "BDT_score>="+to_string(bestBDT[1])+" && photon_pt[llphoton_iph[0]] <= 35", "BDT_score>="+to_string(bestBDT[2])+" && BDT_score<"+to_string(bestBDT[3])+"&& photon_pt[llphoton_iph[0]] > 35", "BDT_score>="+to_string(bestBDT[3])+"&& photon_pt[llphoton_iph[0]] > 35"};
  //  int bin = 0;
  //  RDataFrame d_bkg_smzg("outtree", files_bkg_smzg);
  //  RDataFrame d_sig_1("outtree", files_sig);
  
  RooRealVar x("x", "m_llg", 100, 180);
  RooRealVar y("y", "photon pT", 15., 80.);
  RooRealVar w("w", "weight", -40., 40.);
  RooRealVar bdt("bdt", "bdt", -1, 1);
  RooRealVar year("year", "year", 2015, 2019);
  RooRealVar lep("lep", "lep", 0, 1); //0 = electron, 1 = muon
  RooRealVar ph_eta("ph_eta", "ph_eta", -3, 3);
  RooRealVar nlep("nlep", "nlep", 0, 10);
  RooRealVar njet("njet", "njet", 0, 10);
  //  string weight = "w";
  RooDataSet rds1("rds1","rds1",RooArgSet(x,y, bdt, w, year, lep, ph_eta, nlep, njet));//,"w");//,weight);
  RooDataSet rds2("rds2","rds2",RooArgSet(x,y, bdt, w, year, lep, ph_eta, nlep, njet));//,"w");//,weight);
  RooDataSet rds3("rds3","rds3",RooArgSet(x,y, bdt, w, year, lep, ph_eta, nlep, njet));//,"w");//,weight);
  RooDataSet rds4("rds4","rds4",RooArgSet(x,y, bdt, w, year, lep, ph_eta, nlep, njet));//,"w");//,weight);

  TChain *treeSig = new TChain("outtree");
  cout << "L93" << endl;
  for (auto file_sig : files) {
    cout << file_sig << endl;
    treeSig->AddFile(file_sig.c_str());
  }
  
  
  // TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1000);
  // c1->cd();

  // TH1F *h_test = new TH1F("h_test","h_test", 50, 100, 180); 
  // //  treeSig->Draw("llphoton_m[0] >> h_test", "2*w_lumi*w_year*(((nel >= 2 && el_pt[0] > 25. && el_pt[1] > 15.) ||(nmu >= 2 && mu_pt[0] > 20. && mu_pt[1]> 10.)) && nphoton > 0 && ll_m[0] > 50 && (nel+nmu) ==2 && llphoton_m[0]+ll_m[0]>185 && photon_drmin[llphoton_iph[0]] > 0.4 && photon_pt[llphoton_iph[0]]> 15 && ll_m[0] > 80 && ll_m[0] < 100 && is_test==1  && ((abs(photon_eta[llphoton_iph[0]]) < 1.4442 &&  photon_idmva[llphoton_iph[0]] > 0.42) || (abs(photon_eta[llphoton_iph[0]]) > 1.566 && abs(photon_eta[llphoton_iph[0]]) < 2.5 && photon_idmva[llphoton_iph[0]]>0.14)))");
  // treeSig->Draw("llphoton_m[0] >> h_test", "2*w_lumi*w_year*(((nmu >= 2 && mu_pt[0] > 20.) ))");//||(nmu >= 2 && mu_pt[0] < 20.) || (nmu >= 2 && mu_pt[0]> 20.)))");
  // //  TH1F *h_testel = new TH1F("el_pt","el_pt", 50, 0, 60);
  // //  treeSig->Draw("el_pt[1] >> el_pt", "2*w_lumi*w_year*(@el_pt.size() >= 2)");
  // c1->Draw();
  // c1->SaveAs("test1_1.png");

  // TCanvas *c4 = new TCanvas("c4", "c4", 1000, 1000);
  // c4->cd();

  // TH1F *h_testmu = new TH1F("mu_pt","mu_pt", 50, 0, 60);
  // treeSig->Draw("mu_pt[1] >> mu_pt", "2*w_lumi*w_year*(nmu >= 2)");
  // c4->Draw();
  // c4->SaveAs("test1_mupt_sub.png");

  
  TFile *file = new TFile("temp.root","RECREATE");
  file->cd();
  TTree *temp_tree;
  temp_tree = treeSig->CopyTree(allcuts.c_str());

  // TCanvas *c2 = new TCanvas("c2", "c2", 1000, 1000);
  // c2->cd();

  // TH1F *h_test1 = new TH1F("h_test1","h_test1", 50, 100, 180);
  // temp_tree->Draw("llphoton_m[0] >> h_test1","2*w_lumi*w_year");
  // c2->Draw();
  // c2->SaveAs("test2.png");


  std::vector<float> *llphoton_m = new std::vector<float>(10);
  std::vector<float> *photon_pt = new std::vector<float>(10);
  std::vector<double> *el_pt = new std::vector<double>;
  std::vector<double> *mu_pt = new std::vector<double>;
  std::vector<float> *photon_eta = new std::vector<float>(10);
  std::vector<int> *ll_lepid = new std::vector<int>(10);
  std::vector<int> *llphoton_iph = new std::vector<int>(10);
  
  float w_lumi, w_year, BDT_score;
  int njet_,nel,nmu;
  
  temp_tree->SetBranchAddress("llphoton_m", &llphoton_m);
  temp_tree->SetBranchAddress("photon_pt", &photon_pt);
  temp_tree->SetBranchAddress("llphoton_iph", &llphoton_iph);
  temp_tree->SetBranchAddress("w_lumi", &w_lumi);
  temp_tree->SetBranchAddress("w_year", &w_year);
  temp_tree->SetBranchAddress("photon_eta", &photon_eta);
  temp_tree->SetBranchAddress("BDT_score", &BDT_score);
  temp_tree->SetBranchAddress("nel", &nel);
  temp_tree->SetBranchAddress("nmu", &nmu);
  temp_tree->SetBranchAddress("njet", &njet_);
  temp_tree->SetBranchAddress("ll_lepid", &ll_lepid);
  temp_tree->SetBranchAddress("el_pt", &el_pt);
  temp_tree->SetBranchAddress("mu_pt", &mu_pt);

  
  //  RooDataSet* sigSample = RooDataSet::read("../fitting/data/FullSig_deathvalley_v3_untagged.dat", RooArgList(x, y, bdt, w, year, lep, ph_eta, nlep, njet));
  //  RooDataSet* bkgSample = RooDataSet::read("../fitting/data/SMZg_deathvalley_v3_untagged.dat, ../fitting/data/DY_deathvalley_v3_untagged.dat", RooArgList(x, y, bdt, w, year, lep, ph_eta, nlep, njet));
  //  d_bkg_smzg.Define("x", "llphoton_m[0]").Define("y", "photon_pt[llphoton_iph[0]]");
  //  RooArgList arglist_vars = RooArgList(x);
  //  RooDataSet ds("SMZg_2017_untagged.dat","SMZg_untagged", RooArgList(x));//, bdt, w, year, lep, ph_eta, nlep, njet));
  // auto save_in_ds = [&arglist_vars,&ds] (const std::vector<double>& rdf_vars) -> void {
  //   auto ctr = 0u;
  //   auto vit = arglist_vars.createIterator();
  //   while(auto rrv = static_cast<RooRealVar*>(vit->Next()))
  //     rrv->setVal(rdf_vars[ctr++]);
  //   ds.add(arglist_vars);
  // };
  
  // d_bkg_smzg.Define("rdf_vars_list", "llphoton_m[0]");//,photon_pt[llphoton_iph[0]]");
  // d_bkg_smzg.Foreach(save_in_ds,{"rdf_vars_list"});//, "BDT_score", "2017", ""});
  // TChain *treeSMZg = new TChain("outtree");
  // for (auto file_bkg : files_bkg_smzg) {
  //   treeSMZg->AddFile(file_bkg.c_str());
  // }
  // TTree *temp_tree;
  // temp_tree = treeSMZg->CopyTree(allcuts_l.c_str());
  //  double bestBDT[4] = {-0.62, -0.1, 0.3, 0.66};

  RooPlot* xframe = x.frame() ;
  xframe->Draw();
  for(int i=0; i<temp_tree->GetEntries(); i++) {
    temp_tree->GetEntry(i);
    x = (*llphoton_m)[0];
    y = (*photon_pt)[(*llphoton_iph)[0]];
    bdt = BDT_score;
    w = w_lumi*w_year*2;
    //    year=stoi(str_year);
    if (w_year < 18)// 16.80)
      year = 2016;
    else if (w_year < 30) // == 19.51)
      year = 2016;
    else if (w_year < 45)//== 41.48)
      year = 2017;
    else if (w_year <60)// 59.83)
      year = 2018;
    else 
      cout << "NO YEAR INFO: " << w_year << endl;
    if ((*ll_lepid)[0] == 11) { 
      if (nel >= 2 and el_pt->at(0) > 25. and el_pt->at(1) > 15.)
	lep = 0;
      else
	continue;
    } else {
      if (nmu >= 2 and mu_pt->at(0) > 20. and mu_pt->at(1) >10.)
	lep = 1;
      else
	continue;
    }     
    //    lep = (*ll_lepid)[0] == 11 ? 0: 1;
    nlep = nel+nmu;
    njet = njet_;
    ph_eta = (*photon_eta)[(*llphoton_iph)[0]];
    if (BDT_score > bestBDT[3])
      rds1.add(RooArgSet(x,y, bdt, w, year, lep, ph_eta, nlep, njet),w_lumi*w_year*2);
    else if (BDT_score > bestBDT[2])
      rds2.add(RooArgSet(x,y, bdt, w, year, lep, ph_eta, nlep, njet),w_lumi*w_year*2);
    else if (BDT_score > bestBDT[1])
      rds3.add(RooArgSet(x,y, bdt, w, year, lep, ph_eta, nlep, njet),w_lumi*w_year*2);
    else if (BDT_score > bestBDT[0])
      rds4.add(RooArgSet(x,y, bdt, w, year, lep, ph_eta, nlep, njet),w_lumi*w_year*2);
    //    cout << i <<< endl;
  }

  TCanvas *c3 = new TCanvas("c3", "c3", 1000, 1000);                                                                                                                                                                                             
  c3->cd(); 
  rds1.plotOn(xframe) ;
  xframe->Draw();
  
  cout << "Rui: " << rds1.sumEntries() << endl;
  c3->Draw();
  c3->SaveAs(((string)("test"+to_string(sample)+".png")).c_str());
  
  file->Close();
  rds1.Print();
  rds2.Print();
  rds3.Print();
  rds4.Print();
  rds1.write(("data/"+sample_name+"_untagged1"+tag+".dat").c_str());
  rds2.write(("data/"+sample_name+"_untagged2"+tag+".dat").c_str());
  rds3.write(("data/"+sample_name+"_untagged3"+tag+".dat").c_str());
  rds4.write(("data/"+sample_name+"_untagged4"+tag+".dat").c_str());
  //  file->Close();
  delete file;
  //    os.system("rm -rf temp.root")

  // RooDataSet rds("rds","rds",treeSMZg, RooArgSet(x));//,treeSMZg);//,"x>100") ;
  
  //  rds.Print();
      //      RooDataSet u_sigSample("f_sigSample", "f_sigSample", sigSample, RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "nlep <= 2 && njet < 2","w");
      //      RooDataSet u_bkgSample("f_bkgSample", "f_bkgSample", bkgSample, RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "nlep <= 2 && njet < 2","w");

    /*
  // cout << "Run2 signi = " <<GetSignificance(run2BDT, *BDT_sig, *BDT_bkg) << endl << "best signi = " << GetSignificance(bestBDT, *BDT_sig, *BDT_bkg) << endl;
  // //cout << "Best significance = " << bestSig << endl;
    double numSig = h_sig->Integral();
    double numBkg = h_bkg->Integral();
    cout << "numSig = " << numSig << endl << "numBkg = " << numBkg << endl;
    h_sig->Scale(1./numSig);
    h_bkg->Scale(1./numBkg);
    //    h_sig->SetMarkerStyle(8);
    h_sig->SetLineColor(kRed);
    h_sig->SetLineWidth(3);
    h_bkg->SetLineColor(kBlue);
    h_bkg->SetLineWidth(3);
    h_bkg->SetTitle("");
    h_sig->SetTitle("");
    //    h_bkg->SetMarkerStyle(8);
    //    h_sig->SetMarkerColor(kRed);
    //    h_bkg->SetMarkerColor(kBlue);
    double  highY = max(h_bkg->GetMaximum(),h_sig->GetMaximum())*1.1;//GetYaxis()->GetXmax();
    double  lowY = 0;//h_bkg->GetYaxis()->GetXmin();
    h_sig->SetTitle("Test BDT Categorization");
    TLine *line1 = new TLine(bestBDT[0], lowY, bestBDT[0], highY);
    line1->SetLineColor(kRed);
    TLine *line2 = new TLine(bestBDT[1], lowY, bestBDT[1], highY);
    line2->SetLineColor(kRed);
    TLine *line3 = new TLine(bestBDT[2], lowY, bestBDT[2], highY);
    line3->SetLineColor(kRed);
    TLine *line4 = new TLine(bestBDT[3], lowY, bestBDT[3], highY);
    line4->SetLineColor(kRed);
  
  // TLine *line11 = new TLine(run2BDT[0], lowY, run2BDT[0], 0.015);
  // line11->SetLineColor(kGreen);
  // TLine *line21 = new TLine(run2BDT[1], lowY, run2BDT[1], 0.015);
  // line21->SetLineColor(kGreen);
  // TLine *line31 = new TLine(run2BDT[2], lowY, run2BDT[2], 0.015);
  // line31->SetLineColor(kGreen);
  // TLine *line41 = new TLine(run2BDT[3], lowY, run2BDT[3], 0.015);
  // line41->SetLineColor(kGreen);

  // cout << "(" <<  bestBDT[0] << ", " << bestBDT[1] << ", " << bestBDT[2] << ", " << bestBDT[3] << ")" << endl;
    TCanvas *c3 = new TCanvas("c3", "c3", 1000, 1000);
    c3->cd();
    //h_sig->Draw();
    double rangeminimum = h_bkg->GetXaxis()->GetXmin();
    double rangemaximum = h_sig->GetXaxis()->GetXmax();
    h_bkg->Draw("HIST");
    h_sig->Draw("HIST same");
    h_bkg->GetYaxis()->SetRangeUser(lowY, highY);
    h_bkg->GetXaxis()->SetRangeUser(rangeminimum, rangemaximum);
    h_sig->GetXaxis()->SetRangeUser(rangeminimum, rangemaximum);
    line1->Draw("same");
    line2->Draw("same");
    line3->Draw("same");
    line4->Draw("same");
    c3->Draw();
    c3->SaveAs(((string)("BDT_"+dir+name+tag+".png")).c_str());
    */
  
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<seconds>(stop - start);

  cout << "processing time: " << duration.count() << "s" << endl;

}

