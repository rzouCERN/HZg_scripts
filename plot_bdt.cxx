#include <chrono>
#include <thread>
using namespace std::chrono;
using namespace std;
using namespace ROOT; 
//TH1::SetDefaultSumw2();

//plot bdt score distribution, scores and mllg distribution in each bin

void plot_bdt(bool test=true, string tag="",bool plotBDT=false, bool smooth=false, bool plotscore=false, bool plotmllg=false){
  auto start = high_resolution_clock::now();
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  // string dir = "Output_Input_G_narrow_strict_stitch";
  // string dir = "Output_Input_G_narrow_strict_ratio_stitch";
  string dir = "Output_Input_Run2_stitch";
  //  string dir = "Output_Input_G_strict_narrow_stitch_pt35";
  //  double bestBDT[4] = {-0.46, 0.18, -0.58, 0.06}; // pt35
  double bestBDT[4] = {-0.3, -0.06, 0.1, 0.18}; //run2
  //  double bestBDT[4] = {-0.58, 0.02, 0.38, 0.66}; // ratio
  //double bestBDT[4] = {-0.62, -0.1, 0.3, 0.66}; //strict
  //  string dir = "Output_Run2_fixed_photon_pt_ratio";
  //  string dir = "Output_Run2_fixed_G_narrow_strict";
  //  string dir = "Output_Run2_fixed";
  //  string dir = "Output_Run2_fixed_G_narrow_strict_photon_pt_ratio";
  std::string years[4] = {"2016APV","2016","2017","2018"};
  std::string bkg_samples[2] = {"DY","SMZg"};
  std::string sig_samples[2] {"GGF","VBF"};//,"VBF","VH","ttH"};
  std::vector<std::string> files_bkg, files_sig;
  std::vector<std::string> files_bkg_dy, files_bkg_smzg;
  for (auto year : years) {
    for (auto bkg_sample : bkg_samples){
      files_bkg.push_back(((string)(dir + "/"+bkg_sample+"_"+year+"_output.root")).c_str());
    }
    for (auto sig_sample : sig_samples){
      files_sig.push_back(((string)(dir + "/"+sig_sample+"_"+year+"_output.root")).c_str());
    }
    files_bkg_dy.push_back(((string)(dir + "/"+bkg_samples[0]+"_"+year+"_output.root")).c_str());
    files_bkg_smzg.push_back(((string)(dir + "/"+bkg_samples[1]+"_"+year+"_output.root")).c_str());
  }

  //  std::vector<std::string> files_bkg = {((string)(dir + "/DY_2017_output.root")).c_str(),((string)(dir + "/SMZg_2017_output.root")).c_str()};
  //  std::vector<std::string> files_sig = {((string)(dir + "/GGF_2017_output.root")).c_str()};
  //  ROOT::EnableImplicitMT(); 
  RDataFrame d_bkg("outtree", files_bkg);
  RDataFrame d_sig("outtree", files_sig);
  string is_test = "";
  string name = "";
  if (test) {
    is_test = "is_test==1";
    name = "_test";
  } else {
    is_test = "is_test==0";
    name = "_train";
  }
    
  string el_cuts = "ll_lepid[0] == 11 && nel >= 2 && el_pt[0] > 25 && el_pt[1] > 15";                                                                                                                                                          
  string mu_cuts = "ll_lepid[0] == 13 && nmu >= 2 && mu_pt[0] > 20 && mu_pt[1]> 10";                                                                                                                                                             
  //  string baseline = "nphoton > 0 && ll_m[0] > 50 && llphoton_m[0]+ll_m[0]>185 && njet <=1 && (nel+nmu) ==2 && photon_drmin[llphoton_iph[0]] > 0.4 && photon_pt[llphoton_iph[0]]> 15 && ll_m[0] > 80 && ll_m[0] < 100 && "+is_test+"  && ((abs(photon_eta[llphoton_iph[0]]) < 1.4442 &&  photon_idmva[llphoton_iph[0]] > 0.42) || (abs(photon_eta[llphoton_iph[0]]) > 1.566 && abs(photon_eta[llphoton_iph[0]]) < 2.5 && photon_idmva[llphoton_iph[0]]>0.14))";
  //  string baseline = "nphoton > 0 && ll_m[0] > 50 &&  njet <=1 && (nel+nmu) ==2 && photon_drmin[llphoton_iph[0]] > 0.4 && photon_pt[llphoton_iph[0]]> 15 && ll_m[0] > 80 && ll_m[0] < 100 && "+is_test+"  && ((abs(photon_eta[llphoton_iph[0]]) < 1.4442 &&  photon_idmva[llphoton_iph[0]] > 0.42) || (abs(photon_eta[llphoton_iph[0]]) > 1.566 && abs(photon_eta[llphoton_iph[0]]) < 2.5 && photon_idmva[llphoton_iph[0]]>0.14))";
  string baseline = "nphoton > 0 && ll_m[0] > 50 && llphoton_m[0]+ll_m[0]>185 && photon_drmin[llphoton_iph[0]] > 0.4 && photon_pt[llphoton_iph[0]]> 15 && "+is_test+" && njet <=1 && (nel+nmu) == 2 && ((abs(photon_eta[llphoton_iph[0]]) < 1.4442 &&  photon_idmva[llphoton_iph[0]] > -0.4) || (abs(photon_eta[llphoton_iph[0]]) > 1.566 && abs(photon_eta[llphoton_iph[0]]) < 2.5 && photon_idmva[llphoton_iph[0]]>-0.59))";  
  //  string allcuts =  "((el_trigger = 1 &&"+el_cuts+") || (mu_trigger ==1 && "+mu_cuts+")) && ("+baseline+")";
  string narrowrange = "llphoton_m[0] > 120 && llphoton_m[0] < 130";
  //  string fullrange = "llphoton_m[0] > 100 && llphoton_m[0] < 180";
  string allcuts = "(("+el_cuts+") || ("+mu_cuts+")) && ("+baseline+ "&&"+narrowrange+ ")";
  string allcuts_fullrange = "(("+el_cuts+") || ("+mu_cuts+")) && ("+baseline+")";

  const int nbins = 50;
  //  auto h_bkg = d_bkg.Filter(allcuts).Histo1D(("BDT_bkg", "BDT_bkg", nbins, -1., 1.), "BDT_score", "w_lumi");
  //  auto h_bkg = d_bkg.Filter(allcuts).Define("weightY", "w_lumi*w_year*2").Histo1D({"h_bkg","BDT score",nbins,0.,0.},"BDT_score", "weightY");
  RDataFrame d_bkg_dy("outtree", files_bkg_dy);
  RDataFrame d_bkg_smzg("outtree", files_bkg_smzg);
  if (plotBDT) {
  auto h_bkg = d_bkg.Filter(allcuts).Define("weightY", "w_lumi*w_year*2").Histo1D({"h_bkg","BDT score",nbins,-1,1},"BDT_score", "weightY");
  auto h_bkg_dy = d_bkg_dy.Filter(allcuts).Define("weightY", "w_lumi*w_year*2").Histo1D({"h_bkg_dy","BDT score",nbins,-1,1},"BDT_score", "weightY");
  auto h_bkg_smzg = d_bkg_smzg.Filter(allcuts).Define("weightY", "w_lumi*w_year*2").Histo1D({"h_bkg_smzg","BDT score",nbins,-1,1},"BDT_score", "weightY");
  auto h_sig = d_sig.Filter(allcuts).Define("weightY", "w_lumi*w_year*2").Histo1D({"h_sig","BDT score",nbins,-1,1},"BDT_score", "weightY");
  //  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1000);
  //  c1->cd();
  //  h_bkg->Draw();
  //  h_sig->Draw();
  //  c1->Draw();
  //  c1->SaveAs("h_sig.png");
  //  TCanvas *c3 = new TCanvas("c3", "c3", 1000, 1000);
  //  c3->cd();

  // bestBDT[0] = h_bkg_int->GetXaxis()->GetBinCenter(bestBDT[0]);
  // bestBDT[1] = h_bkg_int->GetXaxis()->GetBinCenter(bestBDT[1]);
  // bestBDT[2] = h_bkg_int->GetXaxis()->GetBinCenter(bestBDT[2]);
  // bestBDT[3] = h_bkg_int->GetXaxis()->GetBinCenter(bestBDT[3]);
  // cout << "{" << bestBDT[0] << ", " << bestBDT[1] << ", " << bestBDT[2] << ", "<< bestBDT[3] << "}" << endl; 
  // cout << TMath::Sqrt(bestSig_2) << endl;

  double numSig = h_sig->Integral();
  double numBkg = h_bkg_dy->Integral()+h_bkg_smzg->Integral();
  THStack * hs_bkg_bdt = new THStack("hsbdt","");
  cout << "numSig = " << numSig << endl << "numBkg = " << numBkg << endl;
  h_sig->Scale(1./numSig);
  h_bkg->Scale(1./numBkg);
  h_bkg_dy->Scale(1./numBkg);
  h_bkg_smzg->Scale(1./numBkg);
  //    h_sig->SetMarkerStyle(8);
  h_sig->SetLineColor(kRed);
  h_sig->SetLineWidth(3);
  h_bkg_smzg->SetLineColor(kBlue);
  h_bkg_smzg->SetLineWidth(3);
  h_bkg_dy->SetLineColor(kOrange+1);
  h_bkg_dy->SetLineWidth(3);
  hs_bkg_bdt->Add(h_bkg_smzg.GetPtr());
  hs_bkg_bdt->Add(h_bkg_dy.GetPtr());
  //  hs_bkg->SetTitle("");
  h_sig->SetTitle("");
  //    h_bkg->SetMarkerStyle(8);
    //    h_sig->SetMarkerColor(kRed);
    //    h_bkg->SetMarkerColor(kBlue);
  cout <<  hs_bkg_bdt->GetMaximum() << " " << h_sig->GetMaximum() << endl;
  double  highY = max(hs_bkg_bdt->GetMaximum(),h_sig->GetMaximum())*1.1;//GetYaxis()->GetXmax();
  cout << highY << endl;
  double  lowY = 0;//h_bkg->GetYaxis()->GetXmin();
  h_sig->SetTitle("Test BDT Categorization");
  
  // cout << "(" <<  bestBDT[0] << ", " << bestBDT[1] << ", " << bestBDT[2] << ", " << bestBDT[3] << ")" << endl;
  TCanvas *c3 = new TCanvas("c3", "c3", 1000, 1000);
  c3->cd();
  //h_sig->Draw();
  double rangeminimum = h_bkg->GetXaxis()->GetXmin();
  double rangemaximum = h_sig->GetXaxis()->GetXmax();
  //  h_bkg->Draw("HIST");
  //  hs_bkg_bdt->Draw();
  //  h_sig->Draw("HIST same");
  //  hs_bkg_bdt->GetYaxis()->SetRangeUser(lowY, highY);
  //  hs_bkg_bdt->GetXaxis()->SetRangeUser(rangeminimum, rangemaximum);
  hs_bkg_bdt->SetMaximum(highY);
  hs_bkg_bdt->SetMinimum(lowY);
  hs_bkg_bdt->Draw();
  h_sig->Draw("HIST same");
  h_sig->GetXaxis()->SetRangeUser(rangeminimum, rangemaximum);

  if (plotscore) {
    TLine *line1 = new TLine(bestBDT[0], lowY, bestBDT[0], highY);
    line1->SetLineColor(kRed);
    TLine *line2 = new TLine(bestBDT[1], lowY, bestBDT[1], highY);
    line2->SetLineColor(kRed);
    TLine *line3 = new TLine(bestBDT[2], lowY, bestBDT[2], highY);
    line3->SetLineColor(kRed);
    TLine *line4 = new TLine(bestBDT[3], lowY, bestBDT[3], highY);
    line4->SetLineColor(kRed);  
    line1->Draw("same");
    line2->Draw("same");
    line3->Draw("same");
    line4->Draw("same");
  }
  if (smooth){
    //    char newname = "h_bkg_smooth";
    char const *newname = "h_bkg_smooth";
    TH1D* h_bkg_smooth = (TH1D*)h_bkg->Clone(newname);
    h_bkg_smooth->Smooth(3);
    h_bkg_smooth->SetLineWidth(3);
    h_bkg_smooth->Draw("HIST L same");
    char const *newname_s = "h_sig_smooth";
    TH1D* h_sig_smooth = (TH1D*)h_sig->Clone(newname_s);
    h_sig_smooth->Smooth(3);
    h_sig_smooth->SetLineWidth(3);
    h_sig_smooth->Draw("HIST L same");
  }
  c3->Draw();
  c3->SaveAs(((string)("BDT_dist_"+dir+name+tag+".png")).c_str());
}
  if (plotmllg) {
    int s_binl,s_binh,b_1_binl,b_1_binh,b_2_binl,b_2_binh;
    double N_s, N_b, Sig=0.;

    //    string bestBDTStrings[4] = {"BDT_score>="+to_string(bestBDT[0])+" && BDT_score<"+to_string(bestBDT[1]), "BDT_score>="+to_string(bestBDT[1])+" && BDT_score<"+to_string(bestBDT[2]), "BDT_score>="+to_string(bestBDT[2])+" && BDT_score<"+to_string(bestBDT[3]), "BDT_score>="+to_string(bestBDT[3])};
    string bestBDTStrings[4] = {"BDT_score>="+to_string(bestBDT[0])+" && BDT_score<"+to_string(bestBDT[1]), "BDT_score>="+to_string(bestBDT[1])+" && BDT_score<"+to_string(bestBDT[2]), "BDT_score>="+to_string(bestBDT[2])+" && BDT_score<"+to_string(bestBDT[3]), "BDT_score>="+to_string(bestBDT[3])};
    //    string bestBDTStrings[4] = {"BDT_score>="+to_string(bestBDT[0])+" && BDT_score<"+to_string(bestBDT[1])+"&& photon_pt[llphoton_iph[0]] <= 35", "BDT_score>="+to_string(bestBDT[1])+" && photon_pt[llphoton_iph[0]] <= 35", "BDT_score>="+to_string(bestBDT[2])+" && BDT_score<"+to_string(bestBDT[3])+"&& photon_pt[llphoton_iph[0]] > 35", "BDT_score>="+to_string(bestBDT[3])+"&& photon_pt[llphoton_iph[0]] > 35"};
    RDataFrame d_bkg_dy("outtree", files_bkg_dy);
    RDataFrame d_bkg_smzg("outtree", files_bkg_smzg);
    int bin = 0;
    for (const string &bdt : bestBDTStrings) {
      string allcuts_fullrange = "(("+el_cuts+") || ("+mu_cuts+")) && ("+baseline+ " && "+bdt+")";
      //      cout << allcuts_fullrange << endl;
      auto h_bkg_dy_mllg = d_bkg_dy.Filter(allcuts_fullrange).Define("weightY", "w_lumi*w_year*2").Define("mllg","llphoton_m[0]").Histo1D({("h_bkg_dy"+to_string(bin)).c_str(),"mllg",150,100,165},"mllg", "weightY");
      auto h_bkg_smzg_mllg = d_bkg_smzg.Filter(allcuts_fullrange).Define("weightY", "w_lumi*w_year*2").Define("mllg","llphoton_m[0]").Histo1D({("h_bkg_smzg"+to_string(bin)).c_str(),"mllg",150,100,165},"mllg", "weightY");
      auto h_sig_mllg = d_sig.Filter(allcuts_fullrange).Define("weightY", "w_lumi*w_year*2").Define("mllg","llphoton_m[0]").Histo1D({("h_sig_mllg"+to_string(bin)).c_str(),"mllg",150,100,165},"mllg", "weightY");

      THStack * hs_bkg = new THStack("hs1",bdt.c_str());
      s_binl=h_sig_mllg->FindBin(120);
      s_binh=h_sig_mllg->FindBin(130);
      N_s = h_sig_mllg->Integral(s_binl,s_binh);
      b_1_binl=h_bkg_dy_mllg->FindBin(120);
      b_1_binh=h_bkg_dy_mllg->FindBin(130);
      b_2_binl=h_bkg_smzg_mllg->FindBin(120);
      b_2_binh=h_bkg_smzg_mllg->FindBin(130);
      N_b = h_bkg_dy_mllg->Integral(b_1_binl,b_1_binh)+h_bkg_smzg_mllg->Integral(b_2_binl,b_2_binh);
      Sig += N_s*N_s/N_b;
      cout << "bin"<< bin << ": N_s: " << N_s << " N_b: " << N_b << " Sig: " << N_s/TMath::Sqrt(N_b) << endl; 
      
      //      THStack * hs = new THStack("hs2"," unstacked");
      TCanvas *c1 = new TCanvas("c1", "c1", 1500, 1000);
      c1->cd();
      h_sig_mllg->Scale(100); //Rui
      // h_bkg_smzg_mllg->GetXaxis()->SetLabelSize(0.03);
      // h_bkg_smzg_mllg->GetYaxis()->SetLabelSize(0.03);
      // h_bkg_smzg_mllg->SetTitle("");
      // h_bkg_smzg_mllg->GetYaxis()->SetTitle("");
      // h_bkg_smzg_mllg->GetYaxis()->SetTitleOffset(0.5);
      // h_bkg_smzg_mllg->GetXaxis()->SetTitle("m_{ll#gamma}");
      h_bkg_smzg_mllg->SetMarkerStyle(8);
      h_sig_mllg->SetMarkerStyle(8);
      h_bkg_smzg_mllg->SetMarkerColor(kBlue);
      h_bkg_smzg_mllg->SetFillColor(kBlue+11);
      double  highY1 = max((h_bkg_smzg_mllg->GetMaximum()+h_bkg_dy_mllg->GetMaximum()),h_sig_mllg->GetMaximum())*1.1;
      hs_bkg->Add(h_bkg_smzg_mllg.GetPtr());
      h_bkg_dy_mllg->SetMarkerStyle(8);
      h_bkg_dy_mllg->SetMarkerColor(kOrange);
      h_bkg_dy_mllg->SetFillColor(kOrange+1);
      hs_bkg->Add(h_bkg_dy_mllg.GetPtr());
      hs_bkg->SetMaximum(highY1);
      //      hs_bkg->SetTitle(bdt.c_str());
      h_sig_mllg->SetMarkerColor(kRed);

      hs_bkg->Draw("");//nostack");
      h_sig_mllg->Draw("hist p same");
      cout << "h_sig: " <<  h_sig_mllg->Integral() << endl;
      TLegend *leg = new TLegend(0.5, 0.6, 0.9, 0.85);
      leg->AddEntry(h_bkg_smzg_mllg.GetPtr(), "SM Z#gamma", "p");
      leg->AddEntry(h_bkg_dy_mllg.GetPtr(), "DY+jets (stacked)", "p");
      leg->AddEntry(h_sig_mllg.GetPtr(), "ggF&VBF h#rightarrow Z#gamma (x100)", "p"); // Rui
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->Draw("same");
      c1->Draw();
      c1->SaveAs(((string)("h_mllg_bin"+to_string(bin)+"_"+dir+name+tag+".png")).c_str());
      c1->Close(); gSystem->ProcessEvents(); delete c1;
      bin += 1;
    }
    cout << "Total Sig: " << TMath::Sqrt(Sig) << endl;
  }
  // if (smooth){  //draw integral of smoothed function
    
  //   TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1000);                                                                                                                                                                                           
  //   c1->cd();
  //   h_sig_smooth->Scale(numSig*200);
  //   h_bkg_smooth->Scale(numBkg);
  //   auto h_bkg_int = h_bkg_smooth->GetCumulative(false);
  //   auto h_sig_int = h_sig_smooth->GetCumulative(false);  
  //   h_bkg_int->Draw();                                                                                                                                                                                                                           
  //   h_sig_int->Draw("same");                                                                                                                                                                                                                    
  //   c1->Draw();                                                                                                                                                                                                                                  
  //   c1->SaveAs(((string)"h_integral_"+dir+name+tag+".png").c_str());
  // }


  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<seconds>(stop - start);

  cout << "processing time: " << duration.count() << "s" << endl;

}

