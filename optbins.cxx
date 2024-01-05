#include <chrono>
#include <thread>
using namespace std::chrono;
using namespace std;
using namespace ROOT; 
//TH1::SetDefaultSumw2();

void optbins(bool plot, bool test=true, string tag=""){
  auto start = high_resolution_clock::now();
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  string dir = "Output_Input_G_narrow_strict_ratio_stitch";
  std::string years[4] = {"2016APV","2016","2017","2018"};
  std::string bkg_samples[2] = {"DY","SMZg"};
  std::string sig_samples[2] {"GGF","VBF"};
  std::vector<std::string> files_bkg, files_sig;
  for (auto year : years) {
    for (auto bkg_sample : bkg_samples){
      files_bkg.push_back(((string)(dir + "/"+bkg_sample+"_"+year+"_output.root")).c_str());
    }
    for (auto sig_sample : sig_samples){
      files_sig.push_back(((string)(dir + "/"+sig_sample+"_"+year+"_output.root")).c_str());
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

  string el_cuts = "ll_lepid[0] == 11 && nel >= 2 && el_pt[0] > 25 && el_pt[1] > 15";// && (nel==2 || el_pt[3] <7)";                                                                                                                                                          
  string mu_cuts = "ll_lepid[0] == 13 && nmu >= 2 && mu_pt[0] > 20 && mu_pt[1]> 10";// && (nmu==2 || mu_pt[3] <7)";                                                                                                                               //'New' tight baseline                              
  string baseline = "nphoton > 0 && ll_m[0] > 50 && njet <=1 && (nel+nmu) ==2 && llphoton_m[0]+ll_m[0]>185 && photon_drmin[llphoton_iph[0]] > 0.4 && photon_pt[llphoton_iph[0]]> 15 && ll_m[0] > 80 && ll_m[0] < 100 && "+is_test+"  && ((abs(photon_eta[llphoton_iph[0]]) < 1.4442 &&  photon_idmva[llphoton_iph[0]] > 0.42) || (abs(photon_eta[llphoton_iph[0]]) > 1.566 && abs(photon_eta[llphoton_iph[0]]) < 2.5 && photon_idmva[llphoton_iph[0]]>0.14))";
  //Run 2 loose baseline
  //  string baseline = "nphoton > 0 && ll_m[0] > 50 && llphoton_m[0]+ll_m[0]>185 && photon_drmin[llphoton_iph[0]] > 0.4 && photon_pt[llphoton_iph[0]]> 15 && "+is_test+" && njet <=1 && (nel+nmu) == 2 && ((abs(photon_eta[llphoton_iph[0]]) < 1.4442 &&  photon_idmva[llphoton_iph[0]] > -0.4) || (abs(photon_eta[llphoton_iph[0]]) > 1.566 && abs(photon_eta[llphoton_iph[0]]) < 2.5 && photon_idmva[llphoton_iph[0]]>-0.59))";
  //  string allcuts =  "((el_trigger = 1 &&"+el_cuts+") || (mu_trigger ==1 && "+mu_cuts+")) && ("+baseline+")";
  string narrowrange = "llphoton_m[0] > 120 && llphoton_m[0] < 130";
  //  string narrowrange = "llphoton_m[0] > 120 && llphoton_m[0] < 130"; 
  string allcuts = "(("+el_cuts+") || ("+mu_cuts+")) && ("+baseline+ "&&"+narrowrange+ ")";
  string allcuts_fullrange = "(("+el_cuts+") || ("+mu_cuts+")) && ("+baseline+ ")";//&&"+fullrange+ ")";
  
  RDataFrame d_bkg("outtree", files_bkg);
  RDataFrame d_sig("outtree", files_sig);           
  const int nbins = 50;
  auto h_bkg = d_bkg.Filter(allcuts).Define("weightY", "w_lumi*w_year*2").Histo1D({"h_bkg","BDT score",nbins,-1,1},"BDT_score", "weightY");
  auto h_sig = d_sig.Filter(allcuts).Define("weightY", "w_lumi*w_year*2").Histo1D({"h_sig","BDT score",nbins,-1,1},"BDT_score", "weightY"); 
  h_bkg->Smooth(2);
  h_sig->Smooth(2);
  //  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1000);
  //  c1->cd();
  //  h_bkg->Draw();
  //  h_sig->Draw();
  //  c1->Draw();
  //  c1->SaveAs("h_sig.png");
  //  TCanvas *c3 = new TCanvas("c3", "c3", 1000, 1000);
  //  c3->cd();
  auto h_bkg_int = h_bkg->GetCumulative();
  //  cout << h_bkg_int->GetBinContent(1) << " " << h_bkg_int->GetBinContent(nbins) << endl;
  //  cout << h_bkg->GetNbinsX() << " " << h_bkg_int->GetNbinsX() << endl;
  auto h_sig_int = h_sig->GetCumulative();
  //  cout << h_sig_int->GetBinContent(1) << " " <<h_sig_int->GetBinContent(nbins) << endl;
  //  h_bkg_int->Smooth();
  //  h_sig_int->Smooth();
  //  h_bkg_int->Draw();
  //  h_sig_int->Draw("same");
  //  c3->Draw();
  //  c3->SaveAs("h_sig_int.png");

  //  int min_BDT = h_bkg_int->GetXaxis()->GetXmin()*100;
  //  int max_BDT = h_sig_int->GEtXaxis()->GetXmax()*100;
  
  // float whole = floor(n);
  // float decimal = n - whole;
  // std::cout << whole << "\n";
  // std::cout << decimal << "\n";

  // std::cin.get();
  //  delete d_bkg, d_sig, h_bkg, h_sig;
  
  double NBkg3, NBkg2, NBkg1, NBkg0;
  double NSig3, NSig2, NSig1, NSig0;
  double Sig3_2, Sig2_2, Sig1_2, Sig0_2, SigTot_2;
  double bestSig_2 = 0.;
  double bestSigs_2[4]={0.,0.,0.,0.};
  double NSigs[4];
  double NBkgs[4];
  double bestBDT[4] = {};
  double arr[nbins][2] = {{}};
  //  double* arr_sig=new double [nbins][2];
  for (int bin(1); bin < nbins+1; bin++){
    arr[bin-1][0] = h_bkg_int->GetBinContent(bin);
    arr[bin-1][1] = h_sig_int->GetBinContent(bin);
    //    arr_bkg[bin] = h_bkg_int->GetBinContent(bin);
    //    arr_sig[bin] = h_sig_int->GetBinContent(bin);
  }
  cout << "nsig_integral:" << h_sig->Integral() << " h_sig_int:" << h_sig_int->GetBinContent(nbins) << endl;
  cout << "nbkg_integral:" << h_bkg->Integral() << " h_bkg_int:" << h_bkg_int->GetBinContent(nbins) << endl;
  //  std::vector<double> v_bkg(std::begin(arr),std::end(arr));
  //  delete arr;
    //  double* arr_sig=new double [nbins];
    //  arr_sig=h_sig_int->GetArray();
  //  std::vector<double> v_sig(std::begin(arr),std::end(arr));
  //  delete arr;

  //  double run2BDT[4] = {-0.038, 0.0233, 0.0628, 0.0766};
  //  volatile int i = 3;
  //  std::mutex i_mutex; 
  //  const std::lock_guard<std::mutex> lock(i_mutex);
  double nsig_tot = h_sig_int->GetBinContent(nbins);
  double nbkg_tot = h_bkg_int->GetBinContent(nbins);
  for (int i(3); i < nbins; i++){
    //  while (i < nbins) {
    auto pair = arr[i];
    //    cout << "i: "<< i << " nsig_tot: " << nsig_tot << " arr[i][1]: "   <<pair[1] << endl; 
    //    cout << "i: "<< i << " nbkg_tot: " << nbkg_tot << " arr[i][1]: "   <<pair[0] << endl;
    NBkg3 = nbkg_tot-pair[0];//arr_bkg[i];
    NSig3 = nsig_tot-pair[1];//arr_sig[i];
    //    cout << "3: " <<NBkg3 << " " << NSig3 << endl;
    if (NSig3 < 2) continue;
    if (NBkg3 == 0) Sig3_2 = 0.;
    else Sig3_2 =  pow(NSig3,2)/(NBkg3);//+NSig3);
    for (int j(2); j < i; j++){
      pair = arr[j];
      NBkg2 = arr[i][0]-pair[0];//arr_bkg[j]-NBkg3;
      NSig2 = arr[i][1]-pair[1]; //arr_sig[j]-NSig3;
      //      cout << "2: "<< NBkg2 << " " << NSig2 << " " << pair[1] << endl;
      if (NSig2 < 2) continue;
      if (NBkg2 == 0) Sig2_2 = 0.;
      else Sig2_2 = pow(NSig2,2)/(NBkg2);//+NSig2);
      for (int k(1); k < j; k++){
	//      while (k < j) {
	pair = arr[k];
	NBkg1 = arr[j][0]-pair[0];
	NSig1 = arr[j][1]-pair[1];
	//	cout << "1: "<< NBkg1 << " " << NSig1 << " " << pair[1] << endl;
	if (NSig1 < 2) continue;
  	if (NBkg1 == 0) Sig1_2 = 0.; 
  	else Sig1_2 = pow(NSig1,2)/(NBkg1);//+NSig1);
	//	volatile int m = 0;
	//	std::mutex m_mutex;
	//	const std::lock_guard<std::mutex> lock(m_mutex);
	for (int m(0); m < k; m++){
	  pair = arr[m];
	  NBkg0 = arr[k][0]-pair[0];
	  NSig0 = arr[k][1]-pair[1];
	  //	  cout << "0: "<< NBkg0 << " " << NSig0 << " " << pair[1] << endl;
	  //	  cout << "full arr: " << arr[m][1] << " " << arr[k][1] << " " << arr[j][1] << " " << arr[i][1] << endl;
	  if (NSig0 < 2) continue;
  	  if (NBkg0 == 0) Sig0_2 = 0.;
  	  else Sig0_2 = pow(NSig0,2)/(NBkg0);//+NSig0);
	  SigTot_2 = Sig3_2+Sig2_2+Sig1_2+Sig0_2;
	  //  	  SigTot = TMath::Sqrt(Sig3*Sig3 + Sig2*Sig2 + Sig1*Sig1 + Sig0*Sig0);
	  //	  cout << "(" << i << ", " << j << ", "<< k << ", " << m <<")" << endl;
  	  if (SigTot_2 > bestSig_2){ 
  	    bestSig_2 = SigTot_2;
	    bestSigs_2[0] = Sig0_2;
	    bestSigs_2[1] = Sig1_2;
	    bestSigs_2[2] = Sig2_2;
	    bestSigs_2[3] = Sig3_2;
	    NSigs[0]=NSig0;
            NSigs[1]=NSig1;
            NSigs[2]=NSig2;
            NSigs[3]=NSig3;
  	    bestBDT[0] = m;
  	    bestBDT[1] = k;
  	    bestBDT[2] = j;
  	    bestBDT[3] = i;
  	  }
	  //	  m = m+1;
  	}
	//	k = k+1;
      }
      //      j = j+1;
    }
    //    i = i+1;
  }

  //  delete arr;//_bkg, arr_sig;
  cout << "bdt bin: ";
  for (auto best: bestBDT){
    cout << best << ",";
  }
  cout << endl;
  bestBDT[0] = h_bkg_int->GetXaxis()->GetBinCenter(bestBDT[0]);
  bestBDT[1] = h_bkg_int->GetXaxis()->GetBinCenter(bestBDT[1]);
  bestBDT[2] = h_bkg_int->GetXaxis()->GetBinCenter(bestBDT[2]);
  bestBDT[3] = h_bkg_int->GetXaxis()->GetBinCenter(bestBDT[3]);
  cout << "{" << bestBDT[0] << ", " << bestBDT[1] << ", " << bestBDT[2] << ", "<< bestBDT[3] << "}" << endl; 
  for (auto NSig: NSigs){
    cout << NSig << ",";
  }
  cout << endl;
  cout << TMath::Sqrt(bestSigs_2[0]) << ","<< TMath::Sqrt(bestSigs_2[1]) << ","<< TMath::Sqrt(bestSigs_2[2]) << ","<< TMath::Sqrt(bestSigs_2[3]) << endl;
  cout << TMath::Sqrt(bestSig_2) << endl;
  // BDT_sig->Rebin(5);
  // BDT_bkg->Rebin(5);
  // BDT_sig->GetXaxis()->SetRangeUser(-0.1, 0.2);
  // BDT_bkg->GetXaxis()->SetRangeUser(-0.1, 0.2);
  // BDT_bkg->GetYaxis()->SetRangeUser(0., 0.028);
  //  */
  if (plot) {
    string name = "";
    if (test)
      name = "_test";
    else
      name = "_train";
    //    double bestBDT[4] = {-0.62, -0.1, 0.3, 0.66};
    string bestBDTStrings[4] = {"BDT_score>="+to_string(bestBDT[0])+" && BDT_score<"+to_string(bestBDT[1]), "BDT_score>="+to_string(bestBDT[1])+" && BDT_score<"+to_string(bestBDT[2]), "BDT_score>="+to_string(bestBDT[2])+" && BDT_score<"+to_string(bestBDT[3]), "BDT_score>="+to_string(bestBDT[3])};
    int bin = 0;
    RDataFrame d_bkg_1("outtree", files_bkg);
    RDataFrame d_sig_1("outtree", files_sig);
    for (const string &bdt : bestBDTStrings) {

      string allcuts_fullrange = "(("+el_cuts+") || ("+mu_cuts+")) && ("+baseline+" && "+bdt+")";
      cout << allcuts_fullrange << endl;
      auto h_bkg_mllg = d_bkg_1.Filter(allcuts_fullrange).Define("weightY", "w_lumi*w_year*2").Define("mllg","llphoton_m[0]").Histo1D({("h_bkg_dy"+to_string(bin)).c_str(),"mllg",50,95,185},"mllg", "weightY");
      auto h_sig_mllg = d_sig_1.Filter(allcuts_fullrange).Define("weightY", "w_lumi*w_year*2").Define("mllg","llphoton_m[0]").Histo1D({("h_sig"+to_string(bin)).c_str(),"mllg",80,0.,0.},"mllg", "weightY");

      TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1000);  
      c1->cd();
      h_sig_mllg->Scale(100);
      h_bkg_mllg->GetXaxis()->SetLabelSize(0.03);
      h_bkg_mllg->GetYaxis()->SetLabelSize(0.03);
      h_bkg_mllg->SetTitle("");
      h_bkg_mllg->GetYaxis()->SetTitle("");//Number of Events");
      h_bkg_mllg->GetYaxis()->SetTitleOffset(0.5);
      h_bkg_mllg->GetXaxis()->SetTitle("m_{ll#gamma}");
      h_bkg_mllg->SetMarkerStyle(8);
      h_sig_mllg->SetMarkerStyle(8);
      h_bkg_mllg->SetMarkerColor(kBlack);
      h_sig_mllg->SetMarkerColor(kRed);
      h_bkg_mllg->SetTitle(bdt.c_str());
      double  highY1 = max(h_bkg_mllg->GetMaximum(),h_sig_mllg->GetMaximum())*1.1;
      h_bkg_mllg->GetYaxis()->SetRangeUser(0, highY1);
      h_bkg_mllg->Draw("");
      h_sig_mllg->Draw("same");
      cout << "h_sig: " <<  h_sig_mllg->Integral() << endl;
      TLegend *leg = new TLegend(0.6, 0.6, 0.9, 0.85);
      leg->AddEntry(h_bkg_mllg.GetPtr(), "SM Z#gamma & DY", "p");
      leg->AddEntry(h_sig_mllg.GetPtr(), "ggF&VBF h#rightarrow Z#gamma (x100)", "p");
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->Draw("same");
      c1->Draw();
      c1->SaveAs(((string)("h_mllg_bin"+to_string(bin)+"_"+dir+name+tag+".png")).c_str());
      c1->Close(); gSystem->ProcessEvents(); delete c1;
      bin += 1;
    }
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
  }
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<seconds>(stop - start);

  cout << "processing time: " << duration.count() << "s" << endl;

}

