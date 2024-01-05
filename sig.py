#!/usr/bin/env python

import ROOT
import math

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.TH1.SetDefaultSumw2(True)

Ch = ROOT.TChain("outtree")

Ch.AddFile("Output_30_high/SMZg_output.root")
Ch.AddFile("Output_30_high/DY_output.root")
#Ch.AddFile("Output_photonpT/SMZg_output.root")   
#Ch.AddFile("Output_nocosT/SMZg_output.root")

Ch1 = ROOT.TChain("outtree")
Ch1.AddFile("Output_30_high/GGF_output.root")
#Ch1.AddFile("Output_30_low_nopt/SMZg_output.root")
#Ch.AddFile("Output_30_low/SMZg_output.root") 
#Ch.AddFile("Output/DY_output.root")
#Ch.AddFile("Output/SMZg_output.root") 

#Ch = ROOT.TChain("TreeB")
#Ch.AddFile("../output/BDTnew_mllg.root")

ROOT.TH1.SetDefaultSumw2()

can = ROOT.TCanvas("can", "can")

def PlotHist(log, title, xlabel, ylabel, YlabelRatio, nbins, rangemin, rangemax, filename, thelist):
	leg = ROOT.TLegend(0.6, 0.6, 0.9, 0.85)
#	leg = ROOT.TLegend(0.6, 0.3, 0.9, 0.55)
	if len(thelist) > 2:
		legr = ROOT.TLegend(0.6, 0.55, 0.77, 0.85)
	histchain = []
	legend = []
	c = 0
	color = [ROOT.kRed, ROOT.kRed+2, ROOT.kBlue, ROOT.kBlack, ROOT.kGreen+2, ROOT.kViolet, ROOT.kGreen, ROOT.kPink-5, ROOT.kOrange, ROOT.kTeal, ROOT.kOrange]
#	color = [ROOT.kRed, ROOT.kRed+2, ROOT.kOrange, ROOT.kOrange+2, ROOT.kViolet, ROOT.kViolet+2, ROOT.kGreen, ROOT.kGreen+2, ROOT.kPink, ROOT.kTeal, ROOT.kOrange]
        rangeminimum = rangemin
        rangemaximum = rangemax
	a = 0

	for alist in thelist:
		l = len(alist)

		Chains = alist[l-5]
		Chain070 = Chains[0]
#		Chain70140 = Chains[1]
		data       = alist[l-4]
		cut        = alist[l-3]
		name       = alist[l-2]
		legname    = alist[l-1]
		print cut
		
		legend.append(legname)
#		hist  = ROOT.TH2F(name,";"+xlabel+";;"+ylabel,nbins, 0, 50, nbins,rangemin,rangemax)
		hist = ROOT.TH1F(name,";"+xlabel+";;"+ylabel,nbins,rangemin,rangemax)
		Chain070.Draw(data +" >>"+name, cut)
		print data, " ", cut
#		for Chain in Chains[1:]:
#			Chain.Draw(data +" >>+ "+name, cut)
		histchain.append(hist)


	c = 0
#	histchain[0].SetLineColor(color[c])
#	histchain[0].Draw("AP")
	histchain[0].GetXaxis().SetLabelSize(0.06)
        histchain[0].GetXaxis().SetTitleSize(0.06)
        histchain[0].GetXaxis().SetTitleOffset(0.7)
        histchain[0].GetXaxis().SetRangeUser(rangeminimum, rangemaximum)
        histchain[0].GetYaxis().SetLabelSize(0.05)
        histchain[0].SetTitle(title)
	histchain[0].SetMarkerColor(color[c])
	histchain[0].SetMarkerStyle(21)#SetFillStyle(3004)
#	histchain[0].SetFillColor(color[c])
#	histchain[0].GetYaxis().SetRangeUser(0.0002, 0.1)
	histchain[0].Draw()#AP")
	for histo in histchain:
#		histo.SetFillStyle(1001)
		histo.SetMarkerStyle(21)
		histo.SetMarkerColor(color[c])
		histo.SetLineColor(color[c])
#		histo.SetFillColor(color[c])
                histo.SetLineWidth(3)
		histo.Draw("sameP")
		print legend[c]
		leg.AddEntry(histo,legend[c],"l")
		c += 1

	
	leg.SetFillStyle(0)
	leg.SetBorderSize(0)
	leg.Draw("same")
# 	pad1.SetLogy(log)

# #	legr.AddEntry(ratio1,"W #rightarrow #mu#nu","l")
# 	if len(thelist) > 2:
# 		legr.AddEntry(ratio1,"reconstructed","l")
# #		legr.AddEntry(ratio1, "11", "l")
#  	line1 = ROOT.TLine(rangeminimum, 1, rangemaximum, 1)
# 	line11 = ROOT.TLine(rangeminimum, 1.1, rangemaximum, 1.1)
# 	line09 = ROOT.TLine(rangeminimum, 0.9, rangemaximum, 0.9)
# 	line1.SetLineStyle(1)
# 	line11.SetLineStyle(2)
# 	line09.SetLineStyle(2)
	


	can.cd()
	can.SetLogy(log)
	can.SaveAs(filename+".pdf")
	can.SaveAs(filename+".png")

	return


el_cuts = "ll_lepid[0] == 11 && el_pt[0] > 25 && el_pt[1] > 15"
mu_cuts = "ll_lepid[0] == 13 && mu_pt[0] > 20 && mu_pt[1] > 10"
baseline = "nphoton > 0 && ll_m[0] > 50 && llphoton_m[0]+ll_m[0]>185 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && photon_drmin[0] > 0.4"
#PlotHist(0, "", "mllg", "Number of Events", "mllg", 40, 100, 180, "SM_mllg_30_high_nopt_BDT0p05", [[[Ch], "llphoton_m[0]", "weight*(nphoton > 0 && ll_m[0] > 50 && llphoton_m[0]+ll_m[0]>185 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && photon_drmin[0] > 0.4 && BDT_score > 0.05)", "llg", "SM"]])
#PlotHist(0, "", "mllg", "Number of Events", "mllg", 40, 100, 180, "SM_mllg_30_high_zoom_highBDT", [[[Ch], "llphoton_m[0]", "weight*(BDT_score > 0.05 && "+baseline+" )","llg", "SM Zg"]]) 
#&& (("+el_cuts+") || ("+mu_cuts+")))", "llg", "SM Zg"]])
#PlotHist(0, "", "mllg", "Number of Events", "mllg", 40, 100, 180, "SM_mllg_30_low_nopt_zoom_highBDT", [[[Ch1], "llphoton_m[0]", "weight*(BDT_score > 0.01 && "+baseline+" )","llg", "SM Zg"]])
# && (("+el_cuts+") || ("+mu_cuts+")))", "llg", "SM Zg"]])

#PlotHist(0, "", "mllg", "Number of Events", "mllg", 25, 100, 180, "SM_mllg_original_nocosT",[[[Ch], "llphoton_m[0]", "weight*(BDT_score > 0.05 && "+baseline+")","llg", "SM Zg"]])
# && (("+el_cuts+") || ("+mu_cuts+")))", "llg", "SM Zg"]]) 

#PlotHist(0, "", "mllg", "Number of Events", "mllg", 40, 100, 180, "DY_mllg_in_pt_30_high_zoom", [[[Ch], "mllg", "photon_pt > 30", "llg", "DY+Zg"]]) 
#PlotHist(0, "", "mllg", "Number of Events", "mllg", 40, 100, 180, "DY_mllg_in_pt_30_low_zoom", [[[Ch], "mllg", "photon_pt <= 30", "llg", "DY+Zg"]])
