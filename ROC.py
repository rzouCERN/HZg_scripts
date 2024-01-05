#!/usr/bin/env python
import math
import matplotlib
#%matplotlib inline
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import ROOT
ROOT.gStyle.SetOptStat(1111)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.TH1.SetDefaultSumw2(True)
#ROOT.gStyle.SetErrorY(0)

import optparse
import os
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
#dir = os.path.dirname(os.path.abspath(__file__))
dir = "Output_Run2_fixed_photon_pt_ratio"

# inTree = ROOT.TChain("outtree")
# inTrees = ROOT.TChain("outtree")

# inTree.AddFile(dir+"/GGF_2016_output.root")
# inTree.AddFile(dir+"/GGF_2017_output.root")
# inTree.AddFile(dir+"/GGF_2018_output.root")
#histFileName1 = dir+"/GGF_output.root"
histFileName1 = dir+"/GGF_2017_output.root"
histFileName2 = dir+"/BKG_2017_output.root"
infile = ROOT.TFile.Open(histFileName1 ,"READ")
infile2 = ROOT.TFile.Open(histFileName2 ,"READ")

inTree = infile.Get("outtree") #HiggsPt
inTree2 = infile2.Get("outtree") #HiggsPt

entries = inTree.GetEntriesFast()
entries2 = inTree2.GetEntriesFast()
print entries
print entries2
name = "_Run2_fixed_photon_pt_ratio"


histchain = []
can = ROOT.TCanvas()
leg = ROOT.TLegend(0.75,0.7,0.88,0.85)
color = [ROOT.kOrange, ROOT.kBlue, ROOT.kRed+4, ROOT.kRed, ROOT.kViolet+4, ROOT.kViolet, ROOT.kBlack, ROOT.kGray+3, ROOT.kPink, ROOT.kTeal, ROOT.kOrange]


nbins = 50
BDT_signal = ROOT.TH1F("BDT Signal",";Number of events;BDT;",nbins,-1,0.5)
BDT_signal.Sumw2()
BDT_background = ROOT.TH1F("BDT Background",";Number of events;BDT;",nbins,-1,0.5)
BDT_background.Sumw2()


#el_trigger = " && HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL == 1"
#mu_trigger = "&& (HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 == 1 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 == 1)"
el_cuts = "ll_lepid[0] == 11 && el_pt[0] > 25 && el_pt[1] > 15"
mu_cuts = "ll_lepid[0] == 13 && mu_pt[0] > 20 && mu_pt[1]> 10"
baseline = "nphoton > 0 && ll_m[0] > 50 && llphoton_m[0]+ll_m[0]>185 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && photon_drmin[0] > 0.4 && photon_pt[0] > 15 && ll_m[0] > 80 && ll_m[0] < 100"
narrowrange = "llphoton_m[0] > 120 && llphoton_m[0] < 130"

# for jentry in range(entries): #entries
#     ientry = inTree.LoadTree( jentry)
#     if ientry < 0:
#         break
#     nb = inTree.GetEntry( jentry)

#     if nb <= 0:
#         continue

#     #print inTree_op.h[0]
#     if (inTree.el_trig == 0) & (inTree.mu_trig == 0):
#         continue

#     if not (inTree.ll_lepid[0] == 11 or inTree.ll_lepid[0] == 13):
#         continue
#     if (inTree.ll_lepid[0] == 11):
#         if not (inTree.el_pt[inTree.ll_i1[0]]> 25 and inTree.el_pt[inTree.ll_i2[0]] > 15):
#             continue
#     if (inTree.ll_lepid[0] == 13):
#         if not (inTree.mu_pt[inTree.ll_i1[0]]> 20 and inTree.mu_pt[inTree.ll_i2[0]] > 10):
#             continue
#     if not (inTree.nphoton > 0):
#         continue
#     if not (inTree.ll_m[0] > 50 and inTree.llphoton_m[0]+inTree.ll_m[0]>130 and inTree.llphoton_m[0] > 120 and inTree.llphoton_m[0] < 180 and inTree.photon_drmin[0] > 0.4 and inTree.photon_pt[0] > 15):
#         continue
#     BDT_signal.Fill(inTree.BDT_score) #h[0] for "op"

# for jentry in range(entries2): #entries
#     ientry = inTree2.LoadTree( jentry)
#     if ientry < 0:
#         break
#     nb = inTree2.GetEntry( jentry)

#     if nb <= 0:
#         continue

#     if (inTree2.el_trig == 0) and (inTree2.mu_trig == 0):
#         continue
#     if not (inTree2.ll_lepid[0] == 11 or inTree2.ll_lepid[0] == 13):
#         continue
#     if (inTree2.ll_lepid[0] == 11):
#         if not (inTree2.el_pt[inTree2.ll_i1[0]]> 25 and inTree2.el_pt[inTree2.ll_i2[0]] > 15):
#             continue
#     if (inTree2.ll_lepid[0] == 13):
#         if not (inTree2.mu_pt[inTree2.ll_i1[0]]> 20 and inTree2.mu_pt[inTree2.ll_i2[0]] > 10):
#             continue
#     if not (inTree2.nphoton > 0):
#         continue    
#     if not (inTree2.ll_m[0] > 50 and inTree2.llphoton_m[0]+inTree2.ll_m[0]>130 and inTree2.llphoton_m[0] > 120 and inTree2.llphoton_m[0] < 180 and inTree2.photon_drmin[0] > 0.4 and inTree2.photon_pt[0] > 15):
#         continue
#     #print inTree_op.h[0]
#     BDT_background.Fill(inTree2.BDT_score)

    
#inTree.Draw("BDT_score >> BDT Signal","w_lumi*w_year*(((el_trig == 1 && "+el_cuts+") || (mu_trig ==1 && "+mu_cuts+")) && ("+baseline+") && ("+narrowrange+"))")
#inTree2.Draw("BDT_score >> BDT Background","w_lumi*w_year*(((el_trig == 1 && "+el_cuts+") || (mu_trig ==1 && "+mu_cuts+")) && ("+baseline+") && ("+narrowrange+"))")
inTree.Draw("BDT_score >> BDT Signal","w_lumi*w_year*( ("+baseline+") && ("+narrowrange+"))")
inTree2.Draw("BDT_score >> BDT Background","w_lumi*w_year*(("+baseline+") && ("+narrowrange+"))")

norm = BDT_signal.Integral()
BDT_signal.Scale(1/norm)        
norm2 = BDT_background.Integral()
BDT_background.Scale(1/norm2)
histchain.append(BDT_signal)
histchain.append(BDT_background)
#print BDT_signal.GetNbinsX()
#print BDT_background.GetNbinsX()
#print BDT_signal.Integral(0,50)


FSR = [] #False signal rate
BR = [] #Background rate
TSR = [] #True signal rate
threshold = []
FSR_ref = []
BR_ref = []
for ibin in range(nbins):
    #print "bin center ttH ",BDT_signal.GetBinCenter(ibin+1)
    #print "bin center ttW ",BDT_background.GetBinCenter(ibin+1)
    nSig = BDT_signal.Integral(ibin,nbins)
    #print "ttH ",nttH
    nBkg = BDT_background.Integral(0,ibin+1)
#    nBkg = BDT_background.Integral(ibin,nbins)
    #print "ttW ",nttW

    BR.append(nBkg)
    TSR.append(nSig)
    BR_ref.append(ibin/(nbins-1))
    threshold.append(BDT_signal.GetBinCenter(ibin+1))

TSR_ref = BR_ref

print np.trapz(y=BR,x=TSR)
plotname = 'AUC : '+str(abs(round(np.trapz(y=BR,x=TSR),3)))
print plotname

fig = plt.figure()
ax = fig.add_subplot(111)
#ax.plot(BR,TSR,'r',label=plotname)
#ax.plot(BR_ref,TSR_ref,'k-.')
ax.plot(TSR,BR,'r',label=plotname)
#ax.plot(TSR_ref,BR_ref,'k-.')
plt.ylabel("Background Rejection")
plt.xlabel("Signal Efficiency")
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.legend()
fig.savefig('roc'+name+'.png')

# plt.plot(FSR,TSR,'r',label=plotname)
# plt.plot(FSR_ref,TSR_ref,'k-.')
# plt.xlabel("Fake Signal Rate")
# plt.ylabel("True Signal Rate")
# plt.legend()
# plt.close()
# #plt.show()





c = 1
print len(histchain)
histchain[0].Draw("HIST")
histchain[0].SetFillColor(color[0])
histchain[0].GetXaxis().SetLabelSize(0.05)
histchain[0].GetXaxis().SetTitleSize(0.05)
histchain[0].GetYaxis().SetLabelSize(0.05)
histchain[0].GetYaxis().SetTitleSize(0.05)

leg.AddEntry(histchain[0],histchain[0].GetName(),"f")
if len(histchain)>1:
    for histo2 in histchain[1:]:
        histo2.Draw("HIST l same")
        histo2.SetLineColor(color[c])
        histo2.SetLineWidth(3)
        leg.AddEntry(histo2,histo2.GetName(),"l")
        c+=1
leg.Draw("same")    
can.SaveAs("BDT_dist"+name+".png")


