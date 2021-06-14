import ROOT
import numpy as np
import iminuit as imin
import math
import os
import sys
from array import array
from Plotting import infofile
from Plotting import PlotStyle as PS

#quick and dirty scaling and adding of the other MC-Files
scaleDB  = infofile.infos
# vorher
#hname = "WtMass"
#SemilepWMass, W_mass
hname = "nn"

def getScaleFactor(scalingkey):
    entry = scaleDB[scalingkey]
    return 1000*entry["xsec"]/(entry["sumw"]*entry["red_eff"]) 

indir = "output/"


######### signal
fele = ROOT.TFile(indir +'h_data_Egamma.root')
hdata = fele.Get(hname)
hdata.Draw()
fmuo = ROOT.TFile(indir +'h_data_Muons.root')
hdata.Add(fmuo.Get(hname))


######### signal
fsig = ROOT.TFile(indir +'h_ttbar_lep.root')
hsig = fsig.Get(hname)
hsig.Scale(getScaleFactor("ttbar_lep"))

######### background

## Make an empty histogram with the axis properties of the data histogram
hbkg = ROOT.TH1F("bkg","bkg",hdata.GetNbinsX(),hdata.GetXaxis().GetXmin(),hdata.GetXaxis().GetXmax())

bkg = ["WenuWithB", "WmunuWithB", "WtaunuWithB", "stop_schan", "stop_tchan_antitop", "stop_tchan_top", "stop_wtchan", "WenuJetsBVeto", "WenuNoJetsBVeto", "WmunuJetsBVeto", "WmunuNoJetsBVeto", "WtaunuJetsBVeto", "WtaunuNoJetsBVeto", "WW", "WZ", "Zee", "Zmumu", "Ztautau", "ZZ"] #TODO add all other backgrounds  


for b in bkg: 
  f = ROOT.TFile(indir + "h_" + b + ".root") 
  h = f.Get(hname) 
  h.Scale(getScaleFactor(b))
  hbkg.Add(h)


nbins = hdata.GetNbinsX()+1
data = []
sig  = []
bkg = []
for i in range (1, nbins) :
    data.append(hdata.GetBinContent(i))
    sig.append(hsig.GetBinContent(i))
    bkg.append(hbkg.GetBinContent(i))
   
adata = np.array(data[0:])
asig = np.array(sig[0:])
abkg = np.array(bkg[0:])

# Calculate expected number of events
sig_exp = asig.sum()
bkg_exp = abkg.sum()

# Normalise templates to unity
asig = asig/asig.sum()
abkg = abkg/abkg.sum()


print("Number of events in data:", adata.sum())
print("Signal expectation: ", sig_exp)
print("Background expectation:", bkg_exp)

############# Likelihood function #######################################
def bllh(s_sig,s_bkg):
    pred = s_bkg*bkg_exp*abkg + s_sig*sig_exp*asig 
    obs  = adata
    ll = 0
    for i in range(len(adata)):
        if (pred[i] > 0):
            ll += -pred[i] + obs[i]*math.log(pred[i])
    return (-ll)

bllh.errordef=imin.Minuit.LIKELIHOOD

m = imin.Minuit(bllh, s_sig=1, s_bkg=1)

m.migrad()  # run optimiser

print("Signal scale factor {0:3.2f} +- {1:3.2f}".format(m.values[0],m.errors[0]))
print("Background scale factor {0:3.2f} +- {1:3.2f}".format(m.values[1],m.errors[1]))


PS.setStyle()
canvas = ROOT.TCanvas("canvas","canvas",900,900)
canvas.cd()


hsig.Scale(m.values[0])
hsig.SetFillColor(ROOT.kRed)
hbkg.Scale(m.values[0])
hbkg.SetFillColor(ROOT.kYellow)

htotal = ROOT.THStack("stack","Neural network output")
htotal.Add(hbkg)
htotal.Add(hsig)

htotal.Draw("histo")

hdata.SetMarkerStyle(20)
hdata.SetMarkerSize(1)
hdata.Draw("p same")

#canvas.SaveAs("nn_fit.png")

############# Ratio plot #######################################

ROOT.gPad.SetBottomMargin(0.3)
ratiopad = ROOT.TPad("ratio","ratio",0,0,1,0.25)
ratiopad.SetTopMargin(0)
ratiopad.SetBottomMargin(0.3)
ratiopad.Draw()
ratiopad.cd()

ratio = hdata.Clone("ratio")
ratio.Divide(htotal.GetStack().Last())


ratio.SetNdivisions(205,"y")
ratio.SetNdivisions(505,"x")
ratio.SetLabelSize(0.1,"x")
ratio.SetLabelSize(0.1,"y")
ratio.SetTitleSize(0.1,"x")
ratio.SetTitleSize(0.1,"y")

ratio.SetTitle("")
ratio.GetYaxis().SetTitle("data / prediction")
ratio.Draw("p")

oneline = ROOT.TLine(ratio.GetXaxis().GetXmin(),1,ratio.GetXaxis().GetXmax(),1)
oneline.Draw()

canvas.SaveAs("nn_fit_pimped.png")


