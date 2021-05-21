"""This file defines standard histograms which can be reused in various analyses.
The ranges of these histograms should accomodate most analyses.
"""

import ROOT

def getStandardHistogram(name):
    if (name == "vxp_z"):           return ROOT.TH1D("vxp_z",           "Primary Vertex Position; z_{Vertex}; Events", 40, -200,200)
    if (name == "pvxp_n"):          return ROOT.TH1D("pvxp_n",          "Number of Vertices; N_{vertex}; Events", 30, -0.5,29.5)
    if (name == "etmiss"):          return ROOT.TH1D("etmiss",          "Missing Transverse Momentum;E_{T,Miss} [GeV];Events", 20, 0,200)
                                                                        
    if (name == "n_jets"):          return ROOT.TH1D("n_jets",          "Number of Jets;N_{jets};Events", 10, -0.5, 9.5)
    if (name == "n_btags"):         return ROOT.TH1D("n_btags",         "Number of b-tags;N_{btags};Events", 10, -0.5, 9.5)
    if (name == "jet_pt"):          return ROOT.TH1D("jet_pt",          "Jet Transverse Momentum;p_{T}^{jet} [GeV];Jets", 20, 0, 200)
    if (name == "jet_m"):           return ROOT.TH1D("jet_m",           "Jet Mass; m^{jet} [MeV]; Jets", 20, 0, 20000)
    if (name == "jet_jvf"):         return ROOT.TH1D("jet_jvf",         "Jet Vertex Fraction; JVF ; Jets", 20, 0, 1)
    if (name == "jet_eta"):         return ROOT.TH1D("jet_eta",         "Jet Pseudorapidity; #eta^{jet}; Jets", 30, -3, 3)
    if (name == "jet_MV1"):         return ROOT.TH1D("jet_MV1",         "Jet MV1; MV1 weight ; Jets", 20, 0, 1)
    if (name == "lep_n"):           return ROOT.TH1D("lep_n",           "Number of Leptons; N_{lep} ;Events", 10, -0.5, 9.5)
                                                                        
    #Erwähne ich später bei eigenen Variablen
    #if (name == "lep_eta"):         return ROOT.TH1D("lep_eta",         "Lepton Pseudorapidity; #eta^{lep}; Leptons", 30, -3, 3)
    #if (name == "lep_pt"):          return ROOT.TH1D("lep_pt",          "Lepton Transverse Momentum;p_{T}^{lep} [GeV];Leptons", 20, 0, 200)
    if (name == "lep_E"):           return ROOT.TH1D("lep_E",           "Lepton Energy; E^{lep} [GeV]; Leptons", 30, 0, 300)
    if (name == "lep_phi"):         return ROOT.TH1D("lep_phi",         "Lepton Azimuthal Angle ; #phi^{lep}; Leptons", 32, -3.2, 3.2)
    if (name == "lep_charge"):      return ROOT.TH1D("lep_charge",      "Lepton Charge; Q^{lep}; Leptons", 7, -1.75, 1.75)
    if (name == "lep_type"):        return ROOT.TH1D("lep_type",        "Lepton Absolute PDG ID; |PDG ID|^{lep}; Leptons", 31, -0.5, 30.5)
    if (name == "lep_ptconerel30"): return ROOT.TH1D("lep_ptconerel30", "Lepton Relative Transverse Momentum Isolation; ptconerel30^{lep}; Leptons", 40, -0.1, 0.4)
    if (name == "lep_etconerel20"): return ROOT.TH1D("lep_etconerel20", "Lepton Relative Transverse Energy Isolation; etconerel20^{lep}; Leptons", 40,  -0.1, 0.4)
    if (name == "lep_z0"):          return ROOT.TH1D("lep_z0",          "Lepton z0 impact parameter; z_{0}^{lep} [mm]; Leptons", 40, -1, 1)
    if (name == "lep_d0"):          return ROOT.TH1D("lep_d0",          "Lepton d0 impact parameter; d_{0}^{lep} [mm]; Leptons", 40, -1, 1)

    if (name == "leadlep_pt"):          return ROOT.TH1D("leadlep_pt",          "Leading Lepton Transverse Momentum;p_{T}^{leadlep} [GeV];Leptons", 20, 0, 200)
    if (name == "leadlep_eta"):         return ROOT.TH1D("leadlep_eta",         "Leading Lepton Pseudorapidity; #eta^{leadlep}; Leptons", 30, -3, 3)
    if (name == "leadlep_E"):           return ROOT.TH1D("leadlep_E",           "Leading Lepton Energy; E^{leadlep} [GeV]; Leptons", 30, 0, 300)
    if (name == "leadlep_phi"):         return ROOT.TH1D("leadlep_phi",         "Leading Lepton Azimuthal Angle ; #phi^{leadlep}; Leptons", 32, -3.2, 3.2)
    if (name == "leadlep_charge"):      return ROOT.TH1D("leadlep_charge",      "Leading Lepton Charge; Q^{leadlep}; Leptons", 7, -1.75, 1.75)
    if (name == "leadlep_type"):        return ROOT.TH1D("leadlep_type",        "Leading Lepton Absolute PDG ID; |PDG ID|^{leadlep}; Leptons",  31, -0.5, 30.5)
    if (name == "leadlep_ptconerel30"): return ROOT.TH1D("leadlep_ptconerel30", "Leading Lepton Relative Transverse Momentum Isolation; ptconerel30^{leadlep}; Leptons", 40, -0.1, 0.4)
    if (name == "leadlep_etconerel20"): return ROOT.TH1D("leadlep_etconerel20", "Leading Lepton Relative Transverse Energy Isolation; etconerel20^{leadlep}; Leptons", 40, -0.1, 0.4)
    if (name == "leadlep_z0"):          return ROOT.TH1D("leadlep_z0",          "Leading Lepton z0 impact parameter; z_{0}^{leadlep} [mm]; Leptons", 40, -1, 1)
    if (name == "leadlep_d0"):          return ROOT.TH1D("leadlep_d0",          "Leading Lepton d0 impact parameter; d_{0}^{leadlep} [mm]; Leptons", 40, -1, 1)
                                         
    if (name == "traillep_pt"):          return ROOT.TH1D("traillep_pt",          "Trailing Lepton Transverse Momentum;p_{T}^{traillep} [GeV];Leptons", 20, 0, 200)
    if (name == "traillep_eta"):         return ROOT.TH1D("traillep_eta",         "Trailing Lepton Pseudorapidity; #eta^{traillep}; Leptons", 30, -3, 3)
    if (name == "traillep_E"):           return ROOT.TH1D("traillep_E",           "Trailing Lepton Energy; E^{traillep} [GeV]; Leptons", 30, 0, 300)
    if (name == "traillep_phi"):         return ROOT.TH1D("traillep_phi",         "Trailing Lepton Azimuthal Angle ; #phi^{traillep}; Leptons", 32, -3.2, 3.2)
    if (name == "traillep_charge"):      return ROOT.TH1D("traillep_charge",      "Trailing Lepton Charge; Q^{traillep}; Leptons", 7, -1.75, 1.75)
    if (name == "traillep_type"):        return ROOT.TH1D("traillep_type",        "Trailing Lepton Absolute PDG ID; |PDG ID|^{traillep}; Leptons",  31, -0.5, 30.5)
    if (name == "traillep_ptconerel30"): return ROOT.TH1D("traillep_ptconerel30", "Trailing Lepton Relative Transverse Momentum Isolation; ptconerel30^{traillep} [GeV]; Leptons", 40, -0.1, 0.4)
    if (name == "traillep_etconerel20"): return ROOT.TH1D("traillep_etconerel20", "Trailing Lepton Relative Transverse Energy Isolation; etconerel20^{traillep} [GeV]; Leptons", 40, -0.1, 0.4)
    if (name == "traillep_z0"):          return ROOT.TH1D("traillep_z0",          "Trailing Lepton z0 impact parameter; z_{0}^{traillep} [mm]; Leptons", 40, -1, 1)
    if (name == "traillep_d0"):          return ROOT.TH1D("traillep_d0",          "Trailing Lepton d0 impact parameter; d_{0}^{traillep} [mm]; Leptons", 40, -1, 1)
                                        
    if (name == "WtMass"):            return ROOT.TH1D("WtMass",            "Transverse Mass of the W Candidate; M_{T,W} [GeV]; Events", 40, 0, 200)
    if (name == "invMass"):           return ROOT.TH1D("invMass",           "Invariant Mass of the Z Candidate;M_{ll} [GeV]; Events", 30, 60,120)
    
    #Eigene Variablen 
    #Massen
    if (name == "hadr_topmass"):          return ROOT.TH1D("hadr_topmass",            "Hadronic top quark mass; M(top) [GeV]; Events", 30, 100, 220)
    if (name == "semilep_topmass"):          return ROOT.TH1D("semilep_topmass",            "Semileptonic top quark mass; M(top) [GeV]; Events", 30, 100, 220)
    if (name == "W_mass"):            return ROOT.TH1D("W_mass",            "Semileptonic W boson mass; M(W) [GeV]; Events", 20, 40, 120)
    #if (name == "lep_mass"):            return ROOT.TH1D("lep_mass",            "Lepton mass; M(Lep) [GeV]; Events", 20, 40, 120)
    #eta
    if (name == "hadr_topeta"):    return ROOT.TH1D("hadr_topeta",    "Hadronic Top-Quark Pseudorapidity; #eta^{Hadr_top}; Hadronic Top-Quarks", 20, 40,120)
    if (name == "semilep_topeta"):    return ROOT.TH1D("semilep_topeta",    "Semileptonic Top-Quark Pseudorapidity; #eta^{Slep_top}; Semileptonic Top-Quarks ", 20, 40,120)
    if (name == "semilep_weta"):    return ROOT.TH1D("semilep_weta",    "Semileptonic W-Boson Pseudorapidity; #eta^{Slep_W}; Semileptonic W-Bosons", 20, 40,120)
    if (name == "lepeta"):         return ROOT.TH1D("lepeta",         "Lepton Pseudorapidity; #eta^{lep}; Leptons", 30, -3, 3)
    #Transversal Impuls
    if (name == "hadr_topTM"):    return ROOT.TH1D("hadr_topTM",    "Hadronic Top-Quark Transverse Momentum; p_{T}^{Hadr_top} [GeV];Hadronic Top-Quarks", 20, 40,120)
    if (name == "semilep_topTM"):    return ROOT.TH1D("semilep_topTM",    "Semileptonic Top-Quark Transverse Momentum; p_{T}^{Slep_top} [GeV];Semileptonic Top-Quarks", 20, 40,120)
    if (name == "semilep_wTM"):    return ROOT.TH1D("semilep_wTM",    "Semileptonic W-Boson Transverse Momentum; p_{T}^{Slep_W} [GeV];Semileptonic W-Bosons Top-Quark", 20, 40,120)
    if (name == "lepTM"):          return ROOT.TH1D("lepTM",          "Lepton Transverse Momentum;p_{T}^{lep} [GeV];Leptons", 20, 0, 200)
    #Gesamt TP
    if (name == "totTM"):          return ROOT.TH1D("totTM",          "Total Transverse Momentum;p_{T}^{tot} [GeV];Total Value", 20, 0, 200)
    #Center of mass E
    if (name == "com_otherjets"):    return ROOT.TH1D("com_otherjets",    "Other tagged Jets Center of Mass Energy; COM(Other Jets) [GeV];Other tagged Jets", 20, 40,120)
    if (name == "com_bjets"):    return ROOT.TH1D("com_bjets",    "B-tagged Jets Center of Mass Energy; COM(B-Jets) [GeV]; B-tagged Jets", 20, 40,120)
    if (name == "com_tot"):    return ROOT.TH1D("com_tot",    "Total Center of Mass Energy; COM(Total) [GeV];Total Value", 20, 40,120)
    


    return None
