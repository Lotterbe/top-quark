import ROOT
import itertools
import math
import numpy as np
import array

from Analysis import Analysis
from Analysis import AnalysisHelpers
from Analysis import Constants
from cmath import sqrt as csqrt


#======================================================================

class TTbarAnalysis(Analysis.Analysis):
  """Semileptonic TTbarAnalysis loosely based on the ATLAS analyses of top pair events where
  one W boson decays to leptons and one decays to hadrons.
  """
  def __init__(self, store):
    super(TTbarAnalysis, self).__init__(store)

  def initialize(self):
      self.hist_WtMass      =  self.addStandardHistogram("WtMass")

      self.hist_leptn       =  self.addStandardHistogram("lep_n")
      #Erwähne ich später
      #self.hist_leptpt      =  self.addStandardHistogram("lep_pt")
      #self.hist_lepteta     =  self.addStandardHistogram("lep_eta")
      self.hist_leptE       =  self.addStandardHistogram("lep_E")
      self.hist_leptphi     =  self.addStandardHistogram("lep_phi")
      self.hist_leptch      =  self.addStandardHistogram("lep_charge")
      self.hist_leptID      =  self.addStandardHistogram("lep_type")
      self.hist_leptptc     =  self.addStandardHistogram("lep_ptconerel30")
      self.hist_leptetc     =  self.addStandardHistogram("lep_etconerel20")
      self.hist_lepz0       =  self.addStandardHistogram("lep_z0")
      self.hist_lepd0       =  self.addStandardHistogram("lep_d0")

      self.hist_nbtags      =  self.addStandardHistogram("n_btags")
      self.hist_njets       =  self.addStandardHistogram("n_jets")
      self.hist_jetspt      =  self.addStandardHistogram("jet_pt")
      self.hist_jetm        =  self.addStandardHistogram("jet_m")
      self.hist_jetJVF      =  self.addStandardHistogram("jet_jvf")
      self.hist_jeteta      =  self.addStandardHistogram("jet_eta")
      self.hist_jetmv1      =  self.addStandardHistogram("jet_MV1")

      self.hist_etmiss      = self.addStandardHistogram("etmiss")
      self.hist_vxp_z       = self.addStandardHistogram("vxp_z")
      self.hist_pvxp_n      = self.addStandardHistogram("pvxp_n")

      #heißt bei uns jetzt anders
      #self.hist_topmass     = self.addStandardHistogram("top_mass")

      #Eigene Histogramme
      #Massen
      self.hist_hadr_topmass   =  self.addStandardHistogram("hadr_topmass")
      self.hist_semilep_topmass     =  self.addStandardHistogram("semilep_topmass")
      self.hist_Wmass       = self.addStandardHistogram("W_mass")
      #etas
      self.hist_hadr_topeta   =  self.addStandardHistogram("hadr_topeta")
      self.hist_semilep_topeta=  self.addStandardHistogram("semilep_topeta")
      self.hist_semilep_weta   =  self.addStandardHistogram("semilep_weta")
      self.hist_lepeta  =  self.addStandardHistogram("lepeta")
      #transversal impuls
      self.hist_hadr_topTM    =  self.addStandardHistogram("hadr_topTM")
      self.hist_semilep_topTM =  self.addStandardHistogram("semilep_topTM")
      self.hist_semilep_wTM   =  self.addStandardHistogram("semilep_wTM")
      self.hist_lepTM   =  self.addStandardHistogram("lepTM")
      self.hist_totTM   =  self.addStandardHistogram("totTM")
      #com
      self.hist_com_otherjets =  self.addStandardHistogram("com_otherjets")
      self.hist_com_bjets     =  self.addStandardHistogram("com_bjets")
      self.hist_com_tot =  self.addStandardHistogram("com_tot")
      
      #Variablen Namen für das output tuple und das erstellen dieser
      self.varnames = ["HadronicTopMass", "SemilepTopMass", "SemilepWMass", "HadronicTopEta", "SemilepTopEta", "SemilepWEta", \
      "LepEta", "HadronicTopTransMom", "SemilepTopTransMom", "SemilepWTransMom", "LepTransMom", "TotalTransMom",\
      "COMOtherJets", "COMbJets", "COMTotal"]
      self.varnames.append("weight")
      self.output_tuple = ROOT.TNtuple("tuple", "tuple", ":".join(self.varnames))

#z-Komponente des Neutrino-Impulses
  def NeutrinoZMomentum(self, lepton_vec, neutrino_trans_momentum_x, neutrino_trans_momentum_y):
      neutrino_px = neutrino_trans_momentum_x
      neutrino_py = neutrino_trans_momentum_y
      lepton = lepton_vec
      #Annahme für die die Berechnung des pz vum Neutrino gilt 
      mw = 80.42
      #mw = 0
      mu = mw / 2 + (neutrino_px * lepton.Px() + neutrino_py * lepton.Py())
      #Umformung 
      neutrino_TM_squ = neutrino_px**2 + neutrino_py**2
      #zMomentPos = (mu * lepton.Pz() + lepton.E() * np.lib.scimath.sqrt((lepton.Pz()**2 - lepton.E()**2) * neutrino_TM_squ \
      #      + mu**2)) / (lepton.E()**2 - lepton.Pz()**2) 
      #print('Pos' + str(zMomentPos))
      #zMomentNeg = (mu * lepton.Pz() - lepton.E() * np.lib.scimath.sqrt((lepton.Pz()**2 - lepton.E()**2) * neutrino_TM_squ \
      #      + mu**2)) / (lepton.E()**2 - lepton.Pz()**2) 
      #straight aus dem skript
      zMomentPos = (mu * lepton.Pz()) / (lepton.E()**2 - lepton.Pz()**2) \
      + np.lib.scimath.sqrt(((mu * lepton.Pz()) / (lepton.E()**2 - lepton.Pz()**2))**2 \
      - (lepton.E()**2 * neutrino_TM_squ - mu**2) / (lepton.E()**2 - lepton.Pz()**2))
      zMomentNeg = (mu * lepton.Pz()) / (lepton.E()**2 - lepton.Pz()**2) \
      - np.lib.scimath.sqrt(((mu * lepton.Pz()) / (lepton.E()**2 - lepton.Pz()**2))**2 \
      - (lepton.E()**2 * neutrino_TM_squ - mu**2) / (lepton.E()**2 - lepton.Pz()**2))
      #print('Neg' + str(zMomentNeg))
      #vorher
      #mu = mw**2/2 + math.cos(phi) * (lepton.Px()*neutrino_px + lepton.Py()*neutrino_py)
      #neutrino_trans_squ = neutrino_px**2 + neutrino_py**2
      #zMomentPos = (lepton.Pz()*mu + lepton.E()*((lepton.E()**2+lepton.Pz()**2)*neutrino_trans_squ\
      #            +mu**2)**(1/2)) / (lepton.E()**2-lepton.Pz()**2)
      #zMomentNeg = (lepton.Pz()*mu - lepton.E()*((lepton.E()**2+lepton.Pz()**2)*neutrino_trans_squ\
      #            +mu**2)**(1/2)) / (lepton.E()**2-lepton.Pz()**2)
      #Filtern, welches die 'richtige' Lösung ist
      #print('Pos' + str(type(zMomentPos)))
      #print('Neg' + str(type(zMomentNeg))) 

      if zMomentPos.imag != 0.0:
            #print('Henlo, bin in Fall 1')
            if zMomentNeg.imag != 0.0:
                  if np.abs(zMomentNeg.real) > np.abs(zMomentPos.real):
                        zMoment = zMomentPos.real
                  else:
                        zMoment = zMomentNeg.real
            else:
                  zMoment = zMomentNeg
      elif zMomentNeg.imag != 0.0:
            #print('Henlo, bin in Fall 2')
            zMoment = zMomentPos
      else:
            #print('Henlo, bin in Fall 3')
            if np.abs(zMomentNeg) > np.abs(zMomentPos):
                  zMoment = zMomentPos
            else:
                  zMoment = zMomentNeg
      #print(zMoment)
      return zMoment

# Definierung von chi² Funktion
  def chi_squared(self, hadronic_w, semilep_top, hadronic_top):
      #values
      m_w_jjexp = 79.3 #GeV
      sigma_mwjjexp = 18.0 #GeV
      sigma_mtexp = 32.1 #GeV
      #Calculation of chi squared
      chi_squ = (hadronic_w.M() - m_w_jjexp)**2 / (sigma_mwjjexp**2) +\
             (semilep_top.M() - hadronic_top.M())**2 / (sigma_mtexp**2)
      return chi_squ

#transverse momentum
  def transverse_momentum(self, particle_tlv):
      transmom = (particle_tlv.Px()**2 + particle_tlv.Py()**2)**0.5
      return transmom 
#Pseudorapidität
  def part_eta(self, particle_tlv):
      eta_particle = - np.log(np.tan(np.abs(particle_tlv.Theta()) / 2))
      return eta_particle

  def calc_kin_obs(self, poss_list, leadlepton):
      #Aufschlüsselung Teilchen: 
      #[0]=semilep_top, [1]=lepton, [2]=neutrino, [3]=lep_w, [4]=semilep_b, [5]=hadr_top, [6]=qjet, [7]=qbarjet
      #[8]=hadr_w, [9]=hadr_b
      #14 verschiedene kinematische Observablen
      #Massen
      hadr_topmass = poss_list[5].M()
      semilep_topmass = poss_list[0].M()
      #eta
      hadr_topeta = poss_list[5].PseudoRapidity()
      semilep_topeta = poss_list[0].PseudoRapidity()
      semilep_weta = poss_list[3].PseudoRapidity()
      lepeta = leadlepton.eta()
      #Transversal impuls
      hadr_topTM = self.transverse_momentum(poss_list[5])
      semilep_topTM = self.transverse_momentum(poss_list[0])
      semilep_wTM = self.transverse_momentum(poss_list[3])
      lepTM = leadlepton.pt()
      #gesamt transversal impuls
      tot_part = poss_list[0]
      for particle in poss_list[1:]:
            tot_part += particle
      totTM = self.transverse_momentum(tot_part)
      #center of mass energy
      com_otherjets = poss_list[6].M() + poss_list[7].M()
      com_bjets = poss_list[4].M() + poss_list[9].M()
      com_tot = tot_part.M()
      return hadr_topmass, semilep_topmass, hadr_topeta, semilep_topeta, semilep_weta, lepeta, hadr_topTM, \
                   semilep_topTM, semilep_wTM, lepTM, totTM, com_otherjets, com_bjets, com_tot


  def analyze(self):
      # retrieving objects
      eventinfo = self.Store.getEventInfo()
      #print(eventinfo)
      weight = eventinfo.scalefactor()*eventinfo.eventWeight() if not self.getIsData() else 1
      self.countEvent("all", weight)

      # apply standard event based selection
      if not AnalysisHelpers.StandardEventCuts(eventinfo): return False
      self.countEvent("EventCuts", weight)

      # neutrinos are expected, so cut on missing transverse momentum
      etmiss = self.Store.getEtMiss()
      if not (etmiss.et() > 30.0): return False
      self.countEvent("MET", weight)

      # one good lepton from one of the W boson decays is expected, so require exactly one good lepton
      goodLeptons = AnalysisHelpers.selectAndSortContainer(self.Store.getLeptons(), AnalysisHelpers.isGoodLepton, lambda p: p.pt())
      if not (len(goodLeptons) == 1): return False
      self.countEvent("1 Lepton", weight)

      leadlepton = goodLeptons[0]

      ##########################################################################
      ### Schnitte auf Anzahl der Jets und b-tags

      # two jets from one of the W boson decays as well as two b-jets from the top pair decays are expected
      goodJets = AnalysisHelpers.selectAndSortContainer(self.Store.getJets(), AnalysisHelpers.isGoodJet, lambda p: p.pt())
      if not len(goodJets) == 4: return False # >= 2 vorher
      self.countEvent("Jets", weight)

      # apply the b-tagging requirement using the MV1 algorithm at 80% efficiency
      btags = sum([1 for jet in goodJets if jet.mv1() > 0.7892])
      if not (btags == 1): return False # >=0 vorher
      self.countEvent("btags", weight)


      # apply a cut on the transverse mass of the W boson decaying to leptons
      if not (AnalysisHelpers.WTransverseMass(leadlepton, etmiss) > 30.0): return False

      ##########################################################################
      ## Aufgabe 2
      zmom = self.NeutrinoZMomentum(leadlepton.tlv(), etmiss.tlv().Px(), etmiss.tlv().Py())

      etmiss.tlv().SetPz(zmom) 
      Wleptonic = leadlepton.tlv() + etmiss.tlv()

      # finde heraus welche jets in goodJets der getaggte b-Jet ist 
      other_jets = []
      b_jet_tag = []
      #btags herausfinden 
      for i in range(0, 4):
            if goodJets[i].mv1() > 0.7892:
                  b_jet_tag.append(i)
            else:
                  other_jets.append(i)

      #Zuordnung der einzelnen Jets mithilfe der Permutationen

      #Überprüfung, welches der jets als bjet getaggt wurde
      btagged_jet = goodJets[b_jet_tag[0]]
      first_other_jet = goodJets[other_jets[0]]
      second_other_jet = goodJets[other_jets[1]]
      third_other_jet = goodJets[other_jets[2]]
      #Berechnung der chi² der verschiedenen Permutationen
      #Listen um alle Permutationen abzuspeichern und später darauf zuzugreifen
      #6 Permutationen und 10 teilnehmende Teilchen
      #Aufschlüsselung Teilchen: 
      #[0]=semilep_top, [1]=lepton, [2]=neutrino, [3]=lep_w, [4]=semilep_b, [5]=hadr_top, [6]=qjet, [7]=qbarjet
      #[8]=hadr_w, [9]=hadr_b
      for i in range(0, 6):
            #semileptonischer t-Zerfall
            poss13_semilep_t = Wleptonic + btagged_jet.tlv()
            if i == 0:
                  poss1_hadr_w = first_other_jet.tlv() + second_other_jet.tlv()
                  poss1_hadr_t = poss1_hadr_w + third_other_jet.tlv()
                  poss1_chi_squ = self.chi_squared(poss1_hadr_w, poss13_semilep_t, poss1_hadr_t)
                  poss1_list = [poss13_semilep_t, leadlepton.tlv(), etmiss.tlv(),\
                   Wleptonic, btagged_jet.tlv(), poss1_hadr_t, first_other_jet.tlv(), second_other_jet.tlv(),\
                   poss1_hadr_w, third_other_jet.tlv()]
            if i == 1:
                  poss2_hadr_w = second_other_jet.tlv() + third_other_jet.tlv()
                  poss2_hadr_t = poss2_hadr_w + first_other_jet.tlv()
                  poss2_chi_squ = self.chi_squared(poss2_hadr_w, poss13_semilep_t, poss2_hadr_t)
                  poss2_list = [poss13_semilep_t, leadlepton.tlv(), etmiss.tlv(),\
                   Wleptonic, btagged_jet.tlv(), poss2_hadr_t, second_other_jet.tlv(), third_other_jet.tlv(),\
                   poss2_hadr_w, first_other_jet.tlv()]
            if i == 2:
                  poss3_hadr_w = third_other_jet.tlv() + first_other_jet.tlv()
                  poss3_hadr_t = poss3_hadr_w + second_other_jet.tlv()
                  poss3_chi_squ = self.chi_squared(poss3_hadr_w, poss13_semilep_t, poss3_hadr_t)
                  poss3_list = [poss13_semilep_t, leadlepton.tlv(), etmiss.tlv(),\
                   Wleptonic, btagged_jet.tlv(), poss3_hadr_t, third_other_jet.tlv(), first_other_jet.tlv(),\
                   poss3_hadr_w, second_other_jet.tlv()]
            #hadronischer t-Zerfall
            if i == 3:
                  poss4_semilep_t = Wleptonic + third_other_jet.tlv()
                  poss4_hadr_w = first_other_jet.tlv() + second_other_jet.tlv()
                  poss4_hadr_t = poss4_hadr_w + btagged_jet.tlv()
                  poss4_chi_squ = self.chi_squared(poss4_hadr_w, poss4_semilep_t, poss4_hadr_t)
                  poss4_list = [poss4_semilep_t, leadlepton.tlv(), etmiss.tlv(),\
                   Wleptonic, third_other_jet.tlv(), poss4_hadr_t, first_other_jet.tlv(), second_other_jet.tlv(),\
                   poss4_hadr_w, btagged_jet.tlv()]
            if i == 4:
                  poss5_semilep_t = Wleptonic + first_other_jet.tlv()
                  poss5_hadr_w = second_other_jet.tlv() + third_other_jet.tlv()
                  poss5_hadr_t = poss5_hadr_w + btagged_jet.tlv()
                  poss5_chi_squ = self.chi_squared(poss5_hadr_w, poss5_semilep_t, poss5_hadr_t)
                  poss5_list = [poss5_semilep_t, leadlepton.tlv(), etmiss.tlv(),\
                   Wleptonic, first_other_jet.tlv(), poss5_hadr_t, second_other_jet.tlv(), third_other_jet.tlv(),\
                   poss5_hadr_w, btagged_jet.tlv()]
            if i == 5:
                  poss6_semilep_t = Wleptonic + second_other_jet.tlv()
                  poss6_hadr_w = third_other_jet.tlv() + first_other_jet.tlv()
                  poss6_hadr_t = poss6_hadr_w + btagged_jet.tlv()
                  poss6_chi_squ = self.chi_squared(poss6_hadr_w, poss6_semilep_t, poss6_hadr_t)
                  poss6_list = [poss6_semilep_t, leadlepton.tlv(), etmiss.tlv(),\
                   Wleptonic, second_other_jet.tlv(), poss6_hadr_t, third_other_jet.tlv(), first_other_jet.tlv(),\
                   poss6_hadr_w, btagged_jet.tlv()]

      #Auswahl des chi² mit dem kleinsten Wert, d.h. Auswahl der richtigen Zuordnung der bJets
      all_chi_squ = [poss1_chi_squ, poss2_chi_squ, poss3_chi_squ, poss4_chi_squ, poss5_chi_squ, poss6_chi_squ]
      min_chi_squ = min(all_chi_squ)
      min_index = all_chi_squ.index(min_chi_squ)

      # Leptonic W boson
      Wmass = Wleptonic.M() # vorher = 80
      

      #Kinematische Größen

      #calculation of the kinematic observables with the right chi-assignment
      hadr_topmass = 0
      semilep_topmass = 0
      hadr_topeta = 0 
      semilep_topeta = 0 
      semilep_weta = 0
      lepeta = 0
      #Transversal impuls
      hadr_topTM = 0
      semilep_topTM = 0
      semilep_wTM = 0
      lepTM = 0
      #gesamt transversal impuls
      totTM = 0
      #center of mass energy
      com_otherjets = 0
      com_bjets = 0
      com_tot = 0

      if min_index == 0:
            #poss1
            hadr_topmass, semilep_topmass, hadr_topeta, semilep_topeta, semilep_weta, lepeta, hadr_topTM, \
                   semilep_topTM, semilep_wTM, lepTM, totTM, com_otherjets, com_bjets, com_tot = self.calc_kin_obs(poss1_list, leadlepton)
      if min_index == 1:
            #poss2
            hadr_topmass, semilep_topmass, hadr_topeta, semilep_topeta, semilep_weta, lepeta, hadr_topTM, \
                   semilep_topTM, semilep_wTM, lepTM, totTM, com_otherjets, com_bjets, com_tot = self.calc_kin_obs(poss2_list, leadlepton)
      if min_index == 2:
            #poss3
            hadr_topmass, semilep_topmass, hadr_topeta, semilep_topeta, semilep_weta, lepeta, hadr_topTM, \
                   semilep_topTM, semilep_wTM, lepTM, totTM, com_otherjets, com_bjets, com_tot = self.calc_kin_obs(poss3_list, leadlepton)
      if min_index == 3:
            #poss4
            hadr_topmass, semilep_topmass, hadr_topeta, semilep_topeta, semilep_weta, lepeta, hadr_topTM, \
                   semilep_topTM, semilep_wTM, lepTM, totTM, com_otherjets, com_bjets, com_tot = self.calc_kin_obs(poss4_list, leadlepton)
      if min_index == 4:
            #poss5
            hadr_topmass, semilep_topmass, hadr_topeta, semilep_topeta, semilep_weta, lepeta, hadr_topTM, \
                   semilep_topTM, semilep_wTM, lepTM, totTM, com_otherjets, com_bjets, com_tot = self.calc_kin_obs(poss5_list, leadlepton)
      if min_index == 5:
            #poss6
            hadr_topmass, semilep_topmass, hadr_topeta, semilep_topeta, semilep_weta, lepeta, hadr_topTM, \
                   semilep_topTM, semilep_wTM, lepTM, totTM, com_otherjets, com_bjets, com_tot = self.calc_kin_obs(poss6_list, leadlepton)

      #Variablen Liste für nTuple
      values = [hadr_topmass, semilep_topmass, Wmass, hadr_topeta, semilep_topeta,\
      semilep_weta, lepeta, hadr_topTM, semilep_topTM, semilep_wTM, \
      lepTM, totTM, com_otherjets, com_bjets, com_tot, weight]

      #konvertiere die Variablen Liste in ein float array und fülle diese in das Output Tuple
      #self.output_tuple.Fill(arr_val)
      self.output_tuple.Fill(array.array('f', values))


      ##########################################################################

      ###########################################################################
      ### Ab hier findet das fuellen der Histogramme statt

      # Histograms detailing event information
      self.hist_vxp_z.Fill(eventinfo.primaryVertexPosition(), weight)
      self.hist_pvxp_n.Fill(eventinfo.numberOfVertices(), weight)

      # histograms for the W boson properties
      self.hist_WtMass.Fill(AnalysisHelpers.WTransverseMass(leadlepton, etmiss), weight)

      # histograms for missing et
      self.hist_etmiss.Fill(etmiss.et(), weight)

      # histograms detailing lepton information
      self.hist_leptn.Fill(len(goodLeptons), weight)
      #erwähne ich später
      #self.hist_leptpt.Fill(leadlepton.pt(), weight)
      #self.hist_lepteta.Fill(leadlepton.eta(), weight)
      self.hist_leptE.Fill(leadlepton.e(), weight)
      self.hist_leptphi.Fill(leadlepton.phi(), weight)
      self.hist_leptch.Fill(leadlepton.charge(), weight)
      self.hist_leptID.Fill(leadlepton.pdgId(), weight)
      self.hist_lepz0.Fill(leadlepton.z0(), weight)
      self.hist_lepd0.Fill(leadlepton.d0(), weight)
      self.hist_leptptc.Fill(leadlepton.isoptconerel30(), weight)
      self.hist_leptetc.Fill(leadlepton.isoetconerel20(), weight)

      # histograms detailing jet information
      self.hist_njets.Fill(len(goodJets), weight)
      self.hist_nbtags.Fill(btags, weight)

      [self.hist_jetm.Fill(jet.m(), weight) for jet in goodJets]
      [self.hist_jetspt.Fill(jet.pt(), weight) for jet in goodJets]
      [self.hist_jetJVF.Fill(jet.jvf(), weight) for jet in goodJets]
      [self.hist_jeteta.Fill(jet.eta(), weight) for jet in goodJets]
      [self.hist_jetmv1.Fill(jet.mv1(), weight) for jet in goodJets]


      #eigene Histogramme
      #print('klappt!\n#################')
      #masses 
      self.hist_hadr_topmass.Fill(hadr_topmass, weight)
      self.hist_semilep_topmass.Fill(semilep_topmass, weight)
      self.hist_Wmass.Fill(Wmass, weight)
      #etas
      self.hist_hadr_topeta.Fill(hadr_topeta, weight)
      self.hist_semilep_topeta.Fill(semilep_topeta, weight)
      self.hist_semilep_weta.Fill(semilep_weta, weight)
      self.hist_lepeta.Fill(lepeta, weight)
      #transversal impulse
      self.hist_hadr_topTM.Fill(hadr_topTM, weight)
      self.hist_semilep_topTM.Fill(semilep_topTM, weight)
      self.hist_semilep_wTM.Fill(semilep_wTM, weight)
      self.hist_lepTM.Fill(lepTM, weight)
      self.hist_totTM.Fill(totTM, weight)
      #com
      self.hist_com_otherjets.Fill(com_otherjets, weight)
      self.hist_com_bjets.Fill(com_bjets, weight)
      self.hist_com_tot.Fill(com_tot, weight)
      
      return True

  def finalize(self):
      self.output_tuple.Write()

