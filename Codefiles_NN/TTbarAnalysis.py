import ROOT
import itertools
import math
import numpy as np

from Analysis import Analysis
from Analysis import AnalysisHelpers
from Analysis import Constants

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
      self.hist_leptpt      =  self.addStandardHistogram("lep_pt")
      self.hist_lepteta     =  self.addStandardHistogram("lep_eta")
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

      self.hist_Wmass       = self.addStandardHistogram("W_mass")
      self.hist_topmass     = self.addStandardHistogram("top_mass")


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
      

#Leas Schleife um Jets zu sortieren 
      #Whadronic = goodJets[Whadronic_indices[0]].tlv() + goodJets[Whadronic_indices[1]].tlv()
      #b1 = goodJets[b_indices[0]].tlv()
      #b2 = goodJets[b_indices[1]].tlv()
      #dr1 = Whadronic.DeltaR(b1)
      #dr2 = Whadronic.DeltaR(b2)
      #if dr1 < dr2:
      #	thadronic = Whadronic + b1
      #else:
      #	thadronic = Whadronic + b2

      #Transversalimpuls Neutrino
      #Transversalimpuls-Komponenten x und y des Neutrinos
      #TransMomNeutrino_x = - (goodJets[0].tlv().Px() + goodJets[1].tlv().Px() + goodJets[2].tlv().Px()\
      #+ goodJets[3].tlv().Px() + leadlepton.tlv().Px())
      #TransMomNeutrino_y = - (goodJets[0].tlv().Py() + goodJets[1].tlv().Py() + goodJets[2].tlv().Py()\
      # + goodJets[3].tlv().Py() + leadlepton.tlv().Py())

      #z-Komponente des Neutrino-Impulses
      def NeutrinoZMomentum(lepton_vec, neutrino_trans_momentum_x, neutrino_trans_momentum_y, neutrino_phi):
            phi = neutrino_phi
            neutrino_px = neutrino_trans_momentum_x
            neutrino_py = neutrino_trans_momentum_y
            lepton = lepton_vec
            #Annahme für die die Berechnung des pz vum Neutrino gilt 
            mw = 80.42
            mu = mw**2/2 + math.cos(phi) * (lepton.Px()*neutrino_px + lepton.Py()*neutrino_py)
            neutrino_trans_squ = neutrino_px**2 + neutrino_py**2
            zMomentPos = (lepton.Pz()*mu + lepton.E()*((lepton.E()**2+lepton.Pz()**2)*neutrino_trans_squ\
                  +mu**2)**(1/2)) / (lepton.E()**2-lepton.Pz()**2)
            zMomentNeg = (lepton.Pz()*mu - lepton.E()*((lepton.E()**2+lepton.Pz()**2)*neutrino_trans_squ\
                  +mu**2)**(1/2)) / (lepton.E()**2-lepton.Pz()**2)
            #Filtern, welches die 'richtige' Lösung ist
            if type(zMomentPos) == complex:
                  if type(zMomentNeg) == complex:
                        if zMomentNeg.real > zMomentPos.real:
                              zMoment = zMomentPos.real
                        else:
                              zMoment = zMomentNeg.real
                  else:
                        zMoment = zMomentNeg
            if type(zMomentNeg) == complex:
                  zMoment = zMomentPos
            else:
                  if zMomentNeg > zMomentPos:
                        zMoment = zMomentPos
                  else:
                        zMoment = zMomentNeg
            #print(zMoment)
            return zMoment


      #print('0: ' + str(leadlepton.tlv()[0]) + '\n1: ' + str(leadlepton.tlv()[1]) + '\n\n2: ' \
      #      + str(leadlepton.tlv()[2]) + '\n3: ' + str(leadlepton.tlv()[3]))
      #print('tlv x: ' + str(etmiss.tlv().Px()))
      #print('Mein x: ' + str(TransMomNeutrino_x))
      #print('tlv y: ' + str(etmiss.tlv().Py()))
      #print('Mein y:' + str(TransMomNeutrino_y))
      #print('#########################')
      etmiss.tlv().SetPz(NeutrinoZMomentum(leadlepton.tlv(), etmiss.tlv().Px(), etmiss.tlv().Py(), etmiss.phi())) 
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
      # Definierung von chi² Funktion
      def chi_squared(hadronic_w, semilep_top, hadronic_top):
            #values
            m_w_jjexp = 79.3 #GeV
            sigma_mwjjexp = 18.0 #GeV
            sigma_mtexp = 32.1 #GeV
            #Calculation of chi squared
            chi_squ = (hadronic_w.M() - m_w_jjexp)**2 / (sigma_mwjjexp**2) +\
             (semilep_top.M() - hadronic_top.M())**2 / (sigma_mtexp**2)
            return chi_squ

      #Überprüfung, welches der jets als bjet getaggt wurde
      btagged_jet = goodJets[b_jet_tag[0]].tlv()
      first_other_jet = goodJets[other_jets[0]].tlv()
      second_other_jet = goodJets[other_jets[1]].tlv()
      third_other_jet = goodJets[other_jets[2]].tlv()
      #Berechnung der chi² der verschiedenen Permutationen
      #Listen um alle Permutationen abzuspeichern und später darauf zuzugreifen
      #6 Permutationen und 10 teilnehmende Teilchen
      #Aufschlüsselung Teilchen: 
      #[0]=semilep_top, [1]=lepton, [2]=neutrino, [3]=lep_w, [4]=semilep_b, [5]=hadr_top, [6]=qjet, [7]=qbarjet
      #[8]=hadr_w, [9]=hadr_b
      for i in range(0, 6):
            #semileptonischer t-Zerfall
            poss13_semilep_t = Wleptonic + btagged_jet
            if i == 0:
                  poss1_hadr_w = first_other_jet + second_other_jet
                  poss1_hadr_t = poss1_hadr_w + third_other_jet
                  poss1_chi_squ = chi_squared(poss1_hadr_w, poss13_semilep_t, poss1_hadr_t)
                  poss1_list = [poss13_semilep_t, leadlepton.tlv(), etmiss.tlv(),\
                   Wleptonic, btagged_jet, poss1_hadr_t, first_other_jet, second_other_jet,\
                   poss1_hadr_w, third_other_jet]
            if i == 1:
                  poss2_hadr_w = second_other_jet + third_other_jet
                  poss2_hadr_t = poss2_hadr_w + first_other_jet
                  poss2_chi_squ = chi_squared(poss2_hadr_w, poss13_semilep_t, poss2_hadr_t)
                  poss2_list = [poss13_semilep_t, leadlepton.tlv(), etmiss.tlv(),\
                   Wleptonic, btagged_jet, poss2_hadr_t, second_other_jet, third_other_jet,\
                   poss2_hadr_w, first_other_jet]
            if i == 2:
                  poss3_hadr_w = third_other_jet + first_other_jet
                  poss3_hadr_t = poss3_hadr_w + second_other_jet
                  poss3_chi_squ = chi_squared(poss3_hadr_w, poss13_semilep_t, poss3_hadr_t)
                  poss3_list = [poss13_semilep_t, leadlepton.tlv(), etmiss.tlv(),\
                   Wleptonic, btagged_jet, poss3_hadr_t, third_other_jet, first_other_jet,\
                   poss3_hadr_w, second_other_jet]
            #hadronischer t-Zerfall
            if i == 3:
                  poss4_semilep_t = Wleptonic + third_other_jet
                  poss4_hadr_w = first_other_jet + second_other_jet
                  poss4_hadr_t = poss4_hadr_w + btagged_jet
                  poss4_chi_squ = chi_squared(poss4_hadr_w, poss4_semilep_t, poss4_hadr_t)
                  poss4_list = [poss4_semilep_t, leadlepton.tlv(), etmiss.tlv(),\
                   Wleptonic, third_other_jet, poss4_hadr_t, first_other_jet, second_other_jet,\
                   poss4_hadr_w, btagged_jet]
            if i == 4:
                  poss5_semilep_t = Wleptonic + first_other_jet
                  poss5_hadr_w = second_other_jet + third_other_jet
                  poss5_hadr_t = poss5_hadr_w + btagged_jet
                  poss5_chi_squ = chi_squared(poss5_hadr_w, poss5_semilep_t, poss5_hadr_t)
                  poss5_list = [poss5_semilep_t, leadlepton.tlv(), etmiss.tlv(),\
                   Wleptonic, first_other_jet, poss5_hadr_t, second_other_jet, third_other_jet,\
                   poss5_hadr_w, btagged_jet]
            if i == 5:
                  poss6_semilep_t = Wleptonic + second_other_jet
                  poss6_hadr_w = third_other_jet + first_other_jet
                  poss6_hadr_t = poss6_hadr_w + btagged_jet
                  poss6_chi_squ = chi_squared(poss6_hadr_w, poss6_semilep_t, poss6_hadr_t)
                  poss6_list = [poss6_semilep_t, leadlepton.tlv(), etmiss.tlv(),\
                   Wleptonic, second_other_jet, poss6_hadr_t, third_other_jet, first_other_jet,\
                   poss6_hadr_w, btagged_jet]

      #print('Verschiedene chi² der unterschiedlichen Permutationen des semileptonischen top-Zerfalls: \n 1: '\
      # + str(poss1_chi_squ) + '\n 2: ' + str(poss2_chi_squ) + '\n 3: ' + str(poss3_chi_squ))   

      #Auswahl des chi² mit dem kleinsten Wert, d.h. Auswahl der richtigen Zuordnung der bJets
      all_chi_squ = [poss1_chi_squ, poss2_chi_squ, poss3_chi_squ, poss4_chi_squ, poss5_chi_squ, poss6_chi_squ]
      min_chi_squ = min(all_chi_squ)
      min_index = all_chi_squ.index(min_chi_squ)
      

      # Leptonic W boson
      Wmass = Wleptonic.M() # vorher = 80

      # Hadronic top quark
      topmass = 172.5 # vorher = 172.5


      ##########################################################################

      ###########################################################################
      ### Ab hier findet das fuellen der Histogramme statt


      # Histograms detailing event information
      self.hist_vxp_z.Fill(eventinfo.primaryVertexPosition(), weight)
      self.hist_pvxp_n.Fill(eventinfo.numberOfVertices(), weight)

      # histograms for the W boson properties
      self.hist_WtMass.Fill(AnalysisHelpers.WTransverseMass(leadlepton, etmiss), weight)

      # histograms for missing et
      self.hist_etmiss.Fill(etmiss.et(),weight)

      # histograms detailing lepton information
      self.hist_leptn.Fill(len(goodLeptons), weight)
      self.hist_leptpt.Fill(leadlepton.pt(), weight)
      self.hist_lepteta.Fill(leadlepton.eta(), weight)
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

      self.hist_Wmass.Fill(Wmass, weight)
      self.hist_topmass.Fill(topmass, weight)

      return True

  def finalize(self):
      pass
