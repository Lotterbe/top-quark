import ROOT
import itertools
import math

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
      if not (btags == 2): return False # >=0 vorher
      self.countEvent("btags", weight)


      # apply a cut on the transverse mass of the W boson decaying to leptons
      if not (AnalysisHelpers.WTransverseMass(leadlepton, etmiss) > 30.0): return False

      ##########################################################################
      ## Aufgabe 2
      # finde heraus welche jets in goodJets die hadronischen Jets sind und addiere diese
      Whadronic_indices = []
      b_indices = []
      #for jet in goodJets:
      for i in range(4):
          if goodJets[i].mv1() > 0.7892:
              b_indices.append(i)
          else:
              Whadronic_indices.append(i)

      Whadronic = goodJets[Whadronic_indices[0]].tlv() + goodJets[Whadronic_indices[1]].tlv()
      b1 = goodJets[b_indices[0]].tlv()
      b2 = goodJets[b_indices[1]].tlv()
      dr1 = Whadronic.DeltaR(b1)
      dr2 = Whadronic.DeltaR(b2)
      if dr1 < dr2:
      	thadronic = Whadronic + b1
      else:
      	thadronic = Whadronic + b2
      ## Hadronic W boson
      Wmass = Whadronic.M()	# vorher = 80

      ## Hadronic top quark
      topmass = thadronic.M() # vorher = 172.5



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
