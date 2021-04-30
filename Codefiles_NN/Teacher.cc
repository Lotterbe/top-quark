#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "NeuroBayesTeacher.hh"

#include <string>
#include <algorithm>
#include <vector>
#include <iostream>



// The initialisation of this neural network was not optimised.
// It is suggested that the user, by looking at the text output and at the analysis plots,
// tries to optimise the setup (suggestion: parameters that should be optimised first are 
// learning speed and number of iteration, rather than number of intermediate nodes or
// preprocessing flags).

Double_t GetValue(TChain *t, std::string name)
{
// Function which returns value of asked variable. We assume, that we have
// only Branch.Leaf variables

  int dotPosition=name.find(".");

  std::string branchName=name.substr(0,dotPosition);
  std::string leafName=name.substr(dotPosition+1,name.size()-dotPosition);

//  std::cout<<"Branch: "<<branchName<<", Leaf: "<<leafName<<std::endl;

  TBranch *branch=t->GetBranch(branchName.data());
  TLeaf *leaf=branch->GetLeaf(leafName.data());
  
  return leaf->GetValue();

}


void Teacher()
{

  //create NeuroBayes instance
  NeuroBayesTeacher* nb = NeuroBayesTeacher::Instance();

  //setup network topology

  // =====  TODO Adjust nvar accordingly ===============
  const int  nvar = 15;     //number of input variables
  // ===================================================

  // number of nodes in input layer is set automatically
  nb->NB_DEF_NODE1(nvar+1);       // nodes in input layer 
  nb->NB_DEF_NODE2(nvar);      // nodes in hidden layer 
  nb->NB_DEF_NODE3(1);       // nodes in output layer
  
  nb->NB_DEF_TASK("CLA");    // binominal classification
  
  nb->NB_DEF_PRE(112);             // 612, 112
  nb->NB_DEF_REG("REG");           // 'OFF','REG' (def) ,'ARD','ASR','ALL'
  nb->NB_DEF_LOSS("ENTROPY");      // 'ENTROPY'(def),'QUADRATIC'
  nb->NB_DEF_SHAPE("OFF");        // 'OFF', 'INCL', 'TOTL'

  nb->NB_DEF_RTRAIN(0.8);          // use 70% of events for training
  nb->NB_DEF_EPOCH(200);           // weight update after n events
  
  nb->NB_DEF_SPEED(2.0);           // multiplicative factor to enhance global learning speed
  nb->NB_DEF_MAXLEARN(1.0);        // multiplicative factor to limit the global learning speed in any direction, this number should be smaller than NB_DEF_SPEED

  nb->NB_DEF_ITER(100);             // number of training iteration
  nb->NB_DEF_METHOD("BFGS");
  nb->NB_DEF_DEBUG(0);

  int i= 4701;
  int j=21; 
  nb->NB_RANVIN(i,j,2);            // random number seed initialisation, i has to be an odd number, the third argument is a debugging flag

  //
  //individual preprocessing flags
  //

  for (int iv=0; iv < nvar; iv++)
    nb->SetIndividualPreproFlag(iv,14);


  nb->SetOutputFile("ttbar_expert.nb");  // expert file

  //
  // Open Input Ntuple
  //
  //  signal part
  //  Einfuegen der Datensaetze: SIGNAL
  TChain *signal=new TChain("tuple");
  signal->Add("/home/student/atlas-outreach-PyROOT-framework-13tev/results/ttbar_lep.root");

  //  background upper sideband
  //  Einfuegen der Datensaetze: UNTERGRUND
  TChain *back=new TChain("tuple");
  back->Add("/home/student/atlas-outreach-PyROOT-framework-13tev/results/WenuWithB.root");
  //back->Add("WenuJetsBVeto.root");
  //back->Add("WenuNoJetsBVeto.root");
  //back->Add("WmunuJetsBVeto.root");
  back->Add("/home/student/atlas-outreach-PyROOT-framework-13tev/results/WmunuWithB.root");
  //back->Add("WmunuNoJetsBVeto.root");
  //back->Add("WtaunuJetsBVeto.root");
  back->Add("/home/student/atlas-outreach-PyROOT-framework-13tev/results/WtaunuWithB.root");
  //back->Add("WtaunuNoJetsBVeto.root");

  Int_t NEntriesbg = Int_t(back->GetEntries());

  
  int eventbg=0;

  do 
    {
      Int_t ientry = back->GetEntry(eventbg);
      if(ientry<=0) break;

      float InputArray[nvar];   

      memset(InputArray,-999,nvar*sizeof(float)); //reset
            
      nb->SetWeight(1.0);  //set weight of event

      // set Target
      nb->SetTarget(0.0) ; // event is a BACKGROUND event

      // Vorher
      // InputArray[0] = GetValue(back,"W_mass");
      // InputArray[1] = GetValue(back,"W_pt");

      // ==== TODO add more input variables

      InputArray[0] = GetValue(back,"HadronicTopMass");
      InputArray[1] = GetValue(back,"SemilepTopMass");
      InputArray[2] = GetValue(back,"SemilepWMass");
      InputArray[3] = GetValue(back,"HadronicTopEta");
      InputArray[4] = GetValue(back,"SemilepTopEta");
      InputArray[5] = GetValue(back,"SemilepWEta");
      InputArray[6] = GetValue(back,"LepEta");
      InputArray[7] = GetValue(back,"HadronicTopTransMom");
      InputArray[8] = GetValue(back,"SemilepTopTransMom");
      InputArray[9] = GetValue(back,"SemilepWTransMom");
      InputArray[10] = GetValue(back,"LepTransMom");
      InputArray[11] = GetValue(back,"TotalTransMom");
      InputArray[12] = GetValue(back,"COMOtherJets");
      InputArray[13] = GetValue(back,"COMbJets");
      InputArray[14] = GetValue(back,"COMTotal");



      //===================================

      //pass input to NeuroBayses
      nb->SetNextInput(nvar,InputArray);
   
      // increase counter for total number of events read in
      eventbg++; 
    } while (eventbg <= NEntriesbg); //do  **** check for NB_MAXPATTERN ****


  Int_t NEntries = Int_t(signal->GetEntries());
 
  int event=0;

  do 
    {
      Int_t ientry = signal->GetEntry(event);
      if(ientry<=0) break;

      float InputArray[nvar];   

      memset(InputArray,-999,nvar*sizeof(float)); //reset
      // get min variables
      
      nb->SetWeight((float) NEntriesbg/ (float) NEntries);  //set weight of event

      // set Target
      nb->SetTarget(1.0); // event is a SIGNAL event

      // vorher
      // InputArray[0] = GetValue(signal,"W_mass");
      // InputArray[1] = GetValue(signal,"W_pt");

      // ==== TODO add more input variables

 InputArray[0] = GetValue(signal,"HadronicTopMass");
      InputArray[1] = GetValue(signal,"SemilepTopMass");
      InputArray[2] = GetValue(signal,"SemilepWMass");
      InputArray[3] = GetValue(signal,"HadronicTopEta");
      InputArray[4] = GetValue(signal,"SemilepTopEta");
      InputArray[5] = GetValue(signal,"SemilepWEta");
      InputArray[6] = GetValue(signal,"LepEta");
      InputArray[7] = GetValue(signal,"HadronicTopTransMom");
      InputArray[8] = GetValue(signal,"SemilepTopTransMom");
      InputArray[9] = GetValue(signal,"SemilepWTransMom");
      InputArray[10] = GetValue(signal,"LepTransMom");
      InputArray[11] = GetValue(signal,"TotalTransMom");
      InputArray[12] = GetValue(signal,"COMOtherJets");
      InputArray[13] = GetValue(signal,"COMbJets");
      InputArray[14] = GetValue(signal,"COMTotal");

      //===================================
	  
      //pass input to NeuroBayses
      nb->SetNextInput(nvar,InputArray);

      // increase counter for total number of events read in
      event++; 
    } while (event <= NEntries); //do  **** check for NB_MAXPATTERN ****


  //perform training
  nb->TrainNet();

  return;

} 

int main(int argc, char** argv)
{
  Teacher();
  return 0;
}

