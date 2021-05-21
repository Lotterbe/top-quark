#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TLeaf.h"
#include "NeuroBayesExpert.hh"
#include "TNtuple.h"
#include "TLorentzVector.h"
#include "TVector3.h"

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



int main(int argc, char** argv)
{ 
  
  if (argc < 2)
  {
    std::cout<<"Bitte Dateinamen als Parameter angeben!"<<std::endl;
    return -1; 
  }

  std::string filename=argv[1];

  std::string outname="output/h_" + filename.substr(filename.rfind("/")+1);

  //create NeuroBayes instance
  Expert* nb = new Expert("../train/ttbar_expert.nb");

  TFile input(filename.c_str());
  TFile output(outname.c_str(), "recreate");

  TChain* chain=(TChain*)input.Get("tuple");

  TH1F* h_nn = new TH1F("nn","nn",20,-1,1); 
  
  int event=0;
  
  Int_t NEntries=chain->GetEntries();

  // =====  TODO Adjust nvar accordingly ===============
  int  nvar = 15;     //number of input variables
  // =================================================== 
  
  do 
    { 
      Int_t ientry = chain->GetEntry(event);
      if(ientry<=0) break;


      float InputArray[nvar];   

      memset(InputArray,-999,nvar*sizeof(float)); //reset

      //InputArray[0] = GetValue(chain,"W_mass");
      //InputArray[1] = GetValue(chain,"W_pt");
      // ==== TODO add more input variables
      InputArray[0] = GetValue(chain,"HadronicTopMass");
      InputArray[1] = GetValue(chain,"SemilepTopMass");
      InputArray[2] = GetValue(chain,"SemilepWMass");
      //InputArray[3] = GetValue(chain,"LepMass");
      InputArray[3] = GetValue(chain,"HadronicTopEta");
      InputArray[4] = GetValue(chain,"SemilepTopEta");
      InputArray[5] = GetValue(chain,"SemilepWEta");
      InputArray[6] = GetValue(chain,"LepEta");
      InputArray[7] = GetValue(chain,"HadronicTopTransMom");
      InputArray[8] = GetValue(chain,"SemilepTopTransMom");
      InputArray[9] = GetValue(chain,"SemilepWTransMom");
      InputArray[10] = GetValue(chain,"LepTransMom");
      InputArray[11] = GetValue(chain,"TotalTransMom");
      InputArray[12] = GetValue(chain,"COMOtherJets");
      InputArray[13] = GetValue(chain,"COMbJets");
      InputArray[14] = GetValue(chain,"COMTotal");

      double nn_out =  nb->nb_expert(InputArray);
 //     double weight = GetValue(chain,"weight");
         double weight = 1;
      h_nn->Fill(nn_out,weight);
 
      // increase counter for total number of events read in
      event++; 
  } while (event <= NEntries); //do  **** check for NB_MAXPATTERN ****



  output.cd();
  h_nn->Write();
  output.Close();

  return 0;
}
