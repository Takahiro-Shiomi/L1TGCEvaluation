//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Dec  2 17:14:00 2019 by ROOT version 6.08/06
// from TTree TrigT1TGCNtuple/
// found on file: TrigT1TGCAlgBranch.root
//////////////////////////////////////////////////////////

#ifndef TrigT1TGCNtuple_h
#define TrigT1TGCNtuple_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class TrigT1TGCNtuple {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           eventNumber;
   std::vector<double>  *truthMuon_Eta;
   std::vector<double>  *truthMuon_Phi;
   std::vector<double>  *truthMuon_Pt;
   std::vector<double>  *truthMuon_Eta_onTGC;
   std::vector<double>  *truthMuon_Phi_onTGC;
   std::vector<double>  *truthMuon_Pt_onTGC;
   std::vector<double>  *truthMuon_Z_onTGC;
   std::vector<int>     *roi_side;
   std::vector<int>     *roi_octant;
   std::vector<int>     *roi_module;
   std::vector<int>     *roi_number;
   std::vector<bool>    *roi_3stationFlag;
   std::vector<bool>    *roi_GoodMFFlag;
   std::vector<bool>    *roi_InnerCoincidenceFlag;
   std::vector<int>     *roi_DR;
   std::vector<int>     *roi_DPhi;
   std::vector<double>  *roi_eta;
   std::vector<double>  *roi_phi;

   // List of branches
   TBranch        *b_eventNumber;   //!
   TBranch        *b_truthMuon_Eta;   //!
   TBranch        *b_truthMuon_Phi;   //!
   TBranch        *b_truthMuon_Pt;   //!
   TBranch        *b_truthMuon_Eta_onTGC;   //!
   TBranch        *b_truthMuon_Phi_onTGC;   //!
   TBranch        *b_truthMuon_Pt_onTGC;   //!
   TBranch        *b_truthMuon_Z_onTGC;   //!
   TBranch        *b_roi_side;   //!
   TBranch        *b_roi_octant;   //!
   TBranch        *b_roi_module;   //!
   TBranch        *b_roi_number;   //!
   TBranch        *b_roi_3stationFlag;   //!
   TBranch        *b_roi_GoodMFFlag;   //!
   TBranch        *b_roi_InnerCoincidenceFlag;   //!
   TBranch        *b_roi_DR;   //!
   TBranch        *b_roi_DPhi;   //!
   TBranch        *b_roi_eta;   //!
   TBranch        *b_roi_phi;   //!

   TrigT1TGCNtuple(TTree *tree=0);
   virtual ~TrigT1TGCNtuple();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TrigT1TGCNtuple_cxx
TrigT1TGCNtuple::TrigT1TGCNtuple(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("TrigT1TGCAlgBranch.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("TrigT1TGCAlgBranch.root");
      }
      f->GetObject("TrigT1TGCNtuple",tree);

   }
   Init(tree);
}

TrigT1TGCNtuple::~TrigT1TGCNtuple()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TrigT1TGCNtuple::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TrigT1TGCNtuple::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TrigT1TGCNtuple::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   truthMuon_Eta = 0;
   truthMuon_Phi = 0;
   truthMuon_Pt = 0;
   truthMuon_Eta_onTGC = 0;
   truthMuon_Phi_onTGC = 0;
   truthMuon_Pt_onTGC = 0;
   truthMuon_Z_onTGC = 0;
   roi_side = 0;
   roi_octant = 0;
   roi_module = 0;
   roi_number = 0;
   roi_3stationFlag = 0;
   roi_GoodMFFlag = 0;
   roi_InnerCoincidenceFlag = 0;
   roi_DR = 0;
   roi_DPhi = 0;
   roi_eta = 0;
   roi_phi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("truthMuon_Eta", &truthMuon_Eta, &b_truthMuon_Eta);
   fChain->SetBranchAddress("truthMuon_Phi", &truthMuon_Phi, &b_truthMuon_Phi);
   fChain->SetBranchAddress("truthMuon_Pt", &truthMuon_Pt, &b_truthMuon_Pt);
   fChain->SetBranchAddress("truthMuon_Eta_onTGC", &truthMuon_Eta_onTGC, &b_truthMuon_Eta_onTGC);
   fChain->SetBranchAddress("truthMuon_Phi_onTGC", &truthMuon_Phi_onTGC, &b_truthMuon_Phi_onTGC);
   fChain->SetBranchAddress("truthMuon_Pt_onTGC", &truthMuon_Pt_onTGC, &b_truthMuon_Pt_onTGC);
   fChain->SetBranchAddress("truthMuon_Z_onTGC", &truthMuon_Z_onTGC, &b_truthMuon_Z_onTGC);
   fChain->SetBranchAddress("roi_side", &roi_side, &b_roi_side);
   fChain->SetBranchAddress("roi_octant", &roi_octant, &b_roi_octant);
   fChain->SetBranchAddress("roi_module", &roi_module, &b_roi_module);
   fChain->SetBranchAddress("roi_number", &roi_number, &b_roi_number);
   fChain->SetBranchAddress("roi_3stationFlag", &roi_3stationFlag, &b_roi_3stationFlag);
   fChain->SetBranchAddress("roi_GoodMFFlag", &roi_GoodMFFlag, &b_roi_GoodMFFlag);
   fChain->SetBranchAddress("roi_InnerCoincidenceFlag", &roi_InnerCoincidenceFlag, &b_roi_InnerCoincidenceFlag);
   fChain->SetBranchAddress("roi_DR", &roi_DR, &b_roi_DR);
   fChain->SetBranchAddress("roi_DPhi", &roi_DPhi, &b_roi_DPhi);
   fChain->SetBranchAddress("roi_eta", &roi_eta, &b_roi_eta);
   fChain->SetBranchAddress("roi_phi", &roi_phi, &b_roi_phi);
   Notify();
}

Bool_t TrigT1TGCNtuple::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TrigT1TGCNtuple::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TrigT1TGCNtuple::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TrigT1TGCNtuple_cxx
