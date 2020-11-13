#include "HistogramManager.h"

#include <TH2.h>
#include <TFile.h>
#include <iostream>
#include <iomanip> 
#include <fstream>
#include <TTree.h>

HistogramManager::HistogramManager(std::string FILE){

  m_FileName=FILE;
  std::string FILE_NAME=m_FileName+"_Histgram.root";
  std::cout<<"Create New File="<<FILE_NAME<<std::endl;
  file= new TFile(FILE_NAME.c_str(),"recreate");

  //TTree
  m_tree = new TTree("physics","physics");
  //Event Information
  m_tree->Branch("eventnumber",&eventnumber,"eventnumber/I");
  m_tree->Branch("bcid",&bcid,"bcid/I");
  //Offline
  m_tree->Branch("muon_n",&muon_n,"muon_n/I");
  muon_pt=new std::vector<float>; m_tree->Branch("muon_pt",&muon_pt);
  muon_eta=new std::vector<float>; m_tree->Branch("muon_eta",&muon_eta);
  muon_phi=new std::vector<float>; m_tree->Branch("muon_phi",&muon_phi);
  muon_m=new std::vector<float>; m_tree->Branch("muon_m",&muon_m);
  muon_charge=new std::vector<int>; m_tree->Branch("muon_charge",&muon_charge);
  muon_author=new std::vector<int>; m_tree->Branch("muon_author",&muon_author);
  muon_Type=new std::vector<int>; m_tree->Branch("muon_Type",&muon_Type);
  //Extrapolate
  m_tree->Branch("ext_mu_n",&ext_mu_n,"ext_mu_n/I");
  ext_mu_type=new std::vector<int>; m_tree->Branch("ext_mu_type",&ext_mu_type);
  ext_mu_index=new std::vector<int>; m_tree->Branch("ext_mu_index",&ext_mu_index);
  ext_mu_size=new std::vector<int>; m_tree->Branch("ext_mu_size",&ext_mu_size);
  ext_mu_targetVec=new std::vector<std::vector<int>>; m_tree->Branch("ext_mu_targetVec",&ext_mu_targetVec);
  ext_mu_targetDistanceVec=new std::vector<std::vector<float>>; m_tree->Branch("ext_mu_targetDistanceVec",&ext_mu_targetDistanceVec);
  ext_mu_targetEtaVec=new std::vector<std::vector<float>>; m_tree->Branch("ext_mu_targetEtaVec",&ext_mu_targetEtaVec);
  ext_mu_targetPhiVec=new std::vector<std::vector<float>>; m_tree->Branch("ext_mu_targetPhiVec",&ext_mu_targetPhiVec);
  ext_mu_targetDeltaEtaVec=new std::vector<std::vector<float>>; m_tree->Branch("ext_mu_targetDeltaEtaVec",&ext_mu_targetDeltaEtaVec);
  ext_mu_targetDeltaPhiVec=new std::vector<std::vector<float>>; m_tree->Branch("ext_mu_targetDeltaPhiVec",&ext_mu_targetDeltaPhiVec);
  ext_mu_targetPxVec=new std::vector<std::vector<float>>; m_tree->Branch("ext_mu_targetPxVec",&ext_mu_targetPxVec);
  ext_mu_targetPyVec=new std::vector<std::vector<float>>; m_tree->Branch("ext_mu_targetPyVec",&ext_mu_targetPyVec);
  ext_mu_targetPzVec=new std::vector<std::vector<float>>; m_tree->Branch("ext_mu_targetPzVec",&ext_mu_targetPzVec);
  //TGC Coin Data
  m_tree->Branch("tgc_coin_n",&tgc_coin_n,"tgc_coin_n/I");
  tgc_coin_x_In=new std::vector<float>; m_tree->Branch("tgc_coin_x_In",&tgc_coin_x_In);
  tgc_coin_y_In=new std::vector<float>; m_tree->Branch("tgc_coin_y_In",&tgc_coin_y_In);
  tgc_coin_z_In=new std::vector<float>; m_tree->Branch("tgc_coin_z_In",&tgc_coin_z_In);
  tgc_coin_x_Out=new std::vector<float>; m_tree->Branch("tgc_coin_x_Out",&tgc_coin_x_Out);
  tgc_coin_y_Out=new std::vector<float>; m_tree->Branch("tgc_coin_y_Out",&tgc_coin_y_Out);
  tgc_coin_z_Out=new std::vector<float>; m_tree->Branch("tgc_coin_z_Out",&tgc_coin_z_Out);
  tgc_coin_width_In=new std::vector<float>; m_tree->Branch("tgc_coin_width_In",&tgc_coin_width_In);
  tgc_coin_width_Out=new std::vector<float>; m_tree->Branch("tgc_coin_width_Out",&tgc_coin_width_Out);
  tgc_coin_width_R=new std::vector<float>; m_tree->Branch("tgc_coin_width_R",&tgc_coin_width_R);
  tgc_coin_width_Phi=new std::vector<float>; m_tree->Branch("tgc_coin_width_Phi",&tgc_coin_width_Phi);
  tgc_coin_isAside=new std::vector<int>; m_tree->Branch("tgc_coin_isAside",&tgc_coin_isAside);
  tgc_coin_isForward=new std::vector<int>; m_tree->Branch("tgc_coin_isForward",&tgc_coin_isForward);
  tgc_coin_isStrip=new std::vector<int>; m_tree->Branch("tgc_coin_isStrip",&tgc_coin_isStrip);
  tgc_coin_isInner=new std::vector<int>; m_tree->Branch("tgc_coin_isInner",&tgc_coin_isInner);
  tgc_coin_type=new std::vector<int>; m_tree->Branch("tgc_coin_type",&tgc_coin_type);
  tgc_coin_trackletId=new std::vector<int>; m_tree->Branch("tgc_coin_trackletId",&tgc_coin_trackletId);
  tgc_coin_trackletIdStrip=new std::vector<int>; m_tree->Branch("tgc_coin_trackletIdStrip",&tgc_coin_trackletIdStrip);
  tgc_coin_phi=new std::vector<int>; m_tree->Branch("tgc_coin_phi",&tgc_coin_phi);
  tgc_coin_roi=new std::vector<int>; m_tree->Branch("tgc_coin_roi",&tgc_coin_roi);
  tgc_coin_pt=new std::vector<int>; m_tree->Branch("tgc_coin_pt",&tgc_coin_pt);
  tgc_coin_delta=new std::vector<int>; m_tree->Branch("tgc_coin_delta",&tgc_coin_delta);
  tgc_coin_sub=new std::vector<int>; m_tree->Branch("tgc_coin_sub",&tgc_coin_sub);
  tgc_coin_veto=new std::vector<int>; m_tree->Branch("tgc_coin_veto",&tgc_coin_veto);
  tgc_coin_bunch=new std::vector<int>; m_tree->Branch("tgc_coin_bunch",&tgc_coin_bunch);
  tgc_coin_inner=new std::vector<int>; m_tree->Branch("tgc_coin_inner",&tgc_coin_inner);
  //RoI
  m_tree->Branch("muctpi_ndatawords",&muctpi_ndatawords,"muctpi_ndatawords/I");
  muctpi_eta=new std::vector<float>; m_tree->Branch("muctpi_eta",&muctpi_eta);
  muctpi_phi=new std::vector<float>; m_tree->Branch("muctpi_phi",&muctpi_phi);
  muctpi_source=new std::vector<short>; m_tree->Branch("muctpi_source",&muctpi_source);
  muctpi_hemisphere=new std::vector<short>; m_tree->Branch("muctpi_hemisphere",&muctpi_hemisphere);
  muctpi_bcid=new std::vector<short>; m_tree->Branch("muctpi_bcid",&muctpi_bcid);
  muctpi_sectorID=new std::vector<short>; m_tree->Branch("muctpi_sectorID",&muctpi_sectorID);
  muctpi_thrNumber=new std::vector<short>; m_tree->Branch("muctpi_thrNumber",&muctpi_thrNumber);
  muctpi_roi=new std::vector<short>; m_tree->Branch("muctpi_roi",&muctpi_roi);
  muctpi_veto=new std::vector<short>; m_tree->Branch("muctpi_veto",&muctpi_veto);
  muctpi_charge=new std::vector<short>; m_tree->Branch("muctpi_charge",&muctpi_charge);
  muctpi_candidateVetoed=new std::vector<short>; m_tree->Branch("muctpi_candidateVetoed",&muctpi_candidateVetoed);
  //HLT
  m_tree->Branch("HLT_info_n",&HLT_info_n,"HLT_info_n/I");
  HLT_info_chain=new std::vector<std::string>; m_tree->Branch("HLT_info_chain",&HLT_info_chain);
  HLT_info_isPassed=new std::vector<int>; m_tree->Branch("HLT_info_isPassed",&HLT_info_isPassed);
  HLT_info_typeVec=new std::vector<std::vector<int>>; m_tree->Branch("HLT_info_typeVec",&HLT_info_typeVec);
  HLT_info_ptVec=new std::vector<std::vector<float>>; m_tree->Branch("HLT_info_ptVec",&HLT_info_ptVec);
  HLT_info_etaVec=new std::vector<std::vector<float>>; m_tree->Branch("HLT_info_etaVec",&HLT_info_etaVec);
  HLT_info_phiVec=new std::vector<std::vector<float>>; m_tree->Branch("HLT_info_phiVec",&HLT_info_phiVec);
  //Truth
  m_tree->Branch("mc_n",&mc_n,"mc_n/I");
  mc_pt=new std::vector<float>; m_tree->Branch("mc_pt",&mc_pt);
  mc_eta=new std::vector<float>; m_tree->Branch("mc_eta",&mc_eta);
  mc_phi=new std::vector<float>; m_tree->Branch("mc_phi",&mc_phi);
  mc_m=new std::vector<float>; m_tree->Branch("mc_m",&mc_m);
  mc_charge=new std::vector<int>; m_tree->Branch("mc_charge",&mc_charge);
  //Run-3
  m_tree->Branch("TGC_Run3_n",&TGC_Run3_n,"TGC_Run3_n/I");
  TGC_Run3_pt=new std::vector<int>; m_tree->Branch("TGC_Run3_pt",&TGC_Run3_pt);
  TGC_Run3_type=new std::vector<int>; m_tree->Branch("TGC_Run3_type",&TGC_Run3_type);
  TGC_Run3_station=new std::vector<int>; m_tree->Branch("TGC_Run3_station",&TGC_Run3_station);
  TGC_Run3_DR=new std::vector<int>; m_tree->Branch("TGC_Run3_DR",&TGC_Run3_DR);
  TGC_Run3_DPhi=new std::vector<int>; m_tree->Branch("TGC_Run3_DPhi",&TGC_Run3_DPhi);
  TGC_Run3_TypeDPhi=new std::vector<int>; m_tree->Branch("TGC_Run3_TypeDPhi",&TGC_Run3_TypeDPhi);
  TGC_Run3_TypeDR=new std::vector<int>; m_tree->Branch("TGC_Run3_TypeDR",&TGC_Run3_TypeDR);
  TGC_Run3_Side=new std::vector<int>; m_tree->Branch("TGC_Run3_Side",&TGC_Run3_Side);
  TGC_Run3_RoI=new std::vector<int>; m_tree->Branch("TGC_Run3_RoI",&TGC_Run3_RoI);
  TGC_Run3_PhiSector=new std::vector<int>; m_tree->Branch("TGC_Run3_PhiSector",&TGC_Run3_PhiSector);
  TGC_Run3_IsEndcap=new std::vector<bool>; m_tree->Branch("TGC_Run3_IsEndcap",&TGC_Run3_IsEndcap);
  TGC_Run3_TrackletIdWire=new std::vector<int>; m_tree->Branch("TGC_Run3_TrackletIdWire",&TGC_Run3_TrackletIdWire);
  TGC_Run3_TrackletIdStrip=new std::vector<int>; m_tree->Branch("TGC_Run3_TrackletIdStrip",&TGC_Run3_TrackletIdStrip);
  TGC_Run3_x=new std::vector<float>; m_tree->Branch("TGC_Run3_x",&TGC_Run3_x);
  TGC_Run3_y=new std::vector<float>; m_tree->Branch("TGC_Run3_y",&TGC_Run3_y);
  TGC_Run3_z=new std::vector<float>; m_tree->Branch("TGC_Run3_z",&TGC_Run3_z);
  TGC_Run3_R=new std::vector<float>; m_tree->Branch("TGC_Run3_R",&TGC_Run3_R);
  TGC_Run3_Phi=new std::vector<float>; m_tree->Branch("TGC_Run3_Phi",&TGC_Run3_Phi);
  TGC_Run3_Charge=new std::vector<int>; m_tree->Branch("TGC_Run3_Charge",&TGC_Run3_Charge);
}

void HistogramManager::Clear(){
    muon_pt->clear();
    muon_phi->clear();
    muon_eta->clear();
    muon_m->clear();
    muon_author->clear();
    muon_charge->clear();
    muon_Type->clear();
    ext_mu_type->clear();
    ext_mu_index->clear();
    ext_mu_size->clear();
    ext_mu_targetVec->clear();
    ext_mu_targetDistanceVec->clear();
    ext_mu_targetEtaVec->clear();
    ext_mu_targetPhiVec->clear();
    ext_mu_targetDeltaEtaVec->clear();
    ext_mu_targetDeltaPhiVec->clear();
    ext_mu_targetPxVec->clear();
    ext_mu_targetPyVec->clear();
    ext_mu_targetPzVec->clear();
    tgc_coin_x_In->clear();
    tgc_coin_y_In->clear();
    tgc_coin_z_In->clear();
    tgc_coin_x_Out->clear();
    tgc_coin_y_Out->clear();
    tgc_coin_z_Out->clear();
    tgc_coin_width_In->clear();
    tgc_coin_width_Out->clear();
    tgc_coin_width_R->clear();
    tgc_coin_width_Phi->clear();
    tgc_coin_isAside->clear();
    tgc_coin_isStrip->clear();
    tgc_coin_isForward->clear();
    tgc_coin_isInner->clear();
    tgc_coin_type->clear();
    tgc_coin_trackletId->clear();
    tgc_coin_trackletIdStrip->clear();
    tgc_coin_phi->clear();
    tgc_coin_roi->clear();
    tgc_coin_pt->clear();
    tgc_coin_delta->clear();
    tgc_coin_sub->clear();
    tgc_coin_veto->clear();
    tgc_coin_bunch->clear();
    tgc_coin_inner->clear();
    muctpi_eta->clear();
    muctpi_phi->clear();
    muctpi_source->clear();
    muctpi_hemisphere->clear();
    muctpi_bcid->clear();
    muctpi_sectorID->clear();
    muctpi_thrNumber->clear();
    muctpi_roi->clear();
    muctpi_veto->clear();
    muctpi_charge->clear();
    muctpi_candidateVetoed->clear();
    HLT_info_chain->clear();
    HLT_info_isPassed->clear();
    HLT_info_typeVec->clear();
    HLT_info_ptVec->clear();
    HLT_info_etaVec->clear();
    HLT_info_phiVec->clear();
    mc_pt->clear();
    mc_eta->clear();
    mc_phi->clear();
    mc_m->clear();
    mc_charge->clear();
    TGC_Run3_pt->clear();
    TGC_Run3_type->clear();
    TGC_Run3_station->clear();
    TGC_Run3_DR->clear();
    TGC_Run3_DPhi->clear();
    TGC_Run3_TypeDPhi->clear();
    TGC_Run3_TypeDR->clear();
    TGC_Run3_Side->clear();
    TGC_Run3_RoI->clear();
    TGC_Run3_PhiSector->clear();
    TGC_Run3_IsEndcap->clear();
    TGC_Run3_TrackletIdWire->clear();
    TGC_Run3_TrackletIdStrip->clear();
    TGC_Run3_x->clear();
    TGC_Run3_y->clear();
    TGC_Run3_z->clear();
    TGC_Run3_R->clear();
    TGC_Run3_Phi->clear();
    TGC_Run3_Charge->clear();
    std::cout<<"a"<<std::endl;
}


HistogramManager::~HistogramManager(){
    delete muon_pt; muon_pt=0;
    delete muon_eta; muon_eta=0;
    delete muon_phi; muon_phi=0;
    delete muon_m; muon_m=0;
    delete muon_charge; muon_charge=0;
    delete muon_author; muon_author=0;
    delete muon_Type; muon_Type=0;
    delete ext_mu_type; ext_mu_type=0;
    delete ext_mu_index; ext_mu_index=0;
    delete ext_mu_size; ext_mu_size=0;
    delete ext_mu_targetVec; ext_mu_targetVec=0;
    delete ext_mu_targetDistanceVec; ext_mu_targetDistanceVec=0;
    delete ext_mu_targetEtaVec; ext_mu_targetEtaVec=0;
    delete ext_mu_targetPhiVec; ext_mu_targetPhiVec=0;
    delete ext_mu_targetDeltaEtaVec; ext_mu_targetDeltaEtaVec=0;
    delete ext_mu_targetDeltaPhiVec; ext_mu_targetDeltaPhiVec=0;
    delete ext_mu_targetPxVec; ext_mu_targetPxVec=0;
    delete ext_mu_targetPyVec; ext_mu_targetPyVec=0;
    delete ext_mu_targetPzVec; ext_mu_targetPzVec=0;
    delete tgc_coin_x_In; tgc_coin_x_In=0;
    delete tgc_coin_y_In; tgc_coin_y_In=0;
    delete tgc_coin_z_In; tgc_coin_z_In=0;
    delete tgc_coin_x_Out; tgc_coin_x_Out=0;
    delete tgc_coin_y_Out; tgc_coin_y_Out=0;
    delete tgc_coin_z_Out; tgc_coin_z_Out=0;
    delete tgc_coin_width_In; tgc_coin_width_In=0;
    delete tgc_coin_width_Out; tgc_coin_width_Out=0;
    delete tgc_coin_width_R; tgc_coin_width_R=0;
    delete tgc_coin_width_Phi; tgc_coin_width_Phi=0;
    delete tgc_coin_isAside; tgc_coin_isAside=0;
    delete tgc_coin_isForward; tgc_coin_isForward=0;
    delete tgc_coin_isStrip; tgc_coin_isStrip=0;
    delete tgc_coin_isInner; tgc_coin_isInner=0;
    delete tgc_coin_type; tgc_coin_type=0;
    delete tgc_coin_trackletId; tgc_coin_trackletId=0;
    delete tgc_coin_trackletIdStrip; tgc_coin_trackletIdStrip=0;
    delete tgc_coin_phi; tgc_coin_phi=0;
    delete tgc_coin_pt; tgc_coin_pt=0;
    delete tgc_coin_delta; tgc_coin_delta=0;
    delete tgc_coin_sub; tgc_coin_sub=0;
    delete tgc_coin_veto; tgc_coin_veto=0;
    delete tgc_coin_bunch; tgc_coin_bunch=0;
    delete tgc_coin_inner; tgc_coin_inner=0;
    delete muctpi_eta; muctpi_eta=0;
    delete muctpi_phi; muctpi_phi=0;
    delete muctpi_source; muctpi_source=0;
    delete muctpi_hemisphere; muctpi_hemisphere=0;
    delete muctpi_bcid; muctpi_bcid=0;
    delete muctpi_sectorID; muctpi_sectorID=0;
    delete muctpi_thrNumber; muctpi_thrNumber=0;
    delete muctpi_roi; muctpi_roi=0;
    delete muctpi_veto; muctpi_veto=0;
    delete muctpi_charge; muctpi_charge=0;
    delete muctpi_candidateVetoed; muctpi_candidateVetoed=0;   
    delete HLT_info_chain; HLT_info_chain=0;
    delete HLT_info_isPassed; HLT_info_isPassed=0;
    delete HLT_info_typeVec; HLT_info_typeVec=0;
    delete HLT_info_ptVec; HLT_info_ptVec=0;
    delete HLT_info_etaVec; HLT_info_etaVec=0;
    delete HLT_info_phiVec; HLT_info_phiVec=0;
    delete mc_pt; mc_pt=0;
    delete mc_eta; mc_eta=0;
    delete mc_phi; mc_phi=0;
    delete mc_m; mc_m=0;
    delete mc_charge; mc_charge=0;
    delete TGC_Run3_pt; TGC_Run3_pt=0;
    delete TGC_Run3_type; TGC_Run3_type=0;
    delete TGC_Run3_station; TGC_Run3_station=0;
    delete TGC_Run3_DR; TGC_Run3_DR=0;
    delete TGC_Run3_DPhi; TGC_Run3_DPhi=0;
    delete TGC_Run3_TypeDPhi; TGC_Run3_TypeDPhi=0;
    delete TGC_Run3_TypeDR; TGC_Run3_TypeDR=0;
    delete TGC_Run3_Side; TGC_Run3_Side=0;
    delete TGC_Run3_RoI; TGC_Run3_RoI=0;
    delete TGC_Run3_PhiSector; TGC_Run3_PhiSector=0;
    delete TGC_Run3_IsEndcap; TGC_Run3_IsEndcap=0;
    delete TGC_Run3_TrackletIdWire; TGC_Run3_TrackletIdWire=0;
    delete TGC_Run3_TrackletIdStrip; TGC_Run3_TrackletIdStrip=0;
    delete TGC_Run3_x; TGC_Run3_x=0;
    delete TGC_Run3_y; TGC_Run3_y=0;
    delete TGC_Run3_z; TGC_Run3_z=0;
    delete TGC_Run3_R; TGC_Run3_R=0;
    delete TGC_Run3_Phi; TGC_Run3_Phi=0;
    delete TGC_Run3_Charge; TGC_Run3_Charge=0;

  file->Write();
  file->Close();
  std::cout<<"*************** FILE SAVED ****************"<<std::endl;
}
