#ifndef HISTOGRAMMANAGER_h
#define HISTOGRAMMANAGER_h

#include <TH2.h>
#include <iostream>
#include <TFile.h>
#include <TTree.h>

class HistogramManager{

 public:
  HistogramManager(std::string FILE);
  ~HistogramManager();

  TFile *file;

  //TTree
  TTree* m_tree;

  Int_t eventnumber;
  Int_t bcid;
  //Offline
  Int_t muon_n;
  std::vector<float>   *muon_pt;
  std::vector<float>   *muon_eta;
  std::vector<float>   *muon_phi;
  std::vector<float>   *muon_m;
  std::vector<int>     *muon_charge;
  std::vector<int>     *muon_author;
  std::vector<int>     *muon_Type;
  //Museg
  Int_t           museg_n;
  std::vector<float>   *museg_x;
  std::vector<float>   *museg_y;
  std::vector<float>   *museg_z;
  std::vector<float>   *museg_px;
  std::vector<float>   *museg_py;
  std::vector<float>   *museg_pz;
  std::vector<float>   *museg_t0;
  std::vector<float>   *museg_t0error;
  std::vector<float>   *museg_chi2;
  std::vector<float>   *museg_ndof;
  std::vector<int>     *museg_sector;
  std::vector<int>     *museg_stationName;
  std::vector<int>     *museg_stationEta;
  std::vector<int>     *museg_author;
  //Extrapolate
  Int_t           ext_mu_n;
  std::vector<int>     *ext_mu_type;
  std::vector<int>     *ext_mu_index;
  std::vector<int>     *ext_mu_size;
  std::vector<std::vector<int> > *ext_mu_targetVec;
  std::vector<std::vector<float> > *ext_mu_targetDistanceVec;
  std::vector<std::vector<float> > *ext_mu_targetEtaVec;
  std::vector<std::vector<float> > *ext_mu_targetPhiVec;
  std::vector<std::vector<float> > *ext_mu_targetDeltaEtaVec;
  std::vector<std::vector<float> > *ext_mu_targetDeltaPhiVec;
  std::vector<std::vector<float> > *ext_mu_targetPxVec;
  std::vector<std::vector<float> > *ext_mu_targetPyVec;
  std::vector<std::vector<float> > *ext_mu_targetPzVec;
  //TGC Coin Data
  Int_t           tgc_coin_n;
  std::vector<float>   *tgc_coin_x_In;
  std::vector<float>   *tgc_coin_y_In;
  std::vector<float>   *tgc_coin_z_In;
  std::vector<float>   *tgc_coin_x_Out;
  std::vector<float>   *tgc_coin_y_Out;
  std::vector<float>   *tgc_coin_z_Out;
  std::vector<float>   *tgc_coin_width_In;
  std::vector<float>   *tgc_coin_width_Out;
  std::vector<float>   *tgc_coin_width_R;
  std::vector<float>   *tgc_coin_width_Phi;
  std::vector<int>     *tgc_coin_isAside;
  std::vector<int>     *tgc_coin_isForward;
  std::vector<int>     *tgc_coin_isStrip;
  std::vector<int>     *tgc_coin_isInner;
  std::vector<int>     *tgc_coin_type;
  std::vector<int>     *tgc_coin_trackletId;
  std::vector<int>     *tgc_coin_trackletIdStrip;
  std::vector<int>     *tgc_coin_phi;
  std::vector<int>     *tgc_coin_roi;
  std::vector<int>     *tgc_coin_pt;
  std::vector<int>     *tgc_coin_delta;
  std::vector<int>     *tgc_coin_sub;
  std::vector<int>     *tgc_coin_veto;
  std::vector<int>     *tgc_coin_bunch;
  std::vector<int>     *tgc_coin_inner;
  //TGC_RPD
  Int_t           TGC_prd_n;
  std::vector<float>   *TGC_prd_x;
  std::vector<float>   *TGC_prd_y;
  std::vector<float>   *TGC_prd_z;
  std::vector<float>   *TGC_prd_shortWidth;
  std::vector<float>   *TGC_prd_longWidth;
  std::vector<float>   *TGC_prd_length;
  std::vector<int>     *TGC_prd_isStrip;
  std::vector<int>     *TGC_prd_gasGap;
  std::vector<int>     *TGC_prd_channel;
  std::vector<int>     *TGC_prd_eta;
  std::vector<int>     *TGC_prd_phi;
  std::vector<int>     *TGC_prd_station;
  std::vector<int>     *TGC_prd_bunch;
  //RPC_PRD
  Int_t           RPC_prd_n;
  std::vector<float>   *RPC_prd_x;
  std::vector<float>   *RPC_prd_y;
  std::vector<float>   *RPC_prd_z;
  std::vector<float>   *RPC_prd_x2;
  std::vector<float>   *RPC_prd_y2;
  std::vector<float>   *RPC_prd_z2;
  std::vector<float>   *RPC_prd_time;
  std::vector<int>     *RPC_prd_triggerInfo;
  std::vector<int>     *RPC_prd_ambiguityFlag;
  std::vector<int>     *RPC_prd_measuresPhi;
  std::vector<int>     *RPC_prd_inRibs;
  std::vector<int>     *RPC_prd_station;
  std::vector<int>     *RPC_prd_stationEta;
  std::vector<int>     *RPC_prd_stationPhi;
  std::vector<int>     *RPC_prd_doubletR;
  std::vector<int>     *RPC_prd_doubletZ;
  std::vector<double>  *RPC_prd_stripWidth;
  std::vector<double>  *RPC_prd_stripLength;
  std::vector<int>     *RPC_prd_gasGap;
  std::vector<int>     *RPC_prd_channel;
  //RoI
  Int_t           muctpi_ndatawords;
  std::vector<float>   *muctpi_eta;
  std::vector<float>   *muctpi_phi;
  std::vector<short>   *muctpi_source;
  std::vector<short>   *muctpi_hemisphere;
  std::vector<short>   *muctpi_bcid;
  std::vector<short>   *muctpi_sectorID;
  std::vector<short>   *muctpi_thrNumber;
  std::vector<short>   *muctpi_roi;
  std::vector<short>   *muctpi_veto;
  std::vector<short>   *muctpi_charge;
  std::vector<short>   *muctpi_candidateVetoed;
  //HLT
  Int_t           HLT_info_n;
  std::vector<std::string>  *HLT_info_chain;
  std::vector<int>     *HLT_info_isPassed;
  std::vector<std::vector<int> > *HLT_info_typeVec;
  std::vector<std::vector<float> > *HLT_info_ptVec;
  std::vector<std::vector<float> > *HLT_info_etaVec;
  std::vector<std::vector<float> > *HLT_info_phiVec;
  //Tile
  Int_t           TILE_murcv_trig_n;
  std::vector<int>     *TILE_murcv_trig_mod;
  std::vector<int>     *TILE_murcv_trig_part;
  std::vector<bool>    *TILE_murcv_trig_bit0;
  std::vector<bool>    *TILE_murcv_trig_bit1;
  std::vector<bool>    *TILE_murcv_trig_bit2;
  std::vector<bool>    *TILE_murcv_trig_bit3;
  Int_t           TILE_murcv_raw_n;
  std::vector<float>   *TILE_murcv_raw_count;
  std::vector<float>   *TILE_murcv_raw_energy;
  std::vector<int>     *TILE_murcv_raw_ros;
  std::vector<int>     *TILE_murcv_raw_drawer;
  std::vector<int>     *TILE_murcv_raw_channel;
  Int_t           TILE_murcv_digit_n;
  std::vector<int>     *TILE_murcv_digit_nSamples;
  std::vector<int>     *TILE_murcv_digit_ros;
  std::vector<int>     *TILE_murcv_digit_drawer;
  std::vector<int>     *TILE_murcv_digit_channel;
  std::vector<std::vector<float> > *TILE_murcv_digit_sampleVec;
  //Truth
  uint32_t             mc_n;
  std::vector<float>   *mc_pt;
  std::vector<float>   *mc_eta;
  std::vector<float>   *mc_phi;
  std::vector<float>   *mc_m;
  std::vector<int>     *mc_charge;
  //Run-3
  Int_t                TGC_Run3_n;
  std::vector<int>     *TGC_Run3_pt;
  std::vector<int>     *TGC_Run3_type;
  std::vector<int>     *TGC_Run3_station;
  std::vector<int>     *TGC_Run3_DR;
  std::vector<int>     *TGC_Run3_DPhi;
  std::vector<int>     *TGC_Run3_TypeDPhi;
  std::vector<int>     *TGC_Run3_TypeDR;
  std::vector<int>     *TGC_Run3_Side;
  std::vector<int>     *TGC_Run3_RoI;
  std::vector<int>     *TGC_Run3_PhiSector;
  std::vector<bool>     *TGC_Run3_IsEndcap;
  std::vector<int>     *TGC_Run3_TrackletIdWire;
  std::vector<int>     *TGC_Run3_TrackletIdStrip;
  std::vector<float>   *TGC_Run3_x;
  std::vector<float>   *TGC_Run3_y;
  std::vector<float>   *TGC_Run3_z;
  std::vector<float>   *TGC_Run3_R;
  std::vector<float>   *TGC_Run3_Phi;
  std::vector<int>     *TGC_Run3_Charge;

  void Clear();

 private:

  std::string m_FileName;

};

#endif
