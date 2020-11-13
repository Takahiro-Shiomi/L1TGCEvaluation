//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Dec  3 17:38:02 2019 by ROOT version 6.08/06
// from TTree physics/physics
// found on file: group.det-muon/group.det-muon.19069326.L1TGCNtuple._000011.root
//////////////////////////////////////////////////////////

#ifndef physics_h
#define physics_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "TH2D.h"
// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class RoIObj;
class HistogramManager;

class physics {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           RunNumber;
   Int_t           EventNumber;
   Int_t           timeStamp;
   Int_t           timeStampNSOffset;
   Int_t           lbn;
   Int_t           bcid;
   Int_t           detmask0;
   Int_t           detmask1;
   Float_t         actualIntPerXing;
   Float_t         averageIntPerXing;
   Int_t           pixelFlags;
   Int_t           sctFlags;
   Int_t           trtFlags;
   Int_t           larFlags;
   Int_t           tileFlags;
   Int_t           muonFlags;
   Int_t           fwdFlags;
   Int_t           coreFlags;
   Int_t           pixelError;
   Int_t           sctError;
   Int_t           trtError;
   Int_t           larError;
   Int_t           tileError;
   Int_t           muonError;
   Int_t           fwdError;
   Int_t           coreError;
   Int_t           trig_L1_mu_n;
   std::vector<float>   *trig_L1_mu_eta;
   std::vector<float>   *trig_L1_mu_phi;
   std::vector<std::string>  *trig_L1_mu_thrName;
   std::vector<short>   *trig_L1_mu_thrNumber;
   std::vector<short>   *trig_L1_mu_RoINumber;
   std::vector<short>   *trig_L1_mu_sectorAddress;
   std::vector<int>     *trig_L1_mu_firstCandidate;
   std::vector<int>     *trig_L1_mu_moreCandInRoI;
   std::vector<int>     *trig_L1_mu_moreCandInSector;
   std::vector<short>   *trig_L1_mu_source;
   std::vector<short>   *trig_L1_mu_hemisphere;
   std::vector<short>   *trig_L1_mu_charge;
   std::vector<int>     *trig_L1_mu_vetoed;
   Int_t           mu_n;
   std::vector<float>   *mu_pt;
   std::vector<float>   *mu_eta;
   std::vector<float>   *mu_phi;
   std::vector<float>   *mu_m;
   std::vector<int>     *mu_charge;
   std::vector<int>     *mu_author;
   std::vector<unsigned short> *mu_allAuthors;
   std::vector<int>     *mu_muonType;
   std::vector<float>   *mu_etcone20;
   std::vector<float>   *mu_etcone30;
   std::vector<float>   *mu_etcone40;
   std::vector<float>   *mu_ptcone20;
   std::vector<float>   *mu_ptcone30;
   std::vector<float>   *mu_ptcone40;
   std::vector<float>   *mu_trackfitchi2;
   std::vector<float>   *mu_trackfitndof;
   std::vector<int>     *mu_isPassedMCP;
   std::vector<int>     *mu_quality;
   std::vector<float>   *mu_msInnerMatchChi2;
   std::vector<float>   *mu_msOuterMatchChi2;
   std::vector<int>     *mu_msInnerMatchDOF;
   std::vector<int>     *mu_msOuterMatchDOF;
   std::vector<int>     *mu_nOutliersOnTrack;
   std::vector<int>     *mu_nBLHits;
   std::vector<int>     *mu_nPixHits;
   std::vector<int>     *mu_nSCTHits;
   std::vector<int>     *mu_nTRTHits;
   std::vector<int>     *mu_nTRTHighTHits;
   std::vector<int>     *mu_nBLSharedHits;
   std::vector<int>     *mu_nPixSharedHits;
   std::vector<int>     *mu_nPixHoles;
   std::vector<int>     *mu_nSCTSharedHits;
   std::vector<int>     *mu_nSCTHoles;
   std::vector<int>     *mu_nTRTOutliers;
   std::vector<int>     *mu_nTRTHighTOutliers;
   std::vector<int>     *mu_nGangedPixels;
   std::vector<int>     *mu_nPixelDeadSensors;
   std::vector<int>     *mu_nSCTDeadSensors;
   std::vector<int>     *mu_nTRTDeadStraws;
   std::vector<int>     *mu_expectBLayerHit;
   std::vector<int>     *mu_nPrecisionLayers;
   std::vector<int>     *mu_nPrecisionHoleLayers;
   std::vector<int>     *mu_nPhiLayers;
   std::vector<int>     *mu_nPhiHoleLayers;
   std::vector<int>     *mu_nTrigEtaLayers;
   std::vector<int>     *mu_nTrigEtaHoleLayers;
   std::vector<int>     *mu_primarySector;
   std::vector<int>     *mu_secondarySector;
   std::vector<int>     *mu_nInnerSmallHits;
   std::vector<int>     *mu_nInnerLargeHits;
   std::vector<int>     *mu_nMiddleSmallHits;
   std::vector<int>     *mu_nMiddleLargeHits;
   std::vector<int>     *mu_nOuterSmallHits;
   std::vector<int>     *mu_nOuterLargeHits;
   std::vector<int>     *mu_nExtendedSmallHits;
   std::vector<int>     *mu_nExtendedLargeHits;
   std::vector<int>     *mu_nInnerSmallHoles;
   std::vector<int>     *mu_nInnerLargeHoles;
   std::vector<int>     *mu_nMiddleSmallHoles;
   std::vector<int>     *mu_nMiddleLargeHoles;
   std::vector<int>     *mu_nOuterSmallHoles;
   std::vector<int>     *mu_nOuterLargeHoles;
   std::vector<int>     *mu_nExtendedSmallHoles;
   std::vector<int>     *mu_nExtendedLargeHoles;
   std::vector<int>     *mu_nPhiLayer1Hits;
   std::vector<int>     *mu_nPhiLayer2Hits;
   std::vector<int>     *mu_nPhiLayer3Hits;
   std::vector<int>     *mu_nPhiLayer4Hits;
   std::vector<int>     *mu_nEtaLayer1Hits;
   std::vector<int>     *mu_nEtaLayer2Hits;
   std::vector<int>     *mu_nEtaLayer3Hits;
   std::vector<int>     *mu_nEtaLayer4Hits;
   std::vector<int>     *mu_nPhiLayer1Holes;
   std::vector<int>     *mu_nPhiLayer2Holes;
   std::vector<int>     *mu_nPhiLayer3Holes;
   std::vector<int>     *mu_nPhiLayer4Holes;
   std::vector<int>     *mu_nEtaLayer1Holes;
   std::vector<int>     *mu_nEtaLayer2Holes;
   std::vector<int>     *mu_nEtaLayer3Holes;
   std::vector<int>     *mu_nEtaLayer4Holes;
   std::vector<float>   *mu_cb_d0;
   std::vector<float>   *mu_cb_z0;
   std::vector<float>   *mu_cb_phi0;
   std::vector<float>   *mu_cb_theta;
   std::vector<float>   *mu_cb_qOverP;
   std::vector<float>   *mu_cb_vx;
   std::vector<float>   *mu_cb_vy;
   std::vector<float>   *mu_cb_vz;
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
   Int_t           ext_mu_bias_n;
   std::vector<int>     *ext_mu_bias_type;
   std::vector<int>     *ext_mu_bias_index;
   std::vector<int>     *ext_mu_bias_size;
   std::vector<std::vector<int> > *ext_mu_bias_targetVec;
   std::vector<std::vector<float> > *ext_mu_bias_targetDistanceVec;
   std::vector<std::vector<float> > *ext_mu_bias_targetEtaVec;
   std::vector<std::vector<float> > *ext_mu_bias_targetPhiVec;
   std::vector<std::vector<float> > *ext_mu_bias_targetDeltaEtaVec;
   std::vector<std::vector<float> > *ext_mu_bias_targetDeltaPhiVec;
   std::vector<std::vector<float> > *ext_mu_bias_targetPxVec;
   std::vector<std::vector<float> > *ext_mu_bias_targetPyVec;
   std::vector<std::vector<float> > *ext_mu_bias_targetPzVec;
   Int_t           ext_mu_ubias_n;
   std::vector<int>     *ext_mu_ubias_type;
   std::vector<int>     *ext_mu_ubias_index;
   std::vector<int>     *ext_mu_ubias_size;
   std::vector<std::vector<int> > *ext_mu_ubias_targetVec;
   std::vector<std::vector<float> > *ext_mu_ubias_targetDistanceVec;
   std::vector<std::vector<float> > *ext_mu_ubias_targetEtaVec;
   std::vector<std::vector<float> > *ext_mu_ubias_targetPhiVec;
   std::vector<std::vector<float> > *ext_mu_ubias_targetDeltaEtaVec;
   std::vector<std::vector<float> > *ext_mu_ubias_targetDeltaPhiVec;
   std::vector<std::vector<float> > *ext_mu_ubias_targetPxVec;
   std::vector<std::vector<float> > *ext_mu_ubias_targetPyVec;
   std::vector<std::vector<float> > *ext_mu_ubias_targetPzVec;
   Int_t           trigger_info_n;
   std::vector<std::string>  *trigger_info_chain;
   std::vector<int>     *trigger_info_isPassed;
   std::vector<int>     *trigger_info_nTracks;
   std::vector<std::vector<int> > *trigger_info_typeVec;
   std::vector<std::vector<float> > *trigger_info_ptVec;
   std::vector<std::vector<float> > *trigger_info_etaVec;
   std::vector<std::vector<float> > *trigger_info_phiVec;
   Int_t           vxp_n;
   std::vector<float>   *vxp_x;
   std::vector<float>   *vxp_y;
   std::vector<float>   *vxp_z;
   std::vector<float>   *vxp_cov_x;
   std::vector<float>   *vxp_cov_y;
   std::vector<float>   *vxp_cov_z;
   std::vector<float>   *vxp_cov_xy;
   std::vector<float>   *vxp_cov_xz;
   std::vector<float>   *vxp_cov_yz;
   std::vector<float>   *vxp_chi2;
   std::vector<int>     *vxp_ndof;
   std::vector<int>     *vxp_nTracks;
   std::vector<int>     *vxp_type;
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
   Int_t           MDT_prd_n;
   std::vector<float>   *MDT_prd_x;
   std::vector<float>   *MDT_prd_y;
   std::vector<float>   *MDT_prd_z;
   std::vector<int>     *MDT_prd_adc;
   std::vector<int>     *MDT_prd_tdc;
   std::vector<int>     *MDT_prd_status;
   std::vector<float>   *MDT_prd_drift_radius;
   std::vector<float>   *MDT_prd_drift_radius_error;
   Int_t           TGC_coin_n;
   std::vector<float>   *TGC_coin_x_In;
   std::vector<float>   *TGC_coin_y_In;
   std::vector<float>   *TGC_coin_z_In;
   std::vector<float>   *TGC_coin_x_Out;
   std::vector<float>   *TGC_coin_y_Out;
   std::vector<float>   *TGC_coin_z_Out;
   std::vector<float>   *TGC_coin_width_In;
   std::vector<float>   *TGC_coin_width_Out;
   std::vector<float>   *TGC_coin_width_R;
   std::vector<float>   *TGC_coin_width_Phi;
   std::vector<int>     *TGC_coin_isAside;
   std::vector<int>     *TGC_coin_isForward;
   std::vector<int>     *TGC_coin_isStrip;
   std::vector<int>     *TGC_coin_isInner;
   std::vector<int>     *TGC_coin_isPositiveDeltaR;
   std::vector<int>     *TGC_coin_type;
   std::vector<int>     *TGC_coin_trackletId;
   std::vector<int>     *TGC_coin_trackletIdStrip;
   std::vector<int>     *TGC_coin_phi;
   std::vector<int>     *TGC_coin_roi;
   std::vector<int>     *TGC_coin_pt;
   std::vector<int>     *TGC_coin_delta;
   std::vector<int>     *TGC_coin_sub;
   std::vector<int>     *TGC_coin_veto;
   std::vector<int>     *TGC_coin_bunch;
   std::vector<int>     *TGC_coin_inner;
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
   Int_t           TGC_hierarchy_n;
   std::vector<int>     *TGC_hierarchy_index;
   std::vector<int>     *TGC_hierarchy_dR_hiPt;
   std::vector<int>     *TGC_hierarchy_dPhi_hiPt;
   std::vector<int>     *TGC_hierarchy_dR_tracklet;
   std::vector<int>     *TGC_hierarchy_dPhi_tracklet;
   std::vector<int>     *TGC_hierarchy_isChamberBoundary;
   std::vector<unsigned int> *muctpi_candidateMultiplicities;
   Int_t           muctpi_nDataWords;
   std::vector<unsigned int> *muctpi_dataWords;
   std::vector<float>   *muctpi_dw_eta;
   std::vector<float>   *muctpi_dw_phi;
   std::vector<short>   *muctpi_dw_source;
   std::vector<short>   *muctpi_dw_hemisphere;
   std::vector<short>   *muctpi_dw_bcid;
   std::vector<short>   *muctpi_dw_sectorID;
   std::vector<short>   *muctpi_dw_thrNumber;
   std::vector<short>   *muctpi_dw_roi;
   std::vector<short>   *muctpi_dw_veto;
   std::vector<short>   *muctpi_dw_firstCandidate;
   std::vector<short>   *muctpi_dw_moreCandInRoI;
   std::vector<short>   *muctpi_dw_moreCandInSector;
   std::vector<short>   *muctpi_dw_charge;
   std::vector<short>   *muctpi_dw_candidateVetoed;
   Int_t           TILE_cell_n;
   std::vector<float>   *TILE_cell_E;
   std::vector<float>   *TILE_cell_ene1;
   std::vector<float>   *TILE_cell_eta;
   std::vector<float>   *TILE_cell_phi;
   std::vector<float>   *TILE_cell_sinTh;
   std::vector<float>   *TILE_cell_cosTh;
   std::vector<float>   *TILE_cell_cotTh;
   std::vector<float>   *TILE_cell_x;
   std::vector<float>   *TILE_cell_y;
   std::vector<float>   *TILE_cell_z;
   std::vector<int>     *TILE_cell_badcell;
   std::vector<int>     *TILE_cell_partition;
   std::vector<int>     *TILE_cell_section;
   std::vector<int>     *TILE_cell_side;
   std::vector<int>     *TILE_cell_module;
   std::vector<int>     *TILE_cell_tower;
   std::vector<int>     *TILE_cell_sample;
   //Truth
   Int_t           mc_n;
   std::vector<float>   *mc_pt;
   std::vector<float>   *mc_eta;
   std::vector<float>   *mc_phi;
   std::vector<float>   *mc_m;
   std::vector<float>   *mc_charge;


   // List of branches
   TBranch        *b_runNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_timeStamp;   //!
   TBranch        *b_timeStampNSOffset;   //!
   TBranch        *b_lbn;   //!
   TBranch        *b_bcid;   //!
   TBranch        *b_detmask0;   //!
   TBranch        *b_detmask1;   //!
   TBranch        *b_actualIntPerXing;   //!
   TBranch        *b_averageIntPerXing;   //!
   TBranch        *b_pixelFlags;   //!
   TBranch        *b_sctFlags;   //!
   TBranch        *b_trtFlags;   //!
   TBranch        *b_larFlags;   //!
   TBranch        *b_tileFlags;   //!
   TBranch        *b_muonFlags;   //!
   TBranch        *b_fwdFlags;   //!
   TBranch        *b_coreFlags;   //!
   TBranch        *b_pixelError;   //!
   TBranch        *b_sctError;   //!
   TBranch        *b_trtError;   //!
   TBranch        *b_larError;   //!
   TBranch        *b_tileError;   //!
   TBranch        *b_muonError;   //!
   TBranch        *b_fwdError;   //!
   TBranch        *b_coreError;   //!
   TBranch        *b_trig_L1_mu_n;   //!
   TBranch        *b_trig_L1_mu_eta;   //!
   TBranch        *b_trig_L1_mu_phi;   //!
   TBranch        *b_trig_L1_mu_thrName;   //!
   TBranch        *b_trig_L1_mu_thrNumber;   //!
   TBranch        *b_trig_L1_mu_RoINumber;   //!
   TBranch        *b_trig_L1_mu_sectorAddress;   //!
   TBranch        *b_trig_L1_mu_firstCandidate;   //!
   TBranch        *b_trig_L1_mu_moreCandInRoI;   //!
   TBranch        *b_trig_L1_mu_moreCandInSector;   //!
   TBranch        *b_trig_L1_mu_source;   //!
   TBranch        *b_trig_L1_mu_hemisphere;   //!
   TBranch        *b_trig_L1_mu_charge;   //!
   TBranch        *b_trig_L1_mu_vetoed;   //!
   TBranch        *b_mu_n;   //!
   TBranch        *b_mu_pt;   //!
   TBranch        *b_mu_eta;   //!
   TBranch        *b_mu_phi;   //!
   TBranch        *b_mu_m;   //!
   TBranch        *b_mu_charge;   //!
   TBranch        *b_mu_author;   //!
   TBranch        *b_mu_allAuthors;   //!
   TBranch        *b_mu_muonType;   //!
   TBranch        *b_mu_etcone20;   //!
   TBranch        *b_mu_etcone30;   //!
   TBranch        *b_mu_etcone40;   //!
   TBranch        *b_mu_ptcone20;   //!
   TBranch        *b_mu_ptcone30;   //!
   TBranch        *b_mu_ptcone40;   //!
   TBranch        *b_mu_trackfitchi2;   //!
   TBranch        *b_mu_trackfitndof;   //!
   TBranch        *b_mu_isPassedMCP;   //!
   TBranch        *b_mu_quality;   //!
   TBranch        *b_mu_msInnerMatchChi2;   //!
   TBranch        *b_mu_msOuterMatchChi2;   //!
   TBranch        *b_mu_msInnerMatchDOF;   //!
   TBranch        *b_mu_msOuterMatchDOF;   //!
   TBranch        *b_mu_nOutliersOnTrack;   //!
   TBranch        *b_mu_nBLHits;   //!
   TBranch        *b_mu_nPixHits;   //!
   TBranch        *b_mu_nSCTHits;   //!
   TBranch        *b_mu_nTRTHits;   //!
   TBranch        *b_mu_nTRTHighTHits;   //!
   TBranch        *b_mu_nBLSharedHits;   //!
   TBranch        *b_mu_nPixSharedHits;   //!
   TBranch        *b_mu_nPixHoles;   //!
   TBranch        *b_mu_nSCTSharedHits;   //!
   TBranch        *b_mu_nSCTHoles;   //!
   TBranch        *b_mu_nTRTOutliers;   //!
   TBranch        *b_mu_nTRTHighTOutliers;   //!
   TBranch        *b_mu_nGangedPixels;   //!
   TBranch        *b_mu_nPixelDeadSensors;   //!
   TBranch        *b_mu_nSCTDeadSensors;   //!
   TBranch        *b_mu_nTRTDeadStraws;   //!
   TBranch        *b_mu_expectBLayerHit;   //!
   TBranch        *b_mu_nPrecisionLayers;   //!
   TBranch        *b_mu_nPrecisionHoleLayers;   //!
   TBranch        *b_mu_nPhiLayers;   //!
   TBranch        *b_mu_nPhiHoleLayers;   //!
   TBranch        *b_mu_nTrigEtaLayers;   //!
   TBranch        *b_mu_nTrigEtaHoleLayers;   //!
   TBranch        *b_mu_primarySector;   //!
   TBranch        *b_mu_secondarySector;   //!
   TBranch        *b_mu_nInnerSmallHits;   //!
   TBranch        *b_mu_nInnerLargeHits;   //!
   TBranch        *b_mu_nMiddleSmallHits;   //!
   TBranch        *b_mu_nMiddleLargeHits;   //!
   TBranch        *b_mu_nOuterSmallHits;   //!
   TBranch        *b_mu_nOuterLargeHits;   //!
   TBranch        *b_mu_nExtendedSmallHits;   //!
   TBranch        *b_mu_nExtendedLargeHits;   //!
   TBranch        *b_mu_nInnerSmallHoles;   //!
   TBranch        *b_mu_nInnerLargeHoles;   //!
   TBranch        *b_mu_nMiddleSmallHoles;   //!
   TBranch        *b_mu_nMiddleLargeHoles;   //!
   TBranch        *b_mu_nOuterSmallHoles;   //!
   TBranch        *b_mu_nOuterLargeHoles;   //!
   TBranch        *b_mu_nExtendedSmallHoles;   //!
   TBranch        *b_mu_nExtendedLargeHoles;   //!
   TBranch        *b_mu_nPhiLayer1Hits;   //!
   TBranch        *b_mu_nPhiLayer2Hits;   //!
   TBranch        *b_mu_nPhiLayer3Hits;   //!
   TBranch        *b_mu_nPhiLayer4Hits;   //!
   TBranch        *b_mu_nEtaLayer1Hits;   //!
   TBranch        *b_mu_nEtaLayer2Hits;   //!
   TBranch        *b_mu_nEtaLayer3Hits;   //!
   TBranch        *b_mu_nEtaLayer4Hits;   //!
   TBranch        *b_mu_nPhiLayer1Holes;   //!
   TBranch        *b_mu_nPhiLayer2Holes;   //!
   TBranch        *b_mu_nPhiLayer3Holes;   //!
   TBranch        *b_mu_nPhiLayer4Holes;   //!
   TBranch        *b_mu_nEtaLayer1Holes;   //!
   TBranch        *b_mu_nEtaLayer2Holes;   //!
   TBranch        *b_mu_nEtaLayer3Holes;   //!
   TBranch        *b_mu_nEtaLayer4Holes;   //!
   TBranch        *b_mu_cb_d0;   //!
   TBranch        *b_mu_cb_z0;   //!
   TBranch        *b_mu_cb_phi0;   //!
   TBranch        *b_mu_cb_theta;   //!
   TBranch        *b_mu_cb_qOverP;   //!
   TBranch        *b_mu_cb_vx;   //!
   TBranch        *b_mu_cb_vy;   //!
   TBranch        *b_mu_cb_vz;   //!
   TBranch        *b_TGC_museg_n;   //!
   TBranch        *b_museg_x;   //!
   TBranch        *b_museg_y;   //!
   TBranch        *b_museg_z;   //!
   TBranch        *b_museg_px;   //!
   TBranch        *b_museg_py;   //!
   TBranch        *b_museg_pz;   //!
   TBranch        *b_museg_t0;   //!
   TBranch        *b_museg_t0error;   //!
   TBranch        *b_museg_chi2;   //!
   TBranch        *b_museg_ndof;   //!
   TBranch        *b_museg_sector;   //!
   TBranch        *b_museg_stationName;   //!
   TBranch        *b_museg_stationEta;   //!
   TBranch        *b_museg_author;   //!
   TBranch        *b_ext_mu_bias_n;   //!
   TBranch        *b_ext_mu_bias_type;   //!
   TBranch        *b_ext_mu_bias_index;   //!
   TBranch        *b_ext_mu_bias_size;   //!
   TBranch        *b_ext_mu_bias_targetVec;   //!
   TBranch        *b_ext_mu_bias_targetDistanceVec;   //!
   TBranch        *b_ext_mu_bias_targetEtaVec;   //!
   TBranch        *b_ext_mu_bias_targetPhiVec;   //!
   TBranch        *b_ext_mu_bias_targetDeltaEtaVec;   //!
   TBranch        *b_ext_mu_bias_targetDeltaPhiVec;   //!
   TBranch        *b_ext_mu_bias_targetPxVec;   //!
   TBranch        *b_ext_mu_bias_targetPyVec;   //!
   TBranch        *b_ext_mu_bias_targetPzVec;   //!
   TBranch        *b_ext_mu_ubias_n;   //!
   TBranch        *b_ext_mu_ubias_type;   //!
   TBranch        *b_ext_mu_ubias_index;   //!
   TBranch        *b_ext_mu_ubias_size;   //!
   TBranch        *b_ext_mu_ubias_targetVec;   //!
   TBranch        *b_ext_mu_ubias_targetDistanceVec;   //!
   TBranch        *b_ext_mu_ubias_targetEtaVec;   //!
   TBranch        *b_ext_mu_ubias_targetPhiVec;   //!
   TBranch        *b_ext_mu_ubias_targetDeltaEtaVec;   //!
   TBranch        *b_ext_mu_ubias_targetDeltaPhiVec;   //!
   TBranch        *b_ext_mu_ubias_targetPxVec;   //!
   TBranch        *b_ext_mu_ubias_targetPyVec;   //!
   TBranch        *b_ext_mu_ubias_targetPzVec;   //!
   TBranch        *b_trigger_info_n;   //!
   TBranch        *b_trigger_info_chain;   //!
   TBranch        *b_trigger_info_isPassed;   //!
   TBranch        *b_trigger_info_nTracks;   //!
   TBranch        *b_trigger_info_typeVec;   //!
   TBranch        *b_trigger_info_ptVec;   //!
   TBranch        *b_trigger_info_etaVec;   //!
   TBranch        *b_trigger_info_phiVec;   //!
   TBranch        *b_vxp_n;   //!
   TBranch        *b_vxp_x;   //!
   TBranch        *b_vxp_y;   //!
   TBranch        *b_vxp_z;   //!
   TBranch        *b_vxp_cov_x;   //!
   TBranch        *b_vxp_cov_y;   //!
   TBranch        *b_vxp_cov_z;   //!
   TBranch        *b_vxp_cov_xy;   //!
   TBranch        *b_vxp_cov_xz;   //!
   TBranch        *b_vxp_cov_yz;   //!
   TBranch        *b_vxp_chi2;   //!
   TBranch        *b_vxp_ndof;   //!
   TBranch        *b_vxp_nTracks;   //!
   TBranch        *b_vxp_type;   //!
   TBranch        *b_TGC_prd_n;   //!
   TBranch        *b_TGC_prd_x;   //!
   TBranch        *b_TGC_prd_y;   //!
   TBranch        *b_TGC_prd_z;   //!
   TBranch        *b_TGC_prd_shortWidth;   //!
   TBranch        *b_TGC_prd_longWidth;   //!
   TBranch        *b_TGC_prd_length;   //!
   TBranch        *b_TGC_prd_isStrip;   //!
   TBranch        *b_TGC_prd_gasGap;   //!
   TBranch        *b_TGC_prd_channel;   //!
   TBranch        *b_TGC_prd_eta;   //!
   TBranch        *b_TGC_prd_phi;   //!
   TBranch        *b_TGC_prd_station;   //!
   TBranch        *b_TGC_prd_bunch;   //!
   TBranch        *b_RPC_prd_n;   //!
   TBranch        *b_RPC_prd_x;   //!
   TBranch        *b_RPC_prd_y;   //!
   TBranch        *b_RPC_prd_z;   //!
   TBranch        *b_RPC_prd_x2;   //!
   TBranch        *b_RPC_prd_y2;   //!
   TBranch        *b_RPC_prd_z2;   //!
   TBranch        *b_RPC_prd_time;   //!
   TBranch        *b_RPC_prd_triggerInfo;   //!
   TBranch        *b_RPC_prd_ambiguityFlag;   //!
   TBranch        *b_RPC_prd_measuresPhi;   //!
   TBranch        *b_RPC_prd_inRibs;   //!
   TBranch        *b_RPC_prd_station;   //!
   TBranch        *b_RPC_prd_stationEta;   //!
   TBranch        *b_RPC_prd_stationPhi;   //!
   TBranch        *b_RPC_prd_doubletR;   //!
   TBranch        *b_RPC_prd_doubletZ;   //!
   TBranch        *b_RPC_prd_stripWidth;   //!
   TBranch        *b_RPC_prd_stripLength;   //!
   TBranch        *b_RPC_prd_gasGap;   //!
   TBranch        *b_RPC_prd_channel;   //!
   TBranch        *b_MDT_prd_n;   //!
   TBranch        *b_MDT_prd_x;   //!
   TBranch        *b_MDT_prd_y;   //!
   TBranch        *b_MDT_prd_z;   //!
   TBranch        *b_MDT_prd_adc;   //!
   TBranch        *b_MDT_prd_tdc;   //!
   TBranch        *b_MDT_prd_status;   //!
   TBranch        *b_MDT_prd_drift_radius;   //!
   TBranch        *b_MDT_prd_drift_radius_error;   //!
   TBranch        *b_TGC_coin_n;   //!
   TBranch        *b_TGC_coin_x_In;   //!
   TBranch        *b_TGC_coin_y_In;   //!
   TBranch        *b_TGC_coin_z_In;   //!
   TBranch        *b_TGC_coin_x_Out;   //!
   TBranch        *b_TGC_coin_y_Out;   //!
   TBranch        *b_TGC_coin_z_Out;   //!
   TBranch        *b_TGC_coin_width_In;   //!
   TBranch        *b_TGC_coin_width_Out;   //!
   TBranch        *b_TGC_coin_width_R;   //!
   TBranch        *b_TGC_coin_width_Phi;   //!
   TBranch        *b_TGC_coin_isAside;   //!
   TBranch        *b_TGC_coin_isForward;   //!
   TBranch        *b_TGC_coin_isStrip;   //!
   TBranch        *b_TGC_coin_isInner;   //!
   TBranch        *b_TGC_coin_isPositiveDeltaR;   //!
   TBranch        *b_TGC_coin_type;   //!
   TBranch        *b_TGC_coin_trackletId;   //!
   TBranch        *b_TGC_coin_trackletIdStrip;   //!
   TBranch        *b_TGC_coin_phi;   //!
   TBranch        *b_TGC_coin_roi;   //!
   TBranch        *b_TGC_coin_pt;   //!
   TBranch        *b_TGC_coin_delta;   //!
   TBranch        *b_TGC_coin_sub;   //!
   TBranch        *b_TGC_coin_veto;   //!
   TBranch        *b_TGC_coin_bunch;   //!
   TBranch        *b_TGC_coin_inner;   //!
   TBranch        *b_TILE_murcv_trig_n;   //!
   TBranch        *b_TILE_murcv_trig_mod;   //!
   TBranch        *b_TILE_murcv_trig_part;   //!
   TBranch        *b_TILE_murcv_trig_bit0;   //!
   TBranch        *b_TILE_murcv_trig_bit1;   //!
   TBranch        *b_TILE_murcv_trig_bit2;   //!
   TBranch        *b_TILE_murcv_trig_bit3;   //!
   TBranch        *b_TILE_murcv_raw_n;   //!
   TBranch        *b_TILE_murcv_raw_count;   //!
   TBranch        *b_TILE_murcv_raw_energy;   //!
   TBranch        *b_TILE_murcv_raw_ros;   //!
   TBranch        *b_TILE_murcv_raw_drawer;   //!
   TBranch        *b_TILE_murcv_raw_channel;   //!
   TBranch        *b_TILE_murcv_digit_n;   //!
   TBranch        *b_TILE_murcv_digit_nSamples;   //!
   TBranch        *b_TILE_murcv_digit_ros;   //!
   TBranch        *b_TILE_murcv_digit_drawer;   //!
   TBranch        *b_TILE_murcv_digit_channel;   //!
   TBranch        *b_TILE_murcv_digit_sampleVec;   //!
   TBranch        *b_TGC_hierarchy_n;   //!
   TBranch        *b_TGC_hierarchy_index;   //!
   TBranch        *b_TGC_hierarchy_dR_hiPt;   //!
   TBranch        *b_TGC_hierarchy_dPhi_hiPt;   //!
   TBranch        *b_TGC_hierarchy_dR_tracklet;   //!
   TBranch        *b_TGC_hierarchy_dPhi_tracklet;   //!
   TBranch        *b_TGC_hierarchy_isChamberBoundary;   //!
   TBranch        *b_muctpi_candidateMultiplicities;   //!
   TBranch        *b_muctpi_nDataWords;   //!
   TBranch        *b_muctpi_dataWords;   //!
   TBranch        *b_muctpi_dw_eta;   //!
   TBranch        *b_muctpi_dw_phi;   //!
   TBranch        *b_muctpi_dw_source;   //!
   TBranch        *b_muctpi_dw_hemisphere;   //!
   TBranch        *b_muctpi_dw_bcid;   //!
   TBranch        *b_muctpi_dw_sectorID;   //!
   TBranch        *b_muctpi_dw_thrNumber;   //!
   TBranch        *b_muctpi_dw_roi;   //!
   TBranch        *b_muctpi_dw_veto;   //!
   TBranch        *b_muctpi_dw_firstCandidate;   //!
   TBranch        *b_muctpi_dw_moreCandInRoI;   //!
   TBranch        *b_muctpi_dw_moreCandInSector;   //!
   TBranch        *b_muctpi_dw_charge;   //!
   TBranch        *b_muctpi_dw_candidateVetoed;   //!
   TBranch        *b_TILE_cell_n;   //!
   TBranch        *b_TILE_cell_E;   //!
   TBranch        *b_TILE_cell_ene1;   //!
   TBranch        *b_TILE_cell_eta;   //!
   TBranch        *b_TILE_cell_phi;   //!
   TBranch        *b_TILE_cell_sinTh;   //!
   TBranch        *b_TILE_cell_cosTh;   //!
   TBranch        *b_TILE_cell_cotTh;   //!
   TBranch        *b_TILE_cell_x;   //!
   TBranch        *b_TILE_cell_y;   //!
   TBranch        *b_TILE_cell_z;   //!
   TBranch        *b_TILE_cell_badcell;   //!
   TBranch        *b_TILE_cell_partition;   //!
   TBranch        *b_TILE_cell_section;   //!
   TBranch        *b_TILE_cell_side;   //!
   TBranch        *b_TILE_cell_module;   //!
   TBranch        *b_TILE_cell_tower;   //!
   TBranch        *b_TILE_cell_sample;   //!
   TBranch        *b_mc_n;
   TBranch        *b_mc_pt;
   TBranch        *b_mc_eta;
   TBranch        *b_mc_phi;
   TBranch        *b_mc_m;
   TBranch        *b_mc_charge;

   physics(TTree *tree=0);
   virtual ~physics();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(std::string BWCW_path);
   virtual void     FillNtuple(HistogramManager &hs, std::vector<RoIObj>& RoIs);
   virtual Bool_t   RoISL(std::vector<RoIObj>& RoIs);
   virtual Int_t    getRoI(int wid, int sid);
   virtual void     RoIDRDPhi(std::vector<RoIObj>& RoIs);
   virtual void     RoIHPTLPT(std::vector<RoIObj>& RoIs);
   virtual void     RoITwoToFour(std::vector<RoIObj>& RoIs);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef physics_cxx
physics::physics(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/gpfs/fs7001/shiomi/MCData/Run-2/Residual/LiTGCNtuple/SingleMuonEB/0608/user.shiomi.0601.inTruth.singlemuon.v2.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/gpfs/fs7001/shiomi/MCData/Run-2/Residual/LiTGCNtuple/SingleMuonEB/0608/user.shiomi.0601.inTruth.singlemuon.v2.root");
      }
      f->GetObject("physics",tree);

   }
   Init(tree);
}

physics::~physics()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t physics::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t physics::LoadTree(Long64_t entry)
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

void physics::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   trig_L1_mu_eta = 0;
   trig_L1_mu_phi = 0;
   trig_L1_mu_thrName = 0;
   trig_L1_mu_thrNumber = 0;
   trig_L1_mu_RoINumber = 0;
   trig_L1_mu_sectorAddress = 0;
   trig_L1_mu_firstCandidate = 0;
   trig_L1_mu_moreCandInRoI = 0;
   trig_L1_mu_moreCandInSector = 0;
   trig_L1_mu_source = 0;
   trig_L1_mu_hemisphere = 0;
   trig_L1_mu_charge = 0;
   trig_L1_mu_vetoed = 0;
   mu_pt = 0;
   mu_eta = 0;
   mu_phi = 0;
   mu_m = 0;
   mu_charge = 0;
   mu_author = 0;
   mu_allAuthors = 0;
   mu_muonType = 0;
   mu_etcone20 = 0;
   mu_etcone30 = 0;
   mu_etcone40 = 0;
   mu_ptcone20 = 0;
   mu_ptcone30 = 0;
   mu_ptcone40 = 0;
   mu_trackfitchi2 = 0;
   mu_trackfitndof = 0;
   mu_isPassedMCP = 0;
   mu_quality = 0;
   mu_msInnerMatchChi2 = 0;
   mu_msOuterMatchChi2 = 0;
   mu_msInnerMatchDOF = 0;
   mu_msOuterMatchDOF = 0;
   mu_nOutliersOnTrack = 0;
   mu_nBLHits = 0;
   mu_nPixHits = 0;
   mu_nSCTHits = 0;
   mu_nTRTHits = 0;
   mu_nTRTHighTHits = 0;
   mu_nBLSharedHits = 0;
   mu_nPixSharedHits = 0;
   mu_nPixHoles = 0;
   mu_nSCTSharedHits = 0;
   mu_nSCTHoles = 0;
   mu_nTRTOutliers = 0;
   mu_nTRTHighTOutliers = 0;
   mu_nGangedPixels = 0;
   mu_nPixelDeadSensors = 0;
   mu_nSCTDeadSensors = 0;
   mu_nTRTDeadStraws = 0;
   mu_expectBLayerHit = 0;
   mu_nPrecisionLayers = 0;
   mu_nPrecisionHoleLayers = 0;
   mu_nPhiLayers = 0;
   mu_nPhiHoleLayers = 0;
   mu_nTrigEtaLayers = 0;
   mu_nTrigEtaHoleLayers = 0;
   mu_primarySector = 0;
   mu_secondarySector = 0;
   mu_nInnerSmallHits = 0;
   mu_nInnerLargeHits = 0;
   mu_nMiddleSmallHits = 0;
   mu_nMiddleLargeHits = 0;
   mu_nOuterSmallHits = 0;
   mu_nOuterLargeHits = 0;
   mu_nExtendedSmallHits = 0;
   mu_nExtendedLargeHits = 0;
   mu_nInnerSmallHoles = 0;
   mu_nInnerLargeHoles = 0;
   mu_nMiddleSmallHoles = 0;
   mu_nMiddleLargeHoles = 0;
   mu_nOuterSmallHoles = 0;
   mu_nOuterLargeHoles = 0;
   mu_nExtendedSmallHoles = 0;
   mu_nExtendedLargeHoles = 0;
   mu_nPhiLayer1Hits = 0;
   mu_nPhiLayer2Hits = 0;
   mu_nPhiLayer3Hits = 0;
   mu_nPhiLayer4Hits = 0;
   mu_nEtaLayer1Hits = 0;
   mu_nEtaLayer2Hits = 0;
   mu_nEtaLayer3Hits = 0;
   mu_nEtaLayer4Hits = 0;
   mu_nPhiLayer1Holes = 0;
   mu_nPhiLayer2Holes = 0;
   mu_nPhiLayer3Holes = 0;
   mu_nPhiLayer4Holes = 0;
   mu_nEtaLayer1Holes = 0;
   mu_nEtaLayer2Holes = 0;
   mu_nEtaLayer3Holes = 0;
   mu_nEtaLayer4Holes = 0;
   mu_cb_d0 = 0;
   mu_cb_z0 = 0;
   mu_cb_phi0 = 0;
   mu_cb_theta = 0;
   mu_cb_qOverP = 0;
   mu_cb_vx = 0;
   mu_cb_vy = 0;
   mu_cb_vz = 0;
   museg_x = 0;
   museg_y = 0;
   museg_z = 0;
   museg_px = 0;
   museg_py = 0;
   museg_pz = 0;
   museg_t0 = 0;
   museg_t0error = 0;
   museg_chi2 = 0;
   museg_ndof = 0;
   museg_sector = 0;
   museg_stationName = 0;
   museg_stationEta = 0;
   museg_author = 0;
   ext_mu_bias_type = 0;
   ext_mu_bias_index = 0;
   ext_mu_bias_size = 0;
   ext_mu_bias_targetVec = 0;
   ext_mu_bias_targetDistanceVec = 0;
   ext_mu_bias_targetEtaVec = 0;
   ext_mu_bias_targetPhiVec = 0;
   ext_mu_bias_targetDeltaEtaVec = 0;
   ext_mu_bias_targetDeltaPhiVec = 0;
   ext_mu_bias_targetPxVec = 0;
   ext_mu_bias_targetPyVec = 0;
   ext_mu_bias_targetPzVec = 0;
   ext_mu_ubias_type = 0;
   ext_mu_ubias_index = 0;
   ext_mu_ubias_size = 0;
   ext_mu_ubias_targetVec = 0;
   ext_mu_ubias_targetDistanceVec = 0;
   ext_mu_ubias_targetEtaVec = 0;
   ext_mu_ubias_targetPhiVec = 0;
   ext_mu_ubias_targetDeltaEtaVec = 0;
   ext_mu_ubias_targetDeltaPhiVec = 0;
   ext_mu_ubias_targetPxVec = 0;
   ext_mu_ubias_targetPyVec = 0;
   ext_mu_ubias_targetPzVec = 0;
   trigger_info_chain = 0;
   trigger_info_isPassed = 0;
   trigger_info_nTracks = 0;
   trigger_info_typeVec = 0;
   trigger_info_ptVec = 0;
   trigger_info_etaVec = 0;
   trigger_info_phiVec = 0;
   vxp_x = 0;
   vxp_y = 0;
   vxp_z = 0;
   vxp_cov_x = 0;
   vxp_cov_y = 0;
   vxp_cov_z = 0;
   vxp_cov_xy = 0;
   vxp_cov_xz = 0;
   vxp_cov_yz = 0;
   vxp_chi2 = 0;
   vxp_ndof = 0;
   vxp_nTracks = 0;
   vxp_type = 0;
   TGC_prd_x = 0;
   TGC_prd_y = 0;
   TGC_prd_z = 0;
   TGC_prd_shortWidth = 0;
   TGC_prd_longWidth = 0;
   TGC_prd_length = 0;
   TGC_prd_isStrip = 0;
   TGC_prd_gasGap = 0;
   TGC_prd_channel = 0;
   TGC_prd_eta = 0;
   TGC_prd_phi = 0;
   TGC_prd_station = 0;
   TGC_prd_bunch = 0;
   RPC_prd_x = 0;
   RPC_prd_y = 0;
   RPC_prd_z = 0;
   RPC_prd_x2 = 0;
   RPC_prd_y2 = 0;
   RPC_prd_z2 = 0;
   RPC_prd_time = 0;
   RPC_prd_triggerInfo = 0;
   RPC_prd_ambiguityFlag = 0;
   RPC_prd_measuresPhi = 0;
   RPC_prd_inRibs = 0;
   RPC_prd_station = 0;
   RPC_prd_stationEta = 0;
   RPC_prd_stationPhi = 0;
   RPC_prd_doubletR = 0;
   RPC_prd_doubletZ = 0;
   RPC_prd_stripWidth = 0;
   RPC_prd_stripLength = 0;
   RPC_prd_gasGap = 0;
   RPC_prd_channel = 0;
   MDT_prd_x = 0;
   MDT_prd_y = 0;
   MDT_prd_z = 0;
   MDT_prd_adc = 0;
   MDT_prd_tdc = 0;
   MDT_prd_status = 0;
   MDT_prd_drift_radius = 0;
   MDT_prd_drift_radius_error = 0;
   TGC_coin_x_In = 0;
   TGC_coin_y_In = 0;
   TGC_coin_z_In = 0;
   TGC_coin_x_Out = 0;
   TGC_coin_y_Out = 0;
   TGC_coin_z_Out = 0;
   TGC_coin_width_In = 0;
   TGC_coin_width_Out = 0;
   TGC_coin_width_R = 0;
   TGC_coin_width_Phi = 0;
   TGC_coin_isAside = 0;
   TGC_coin_isForward = 0;
   TGC_coin_isStrip = 0;
   TGC_coin_isInner = 0;
   TGC_coin_isPositiveDeltaR = 0;
   TGC_coin_type = 0;
   TGC_coin_trackletId = 0;
   TGC_coin_trackletIdStrip = 0;
   TGC_coin_phi = 0;
   TGC_coin_roi = 0;
   TGC_coin_pt = 0;
   TGC_coin_delta = 0;
   TGC_coin_sub = 0;
   TGC_coin_veto = 0;
   TGC_coin_bunch = 0;
   TGC_coin_inner = 0;
   TILE_murcv_trig_mod = 0;
   TILE_murcv_trig_part = 0;
   TILE_murcv_trig_bit0 = 0;
   TILE_murcv_trig_bit1 = 0;
   TILE_murcv_trig_bit2 = 0;
   TILE_murcv_trig_bit3 = 0;
   TILE_murcv_raw_count = 0;
   TILE_murcv_raw_energy = 0;
   TILE_murcv_raw_ros = 0;
   TILE_murcv_raw_drawer = 0;
   TILE_murcv_raw_channel = 0;
   TILE_murcv_digit_nSamples = 0;
   TILE_murcv_digit_ros = 0;
   TILE_murcv_digit_drawer = 0;
   TILE_murcv_digit_channel = 0;
   TILE_murcv_digit_sampleVec = 0;
   TGC_hierarchy_index = 0;
   TGC_hierarchy_dR_hiPt = 0;
   TGC_hierarchy_dPhi_hiPt = 0;
   TGC_hierarchy_dR_tracklet = 0;
   TGC_hierarchy_dPhi_tracklet = 0;
   TGC_hierarchy_isChamberBoundary = 0;
   muctpi_candidateMultiplicities = 0;
   muctpi_dataWords = 0;
   muctpi_dw_eta = 0;
   muctpi_dw_phi = 0;
   muctpi_dw_source = 0;
   muctpi_dw_hemisphere = 0;
   muctpi_dw_bcid = 0;
   muctpi_dw_sectorID = 0;
   muctpi_dw_thrNumber = 0;
   muctpi_dw_roi = 0;
   muctpi_dw_veto = 0;
   muctpi_dw_firstCandidate = 0;
   muctpi_dw_moreCandInRoI = 0;
   muctpi_dw_moreCandInSector = 0;
   muctpi_dw_charge = 0;
   muctpi_dw_candidateVetoed = 0;
   TILE_cell_E = 0;
   TILE_cell_ene1 = 0;
   TILE_cell_eta = 0;
   TILE_cell_phi = 0;
   TILE_cell_sinTh = 0;
   TILE_cell_cosTh = 0;
   TILE_cell_cotTh = 0;
   TILE_cell_x = 0;
   TILE_cell_y = 0;
   TILE_cell_z = 0;
   TILE_cell_badcell = 0;
   TILE_cell_partition = 0;
   TILE_cell_section = 0;
   TILE_cell_side = 0;
   TILE_cell_module = 0;
   TILE_cell_tower = 0;
   TILE_cell_sample = 0;
   mc_pt=0;
   mc_eta=0;
   mc_phi=0;
   mc_m=0;
   mc_charge=0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_runNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("timeStamp", &timeStamp, &b_timeStamp);
   fChain->SetBranchAddress("timeStampNSOffset", &timeStampNSOffset, &b_timeStampNSOffset);
   fChain->SetBranchAddress("lbn", &lbn, &b_lbn);
   fChain->SetBranchAddress("bcid", &bcid, &b_bcid);
   /*
   fChain->SetBranchAddress("detmask0", &detmask0, &b_detmask0);
   fChain->SetBranchAddress("detmask1", &detmask1, &b_detmask1);
   fChain->SetBranchAddress("actualIntPerXing", &actualIntPerXing, &b_actualIntPerXing);
   fChain->SetBranchAddress("averageIntPerXing", &averageIntPerXing, &b_averageIntPerXing);
   fChain->SetBranchAddress("pixelFlags", &pixelFlags, &b_pixelFlags);
   fChain->SetBranchAddress("sctFlags", &sctFlags, &b_sctFlags);
   fChain->SetBranchAddress("trtFlags", &trtFlags, &b_trtFlags);
   fChain->SetBranchAddress("larFlags", &larFlags, &b_larFlags);
   fChain->SetBranchAddress("tileFlags", &tileFlags, &b_tileFlags);
   fChain->SetBranchAddress("muonFlags", &muonFlags, &b_muonFlags);
   fChain->SetBranchAddress("fwdFlags", &fwdFlags, &b_fwdFlags);
   fChain->SetBranchAddress("coreFlags", &coreFlags, &b_coreFlags);
   fChain->SetBranchAddress("pixelError", &pixelError, &b_pixelError);
   fChain->SetBranchAddress("sctError", &sctError, &b_sctError);
   fChain->SetBranchAddress("trtError", &trtError, &b_trtError);
   fChain->SetBranchAddress("larError", &larError, &b_larError);
   fChain->SetBranchAddress("tileError", &tileError, &b_tileError);
   fChain->SetBranchAddress("muonError", &muonError, &b_muonError);
   fChain->SetBranchAddress("fwdError", &fwdError, &b_fwdError);
   fChain->SetBranchAddress("coreError", &coreError, &b_coreError);
   */
   fChain->SetBranchAddress("trig_L1_mu_n", &trig_L1_mu_n, &b_trig_L1_mu_n);
   fChain->SetBranchAddress("trig_L1_mu_eta", &trig_L1_mu_eta, &b_trig_L1_mu_eta);
   fChain->SetBranchAddress("trig_L1_mu_phi", &trig_L1_mu_phi, &b_trig_L1_mu_phi);
   fChain->SetBranchAddress("trig_L1_mu_thrName", &trig_L1_mu_thrName, &b_trig_L1_mu_thrName);
   fChain->SetBranchAddress("trig_L1_mu_thrNumber", &trig_L1_mu_thrNumber, &b_trig_L1_mu_thrNumber);
   fChain->SetBranchAddress("trig_L1_mu_RoINumber", &trig_L1_mu_RoINumber, &b_trig_L1_mu_RoINumber);
   fChain->SetBranchAddress("trig_L1_mu_sectorAddress", &trig_L1_mu_sectorAddress, &b_trig_L1_mu_sectorAddress);
   fChain->SetBranchAddress("trig_L1_mu_firstCandidate", &trig_L1_mu_firstCandidate, &b_trig_L1_mu_firstCandidate);
   fChain->SetBranchAddress("trig_L1_mu_moreCandInRoI", &trig_L1_mu_moreCandInRoI, &b_trig_L1_mu_moreCandInRoI);
   fChain->SetBranchAddress("trig_L1_mu_moreCandInSector", &trig_L1_mu_moreCandInSector, &b_trig_L1_mu_moreCandInSector);
   fChain->SetBranchAddress("trig_L1_mu_source", &trig_L1_mu_source, &b_trig_L1_mu_source);
   fChain->SetBranchAddress("trig_L1_mu_hemisphere", &trig_L1_mu_hemisphere, &b_trig_L1_mu_hemisphere);
   fChain->SetBranchAddress("trig_L1_mu_charge", &trig_L1_mu_charge, &b_trig_L1_mu_charge);
   fChain->SetBranchAddress("trig_L1_mu_vetoed", &trig_L1_mu_vetoed, &b_trig_L1_mu_vetoed);
   fChain->SetBranchAddress("mu_n", &mu_n, &b_mu_n);
   fChain->SetBranchAddress("mu_pt", &mu_pt, &b_mu_pt);
   fChain->SetBranchAddress("mu_eta", &mu_eta, &b_mu_eta);
   fChain->SetBranchAddress("mu_phi", &mu_phi, &b_mu_phi);
   fChain->SetBranchAddress("mu_m", &mu_m, &b_mu_m);
   fChain->SetBranchAddress("mu_charge", &mu_charge, &b_mu_charge);
   fChain->SetBranchAddress("mu_author", &mu_author, &b_mu_author);
   fChain->SetBranchAddress("mu_allAuthors", &mu_allAuthors, &b_mu_allAuthors);
   fChain->SetBranchAddress("mu_muonType", &mu_muonType, &b_mu_muonType);
   fChain->SetBranchAddress("mu_etcone20", &mu_etcone20, &b_mu_etcone20);
   fChain->SetBranchAddress("mu_etcone30", &mu_etcone30, &b_mu_etcone30);
   fChain->SetBranchAddress("mu_etcone40", &mu_etcone40, &b_mu_etcone40);
   fChain->SetBranchAddress("mu_ptcone20", &mu_ptcone20, &b_mu_ptcone20);
   fChain->SetBranchAddress("mu_ptcone30", &mu_ptcone30, &b_mu_ptcone30);
   fChain->SetBranchAddress("mu_ptcone40", &mu_ptcone40, &b_mu_ptcone40);
   fChain->SetBranchAddress("mu_trackfitchi2", &mu_trackfitchi2, &b_mu_trackfitchi2);
   fChain->SetBranchAddress("mu_trackfitndof", &mu_trackfitndof, &b_mu_trackfitndof);
   fChain->SetBranchAddress("mu_isPassedMCP", &mu_isPassedMCP, &b_mu_isPassedMCP);
   fChain->SetBranchAddress("mu_quality", &mu_quality, &b_mu_quality);
   /*
   fChain->SetBranchAddress("mu_msInnerMatchChi2", &mu_msInnerMatchChi2, &b_mu_msInnerMatchChi2);
   fChain->SetBranchAddress("mu_msOuterMatchChi2", &mu_msOuterMatchChi2, &b_mu_msOuterMatchChi2);
   fChain->SetBranchAddress("mu_msInnerMatchDOF", &mu_msInnerMatchDOF, &b_mu_msInnerMatchDOF);
   fChain->SetBranchAddress("mu_msOuterMatchDOF", &mu_msOuterMatchDOF, &b_mu_msOuterMatchDOF);
   fChain->SetBranchAddress("mu_nOutliersOnTrack", &mu_nOutliersOnTrack, &b_mu_nOutliersOnTrack);
   fChain->SetBranchAddress("mu_nBLHits", &mu_nBLHits, &b_mu_nBLHits);
   fChain->SetBranchAddress("mu_nPixHits", &mu_nPixHits, &b_mu_nPixHits);
   fChain->SetBranchAddress("mu_nSCTHits", &mu_nSCTHits, &b_mu_nSCTHits);
   fChain->SetBranchAddress("mu_nTRTHits", &mu_nTRTHits, &b_mu_nTRTHits);
   fChain->SetBranchAddress("mu_nTRTHighTHits", &mu_nTRTHighTHits, &b_mu_nTRTHighTHits);
   fChain->SetBranchAddress("mu_nBLSharedHits", &mu_nBLSharedHits, &b_mu_nBLSharedHits);
   fChain->SetBranchAddress("mu_nPixSharedHits", &mu_nPixSharedHits, &b_mu_nPixSharedHits);
   fChain->SetBranchAddress("mu_nPixHoles", &mu_nPixHoles, &b_mu_nPixHoles);
   fChain->SetBranchAddress("mu_nSCTSharedHits", &mu_nSCTSharedHits, &b_mu_nSCTSharedHits);
   fChain->SetBranchAddress("mu_nSCTHoles", &mu_nSCTHoles, &b_mu_nSCTHoles);
   fChain->SetBranchAddress("mu_nTRTOutliers", &mu_nTRTOutliers, &b_mu_nTRTOutliers);
   fChain->SetBranchAddress("mu_nTRTHighTOutliers", &mu_nTRTHighTOutliers, &b_mu_nTRTHighTOutliers);
   fChain->SetBranchAddress("mu_nGangedPixels", &mu_nGangedPixels, &b_mu_nGangedPixels);
   fChain->SetBranchAddress("mu_nPixelDeadSensors", &mu_nPixelDeadSensors, &b_mu_nPixelDeadSensors);
   fChain->SetBranchAddress("mu_nSCTDeadSensors", &mu_nSCTDeadSensors, &b_mu_nSCTDeadSensors);
   fChain->SetBranchAddress("mu_nTRTDeadStraws", &mu_nTRTDeadStraws, &b_mu_nTRTDeadStraws);
   fChain->SetBranchAddress("mu_expectBLayerHit", &mu_expectBLayerHit, &b_mu_expectBLayerHit);
   fChain->SetBranchAddress("mu_nPrecisionLayers", &mu_nPrecisionLayers, &b_mu_nPrecisionLayers);
   fChain->SetBranchAddress("mu_nPrecisionHoleLayers", &mu_nPrecisionHoleLayers, &b_mu_nPrecisionHoleLayers);
   fChain->SetBranchAddress("mu_nPhiLayers", &mu_nPhiLayers, &b_mu_nPhiLayers);
   fChain->SetBranchAddress("mu_nPhiHoleLayers", &mu_nPhiHoleLayers, &b_mu_nPhiHoleLayers);
   fChain->SetBranchAddress("mu_nTrigEtaLayers", &mu_nTrigEtaLayers, &b_mu_nTrigEtaLayers);
   fChain->SetBranchAddress("mu_nTrigEtaHoleLayers", &mu_nTrigEtaHoleLayers, &b_mu_nTrigEtaHoleLayers);
   fChain->SetBranchAddress("mu_primarySector", &mu_primarySector, &b_mu_primarySector);
   fChain->SetBranchAddress("mu_secondarySector", &mu_secondarySector, &b_mu_secondarySector);
   fChain->SetBranchAddress("mu_nInnerSmallHits", &mu_nInnerSmallHits, &b_mu_nInnerSmallHits);
   fChain->SetBranchAddress("mu_nInnerLargeHits", &mu_nInnerLargeHits, &b_mu_nInnerLargeHits);
   fChain->SetBranchAddress("mu_nMiddleSmallHits", &mu_nMiddleSmallHits, &b_mu_nMiddleSmallHits);
   fChain->SetBranchAddress("mu_nMiddleLargeHits", &mu_nMiddleLargeHits, &b_mu_nMiddleLargeHits);
   fChain->SetBranchAddress("mu_nOuterSmallHits", &mu_nOuterSmallHits, &b_mu_nOuterSmallHits);
   fChain->SetBranchAddress("mu_nOuterLargeHits", &mu_nOuterLargeHits, &b_mu_nOuterLargeHits);
   fChain->SetBranchAddress("mu_nExtendedSmallHits", &mu_nExtendedSmallHits, &b_mu_nExtendedSmallHits);
   fChain->SetBranchAddress("mu_nExtendedLargeHits", &mu_nExtendedLargeHits, &b_mu_nExtendedLargeHits);
   fChain->SetBranchAddress("mu_nInnerSmallHoles", &mu_nInnerSmallHoles, &b_mu_nInnerSmallHoles);
   fChain->SetBranchAddress("mu_nInnerLargeHoles", &mu_nInnerLargeHoles, &b_mu_nInnerLargeHoles);
   fChain->SetBranchAddress("mu_nMiddleSmallHoles", &mu_nMiddleSmallHoles, &b_mu_nMiddleSmallHoles);
   fChain->SetBranchAddress("mu_nMiddleLargeHoles", &mu_nMiddleLargeHoles, &b_mu_nMiddleLargeHoles);
   fChain->SetBranchAddress("mu_nOuterSmallHoles", &mu_nOuterSmallHoles, &b_mu_nOuterSmallHoles);
   fChain->SetBranchAddress("mu_nOuterLargeHoles", &mu_nOuterLargeHoles, &b_mu_nOuterLargeHoles);
   fChain->SetBranchAddress("mu_nExtendedSmallHoles", &mu_nExtendedSmallHoles, &b_mu_nExtendedSmallHoles);
   fChain->SetBranchAddress("mu_nExtendedLargeHoles", &mu_nExtendedLargeHoles, &b_mu_nExtendedLargeHoles);
   fChain->SetBranchAddress("mu_nPhiLayer1Hits", &mu_nPhiLayer1Hits, &b_mu_nPhiLayer1Hits);
   fChain->SetBranchAddress("mu_nPhiLayer2Hits", &mu_nPhiLayer2Hits, &b_mu_nPhiLayer2Hits);
   fChain->SetBranchAddress("mu_nPhiLayer3Hits", &mu_nPhiLayer3Hits, &b_mu_nPhiLayer3Hits);
   fChain->SetBranchAddress("mu_nPhiLayer4Hits", &mu_nPhiLayer4Hits, &b_mu_nPhiLayer4Hits);
   fChain->SetBranchAddress("mu_nEtaLayer1Hits", &mu_nEtaLayer1Hits, &b_mu_nEtaLayer1Hits);
   fChain->SetBranchAddress("mu_nEtaLayer2Hits", &mu_nEtaLayer2Hits, &b_mu_nEtaLayer2Hits);
   fChain->SetBranchAddress("mu_nEtaLayer3Hits", &mu_nEtaLayer3Hits, &b_mu_nEtaLayer3Hits);
   fChain->SetBranchAddress("mu_nEtaLayer4Hits", &mu_nEtaLayer4Hits, &b_mu_nEtaLayer4Hits);
   fChain->SetBranchAddress("mu_nPhiLayer1Holes", &mu_nPhiLayer1Holes, &b_mu_nPhiLayer1Holes);
   fChain->SetBranchAddress("mu_nPhiLayer2Holes", &mu_nPhiLayer2Holes, &b_mu_nPhiLayer2Holes);
   fChain->SetBranchAddress("mu_nPhiLayer3Holes", &mu_nPhiLayer3Holes, &b_mu_nPhiLayer3Holes);
   fChain->SetBranchAddress("mu_nPhiLayer4Holes", &mu_nPhiLayer4Holes, &b_mu_nPhiLayer4Holes);
   fChain->SetBranchAddress("mu_nEtaLayer1Holes", &mu_nEtaLayer1Holes, &b_mu_nEtaLayer1Holes);
   fChain->SetBranchAddress("mu_nEtaLayer2Holes", &mu_nEtaLayer2Holes, &b_mu_nEtaLayer2Holes);
   fChain->SetBranchAddress("mu_nEtaLayer3Holes", &mu_nEtaLayer3Holes, &b_mu_nEtaLayer3Holes);
   fChain->SetBranchAddress("mu_nEtaLayer4Holes", &mu_nEtaLayer4Holes, &b_mu_nEtaLayer4Holes);
   fChain->SetBranchAddress("mu_cb_d0", &mu_cb_d0, &b_mu_cb_d0);
   fChain->SetBranchAddress("mu_cb_z0", &mu_cb_z0, &b_mu_cb_z0);
   fChain->SetBranchAddress("mu_cb_phi0", &mu_cb_phi0, &b_mu_cb_phi0);
   fChain->SetBranchAddress("mu_cb_theta", &mu_cb_theta, &b_mu_cb_theta);
   fChain->SetBranchAddress("mu_cb_qOverP", &mu_cb_qOverP, &b_mu_cb_qOverP);
   fChain->SetBranchAddress("mu_cb_vx", &mu_cb_vx, &b_mu_cb_vx);
   fChain->SetBranchAddress("mu_cb_vy", &mu_cb_vy, &b_mu_cb_vy);
   fChain->SetBranchAddress("mu_cb_vz", &mu_cb_vz, &b_mu_cb_vz);
   */
   fChain->SetBranchAddress("museg_n", &museg_n, &b_TGC_museg_n);
   fChain->SetBranchAddress("museg_x", &museg_x, &b_museg_x);
   fChain->SetBranchAddress("museg_y", &museg_y, &b_museg_y);
   fChain->SetBranchAddress("museg_z", &museg_z, &b_museg_z);
   fChain->SetBranchAddress("museg_px", &museg_px, &b_museg_px);
   fChain->SetBranchAddress("museg_py", &museg_py, &b_museg_py);
   fChain->SetBranchAddress("museg_pz", &museg_pz, &b_museg_pz);
   fChain->SetBranchAddress("museg_t0", &museg_t0, &b_museg_t0);
   fChain->SetBranchAddress("museg_t0error", &museg_t0error, &b_museg_t0error);
   fChain->SetBranchAddress("museg_chi2", &museg_chi2, &b_museg_chi2);
   fChain->SetBranchAddress("museg_ndof", &museg_ndof, &b_museg_ndof);
   fChain->SetBranchAddress("museg_sector", &museg_sector, &b_museg_sector);
   fChain->SetBranchAddress("museg_stationName", &museg_stationName, &b_museg_stationName);
   fChain->SetBranchAddress("museg_stationEta", &museg_stationEta, &b_museg_stationEta);
   fChain->SetBranchAddress("museg_author", &museg_author, &b_museg_author);
   fChain->SetBranchAddress("ext_mu_bias_n", &ext_mu_bias_n, &b_ext_mu_bias_n);
   fChain->SetBranchAddress("ext_mu_bias_type", &ext_mu_bias_type, &b_ext_mu_bias_type);
   fChain->SetBranchAddress("ext_mu_bias_index", &ext_mu_bias_index, &b_ext_mu_bias_index);
   fChain->SetBranchAddress("ext_mu_bias_size", &ext_mu_bias_size, &b_ext_mu_bias_size);
   fChain->SetBranchAddress("ext_mu_bias_targetVec", &ext_mu_bias_targetVec, &b_ext_mu_bias_targetVec);
   fChain->SetBranchAddress("ext_mu_bias_targetDistanceVec", &ext_mu_bias_targetDistanceVec, &b_ext_mu_bias_targetDistanceVec);
   fChain->SetBranchAddress("ext_mu_bias_targetEtaVec", &ext_mu_bias_targetEtaVec, &b_ext_mu_bias_targetEtaVec);
   fChain->SetBranchAddress("ext_mu_bias_targetPhiVec", &ext_mu_bias_targetPhiVec, &b_ext_mu_bias_targetPhiVec);
   fChain->SetBranchAddress("ext_mu_bias_targetDeltaEtaVec", &ext_mu_bias_targetDeltaEtaVec, &b_ext_mu_bias_targetDeltaEtaVec);
   fChain->SetBranchAddress("ext_mu_bias_targetDeltaPhiVec", &ext_mu_bias_targetDeltaPhiVec, &b_ext_mu_bias_targetDeltaPhiVec);
   fChain->SetBranchAddress("ext_mu_bias_targetPxVec", &ext_mu_bias_targetPxVec, &b_ext_mu_bias_targetPxVec);
   fChain->SetBranchAddress("ext_mu_bias_targetPyVec", &ext_mu_bias_targetPyVec, &b_ext_mu_bias_targetPyVec);
   fChain->SetBranchAddress("ext_mu_bias_targetPzVec", &ext_mu_bias_targetPzVec, &b_ext_mu_bias_targetPzVec);
   fChain->SetBranchAddress("ext_mu_ubias_n", &ext_mu_ubias_n, &b_ext_mu_ubias_n);
   fChain->SetBranchAddress("ext_mu_ubias_type", &ext_mu_ubias_type, &b_ext_mu_ubias_type);
   fChain->SetBranchAddress("ext_mu_ubias_index", &ext_mu_ubias_index, &b_ext_mu_ubias_index);
   fChain->SetBranchAddress("ext_mu_ubias_size", &ext_mu_ubias_size, &b_ext_mu_ubias_size);
   fChain->SetBranchAddress("ext_mu_ubias_targetVec", &ext_mu_ubias_targetVec, &b_ext_mu_ubias_targetVec);
   fChain->SetBranchAddress("ext_mu_ubias_targetDistanceVec", &ext_mu_ubias_targetDistanceVec, &b_ext_mu_ubias_targetDistanceVec);
   fChain->SetBranchAddress("ext_mu_ubias_targetEtaVec", &ext_mu_ubias_targetEtaVec, &b_ext_mu_ubias_targetEtaVec);
   fChain->SetBranchAddress("ext_mu_ubias_targetPhiVec", &ext_mu_ubias_targetPhiVec, &b_ext_mu_ubias_targetPhiVec);
   fChain->SetBranchAddress("ext_mu_ubias_targetDeltaEtaVec", &ext_mu_ubias_targetDeltaEtaVec, &b_ext_mu_ubias_targetDeltaEtaVec);
   fChain->SetBranchAddress("ext_mu_ubias_targetDeltaPhiVec", &ext_mu_ubias_targetDeltaPhiVec, &b_ext_mu_ubias_targetDeltaPhiVec);
   fChain->SetBranchAddress("ext_mu_ubias_targetPxVec", &ext_mu_ubias_targetPxVec, &b_ext_mu_ubias_targetPxVec);
   fChain->SetBranchAddress("ext_mu_ubias_targetPyVec", &ext_mu_ubias_targetPyVec, &b_ext_mu_ubias_targetPyVec);
   fChain->SetBranchAddress("ext_mu_ubias_targetPzVec", &ext_mu_ubias_targetPzVec, &b_ext_mu_ubias_targetPzVec);
   fChain->SetBranchAddress("trigger_info_n", &trigger_info_n, &b_trigger_info_n);
   fChain->SetBranchAddress("trigger_info_chain", &trigger_info_chain, &b_trigger_info_chain);
   fChain->SetBranchAddress("trigger_info_isPassed", &trigger_info_isPassed, &b_trigger_info_isPassed);
   fChain->SetBranchAddress("trigger_info_nTracks", &trigger_info_nTracks, &b_trigger_info_nTracks);
   fChain->SetBranchAddress("trigger_info_typeVec", &trigger_info_typeVec, &b_trigger_info_typeVec);
   fChain->SetBranchAddress("trigger_info_ptVec", &trigger_info_ptVec, &b_trigger_info_ptVec);
   fChain->SetBranchAddress("trigger_info_etaVec", &trigger_info_etaVec, &b_trigger_info_etaVec);
   fChain->SetBranchAddress("trigger_info_phiVec", &trigger_info_phiVec, &b_trigger_info_phiVec);
   /*
   fChain->SetBranchAddress("vxp_n", &vxp_n, &b_vxp_n);
   fChain->SetBranchAddress("vxp_x", &vxp_x, &b_vxp_x);
   fChain->SetBranchAddress("vxp_y", &vxp_y, &b_vxp_y);
   fChain->SetBranchAddress("vxp_z", &vxp_z, &b_vxp_z);
   fChain->SetBranchAddress("vxp_cov_x", &vxp_cov_x, &b_vxp_cov_x);
   fChain->SetBranchAddress("vxp_cov_y", &vxp_cov_y, &b_vxp_cov_y);
   fChain->SetBranchAddress("vxp_cov_z", &vxp_cov_z, &b_vxp_cov_z);
   fChain->SetBranchAddress("vxp_cov_xy", &vxp_cov_xy, &b_vxp_cov_xy);
   fChain->SetBranchAddress("vxp_cov_xz", &vxp_cov_xz, &b_vxp_cov_xz);
   fChain->SetBranchAddress("vxp_cov_yz", &vxp_cov_yz, &b_vxp_cov_yz);
   fChain->SetBranchAddress("vxp_chi2", &vxp_chi2, &b_vxp_chi2);
   fChain->SetBranchAddress("vxp_ndof", &vxp_ndof, &b_vxp_ndof);
   fChain->SetBranchAddress("vxp_nTracks", &vxp_nTracks, &b_vxp_nTracks);
   fChain->SetBranchAddress("vxp_type", &vxp_type, &b_vxp_type);
   */
   fChain->SetBranchAddress("TGC_prd_n", &TGC_prd_n, &b_TGC_prd_n);
   fChain->SetBranchAddress("TGC_prd_x", &TGC_prd_x, &b_TGC_prd_x);
   fChain->SetBranchAddress("TGC_prd_y", &TGC_prd_y, &b_TGC_prd_y);
   fChain->SetBranchAddress("TGC_prd_z", &TGC_prd_z, &b_TGC_prd_z);
   fChain->SetBranchAddress("TGC_prd_shortWidth", &TGC_prd_shortWidth, &b_TGC_prd_shortWidth);
   fChain->SetBranchAddress("TGC_prd_longWidth", &TGC_prd_longWidth, &b_TGC_prd_longWidth);
   fChain->SetBranchAddress("TGC_prd_length", &TGC_prd_length, &b_TGC_prd_length);
   fChain->SetBranchAddress("TGC_prd_isStrip", &TGC_prd_isStrip, &b_TGC_prd_isStrip);
   fChain->SetBranchAddress("TGC_prd_gasGap", &TGC_prd_gasGap, &b_TGC_prd_gasGap);
   fChain->SetBranchAddress("TGC_prd_channel", &TGC_prd_channel, &b_TGC_prd_channel);
   fChain->SetBranchAddress("TGC_prd_eta", &TGC_prd_eta, &b_TGC_prd_eta);
   fChain->SetBranchAddress("TGC_prd_phi", &TGC_prd_phi, &b_TGC_prd_phi);
   fChain->SetBranchAddress("TGC_prd_station", &TGC_prd_station, &b_TGC_prd_station);
   fChain->SetBranchAddress("TGC_prd_bunch", &TGC_prd_bunch, &b_TGC_prd_bunch);
   fChain->SetBranchAddress("RPC_prd_n", &RPC_prd_n, &b_RPC_prd_n);
   fChain->SetBranchAddress("RPC_prd_x", &RPC_prd_x, &b_RPC_prd_x);
   fChain->SetBranchAddress("RPC_prd_y", &RPC_prd_y, &b_RPC_prd_y);
   fChain->SetBranchAddress("RPC_prd_z", &RPC_prd_z, &b_RPC_prd_z);
   fChain->SetBranchAddress("RPC_prd_x2", &RPC_prd_x2, &b_RPC_prd_x2);
   fChain->SetBranchAddress("RPC_prd_y2", &RPC_prd_y2, &b_RPC_prd_y2);
   fChain->SetBranchAddress("RPC_prd_z2", &RPC_prd_z2, &b_RPC_prd_z2);
   fChain->SetBranchAddress("RPC_prd_time", &RPC_prd_time, &b_RPC_prd_time);
   fChain->SetBranchAddress("RPC_prd_triggerInfo", &RPC_prd_triggerInfo, &b_RPC_prd_triggerInfo);
   fChain->SetBranchAddress("RPC_prd_ambiguityFlag", &RPC_prd_ambiguityFlag, &b_RPC_prd_ambiguityFlag);
   fChain->SetBranchAddress("RPC_prd_measuresPhi", &RPC_prd_measuresPhi, &b_RPC_prd_measuresPhi);
   fChain->SetBranchAddress("RPC_prd_inRibs", &RPC_prd_inRibs, &b_RPC_prd_inRibs);
   fChain->SetBranchAddress("RPC_prd_station", &RPC_prd_station, &b_RPC_prd_station);
   fChain->SetBranchAddress("RPC_prd_stationEta", &RPC_prd_stationEta, &b_RPC_prd_stationEta);
   fChain->SetBranchAddress("RPC_prd_stationPhi", &RPC_prd_stationPhi, &b_RPC_prd_stationPhi);
   fChain->SetBranchAddress("RPC_prd_doubletR", &RPC_prd_doubletR, &b_RPC_prd_doubletR);
   fChain->SetBranchAddress("RPC_prd_doubletZ", &RPC_prd_doubletZ, &b_RPC_prd_doubletZ);
   fChain->SetBranchAddress("RPC_prd_stripWidth", &RPC_prd_stripWidth, &b_RPC_prd_stripWidth);
   fChain->SetBranchAddress("RPC_prd_stripLength", &RPC_prd_stripLength, &b_RPC_prd_stripLength);
   fChain->SetBranchAddress("RPC_prd_gasGap", &RPC_prd_gasGap, &b_RPC_prd_gasGap);
   fChain->SetBranchAddress("RPC_prd_channel", &RPC_prd_channel, &b_RPC_prd_channel);
   /*
   fChain->SetBranchAddress("MDT_prd_n", &MDT_prd_n, &b_MDT_prd_n);
   fChain->SetBranchAddress("MDT_prd_x", &MDT_prd_x, &b_MDT_prd_x);
   fChain->SetBranchAddress("MDT_prd_y", &MDT_prd_y, &b_MDT_prd_y);
   fChain->SetBranchAddress("MDT_prd_z", &MDT_prd_z, &b_MDT_prd_z);
   fChain->SetBranchAddress("MDT_prd_adc", &MDT_prd_adc, &b_MDT_prd_adc);
   fChain->SetBranchAddress("MDT_prd_tdc", &MDT_prd_tdc, &b_MDT_prd_tdc);
   fChain->SetBranchAddress("MDT_prd_status", &MDT_prd_status, &b_MDT_prd_status);
   fChain->SetBranchAddress("MDT_prd_drift_radius", &MDT_prd_drift_radius, &b_MDT_prd_drift_radius);
   fChain->SetBranchAddress("MDT_prd_drift_radius_error", &MDT_prd_drift_radius_error, &b_MDT_prd_drift_radius_error);
   */
   fChain->SetBranchAddress("TGC_coin_n", &TGC_coin_n, &b_TGC_coin_n);
   fChain->SetBranchAddress("TGC_coin_x_In", &TGC_coin_x_In, &b_TGC_coin_x_In);
   fChain->SetBranchAddress("TGC_coin_y_In", &TGC_coin_y_In, &b_TGC_coin_y_In);
   fChain->SetBranchAddress("TGC_coin_z_In", &TGC_coin_z_In, &b_TGC_coin_z_In);
   fChain->SetBranchAddress("TGC_coin_x_Out", &TGC_coin_x_Out, &b_TGC_coin_x_Out);
   fChain->SetBranchAddress("TGC_coin_y_Out", &TGC_coin_y_Out, &b_TGC_coin_y_Out);
   fChain->SetBranchAddress("TGC_coin_z_Out", &TGC_coin_z_Out, &b_TGC_coin_z_Out);
   fChain->SetBranchAddress("TGC_coin_width_In", &TGC_coin_width_In, &b_TGC_coin_width_In);
   fChain->SetBranchAddress("TGC_coin_width_Out", &TGC_coin_width_Out, &b_TGC_coin_width_Out);
   fChain->SetBranchAddress("TGC_coin_width_R", &TGC_coin_width_R, &b_TGC_coin_width_R);
   fChain->SetBranchAddress("TGC_coin_width_Phi", &TGC_coin_width_Phi, &b_TGC_coin_width_Phi);
   fChain->SetBranchAddress("TGC_coin_isAside", &TGC_coin_isAside, &b_TGC_coin_isAside);
   fChain->SetBranchAddress("TGC_coin_isForward", &TGC_coin_isForward, &b_TGC_coin_isForward);
   fChain->SetBranchAddress("TGC_coin_isStrip", &TGC_coin_isStrip, &b_TGC_coin_isStrip);
   fChain->SetBranchAddress("TGC_coin_isInner", &TGC_coin_isInner, &b_TGC_coin_isInner);
   fChain->SetBranchAddress("TGC_coin_isPositiveDeltaR", &TGC_coin_isPositiveDeltaR, &b_TGC_coin_isPositiveDeltaR);
   fChain->SetBranchAddress("TGC_coin_type", &TGC_coin_type, &b_TGC_coin_type);
   fChain->SetBranchAddress("TGC_coin_trackletId", &TGC_coin_trackletId, &b_TGC_coin_trackletId);
   fChain->SetBranchAddress("TGC_coin_trackletIdStrip", &TGC_coin_trackletIdStrip, &b_TGC_coin_trackletIdStrip);
   fChain->SetBranchAddress("TGC_coin_phi", &TGC_coin_phi, &b_TGC_coin_phi);
   fChain->SetBranchAddress("TGC_coin_roi", &TGC_coin_roi, &b_TGC_coin_roi);
   fChain->SetBranchAddress("TGC_coin_pt", &TGC_coin_pt, &b_TGC_coin_pt);
   fChain->SetBranchAddress("TGC_coin_delta", &TGC_coin_delta, &b_TGC_coin_delta);
   fChain->SetBranchAddress("TGC_coin_sub", &TGC_coin_sub, &b_TGC_coin_sub);
   fChain->SetBranchAddress("TGC_coin_veto", &TGC_coin_veto, &b_TGC_coin_veto);
   fChain->SetBranchAddress("TGC_coin_bunch", &TGC_coin_bunch, &b_TGC_coin_bunch);
   fChain->SetBranchAddress("TGC_coin_inner", &TGC_coin_inner, &b_TGC_coin_inner);
   fChain->SetBranchAddress("TILE_murcv_trig_n", &TILE_murcv_trig_n, &b_TILE_murcv_trig_n);
   fChain->SetBranchAddress("TILE_murcv_trig_mod", &TILE_murcv_trig_mod, &b_TILE_murcv_trig_mod);
   fChain->SetBranchAddress("TILE_murcv_trig_part", &TILE_murcv_trig_part, &b_TILE_murcv_trig_part);
   fChain->SetBranchAddress("TILE_murcv_trig_bit0", &TILE_murcv_trig_bit0, &b_TILE_murcv_trig_bit0);
   fChain->SetBranchAddress("TILE_murcv_trig_bit1", &TILE_murcv_trig_bit1, &b_TILE_murcv_trig_bit1);
   fChain->SetBranchAddress("TILE_murcv_trig_bit2", &TILE_murcv_trig_bit2, &b_TILE_murcv_trig_bit2);
   fChain->SetBranchAddress("TILE_murcv_trig_bit3", &TILE_murcv_trig_bit3, &b_TILE_murcv_trig_bit3);
   fChain->SetBranchAddress("TILE_murcv_raw_n", &TILE_murcv_raw_n, &b_TILE_murcv_raw_n);
   fChain->SetBranchAddress("TILE_murcv_raw_count", &TILE_murcv_raw_count, &b_TILE_murcv_raw_count);
   fChain->SetBranchAddress("TILE_murcv_raw_energy", &TILE_murcv_raw_energy, &b_TILE_murcv_raw_energy);
   fChain->SetBranchAddress("TILE_murcv_raw_ros", &TILE_murcv_raw_ros, &b_TILE_murcv_raw_ros);
   fChain->SetBranchAddress("TILE_murcv_raw_drawer", &TILE_murcv_raw_drawer, &b_TILE_murcv_raw_drawer);
   fChain->SetBranchAddress("TILE_murcv_raw_channel", &TILE_murcv_raw_channel, &b_TILE_murcv_raw_channel);
   fChain->SetBranchAddress("TILE_murcv_digit_n", &TILE_murcv_digit_n, &b_TILE_murcv_digit_n);
   fChain->SetBranchAddress("TILE_murcv_digit_nSamples", &TILE_murcv_digit_nSamples, &b_TILE_murcv_digit_nSamples);
   fChain->SetBranchAddress("TILE_murcv_digit_ros", &TILE_murcv_digit_ros, &b_TILE_murcv_digit_ros);
   fChain->SetBranchAddress("TILE_murcv_digit_drawer", &TILE_murcv_digit_drawer, &b_TILE_murcv_digit_drawer);
   fChain->SetBranchAddress("TILE_murcv_digit_channel", &TILE_murcv_digit_channel, &b_TILE_murcv_digit_channel);
   fChain->SetBranchAddress("TILE_murcv_digit_sampleVec", &TILE_murcv_digit_sampleVec, &b_TILE_murcv_digit_sampleVec);
   /*
   fChain->SetBranchAddress("TGC_hierarchy_n", &TGC_hierarchy_n, &b_TGC_hierarchy_n);
   fChain->SetBranchAddress("TGC_hierarchy_index", &TGC_hierarchy_index, &b_TGC_hierarchy_index);
   fChain->SetBranchAddress("TGC_hierarchy_dR_hiPt", &TGC_hierarchy_dR_hiPt, &b_TGC_hierarchy_dR_hiPt);
   fChain->SetBranchAddress("TGC_hierarchy_dPhi_hiPt", &TGC_hierarchy_dPhi_hiPt, &b_TGC_hierarchy_dPhi_hiPt);
   fChain->SetBranchAddress("TGC_hierarchy_dR_tracklet", &TGC_hierarchy_dR_tracklet, &b_TGC_hierarchy_dR_tracklet);
   fChain->SetBranchAddress("TGC_hierarchy_dPhi_tracklet", &TGC_hierarchy_dPhi_tracklet, &b_TGC_hierarchy_dPhi_tracklet);
   fChain->SetBranchAddress("TGC_hierarchy_isChamberBoundary", &TGC_hierarchy_isChamberBoundary, &b_TGC_hierarchy_isChamberBoundary);
   */
   fChain->SetBranchAddress("muctpi_candidateMultiplicities", &muctpi_candidateMultiplicities, &b_muctpi_candidateMultiplicities);
   fChain->SetBranchAddress("muctpi_nDataWords", &muctpi_nDataWords, &b_muctpi_nDataWords);
   fChain->SetBranchAddress("muctpi_dataWords", &muctpi_dataWords, &b_muctpi_dataWords);
   fChain->SetBranchAddress("muctpi_dw_eta", &muctpi_dw_eta, &b_muctpi_dw_eta);
   fChain->SetBranchAddress("muctpi_dw_phi", &muctpi_dw_phi, &b_muctpi_dw_phi);
   fChain->SetBranchAddress("muctpi_dw_source", &muctpi_dw_source, &b_muctpi_dw_source);
   fChain->SetBranchAddress("muctpi_dw_hemisphere", &muctpi_dw_hemisphere, &b_muctpi_dw_hemisphere);
   fChain->SetBranchAddress("muctpi_dw_bcid", &muctpi_dw_bcid, &b_muctpi_dw_bcid);
   fChain->SetBranchAddress("muctpi_dw_sectorID", &muctpi_dw_sectorID, &b_muctpi_dw_sectorID);
   fChain->SetBranchAddress("muctpi_dw_thrNumber", &muctpi_dw_thrNumber, &b_muctpi_dw_thrNumber);
   fChain->SetBranchAddress("muctpi_dw_roi", &muctpi_dw_roi, &b_muctpi_dw_roi);
   fChain->SetBranchAddress("muctpi_dw_veto", &muctpi_dw_veto, &b_muctpi_dw_veto);
   fChain->SetBranchAddress("muctpi_dw_firstCandidate", &muctpi_dw_firstCandidate, &b_muctpi_dw_firstCandidate);
   fChain->SetBranchAddress("muctpi_dw_moreCandInRoI", &muctpi_dw_moreCandInRoI, &b_muctpi_dw_moreCandInRoI);
   fChain->SetBranchAddress("muctpi_dw_moreCandInSector", &muctpi_dw_moreCandInSector, &b_muctpi_dw_moreCandInSector);
   fChain->SetBranchAddress("muctpi_dw_charge", &muctpi_dw_charge, &b_muctpi_dw_charge);
   fChain->SetBranchAddress("muctpi_dw_candidateVetoed", &muctpi_dw_candidateVetoed, &b_muctpi_dw_candidateVetoed);
   /*
   fChain->SetBranchAddress("TILE_cell_n", &TILE_cell_n, &b_TILE_cell_n);
   fChain->SetBranchAddress("TILE_cell_E", &TILE_cell_E, &b_TILE_cell_E);
   fChain->SetBranchAddress("TILE_cell_ene1", &TILE_cell_ene1, &b_TILE_cell_ene1);
   fChain->SetBranchAddress("TILE_cell_eta", &TILE_cell_eta, &b_TILE_cell_eta);
   fChain->SetBranchAddress("TILE_cell_phi", &TILE_cell_phi, &b_TILE_cell_phi);
   fChain->SetBranchAddress("TILE_cell_sinTh", &TILE_cell_sinTh, &b_TILE_cell_sinTh);
   fChain->SetBranchAddress("TILE_cell_cosTh", &TILE_cell_cosTh, &b_TILE_cell_cosTh);
   fChain->SetBranchAddress("TILE_cell_cotTh", &TILE_cell_cotTh, &b_TILE_cell_cotTh);
   fChain->SetBranchAddress("TILE_cell_x", &TILE_cell_x, &b_TILE_cell_x);
   fChain->SetBranchAddress("TILE_cell_y", &TILE_cell_y, &b_TILE_cell_y);
   fChain->SetBranchAddress("TILE_cell_z", &TILE_cell_z, &b_TILE_cell_z);
   fChain->SetBranchAddress("TILE_cell_badcell", &TILE_cell_badcell, &b_TILE_cell_badcell);
   fChain->SetBranchAddress("TILE_cell_partition", &TILE_cell_partition, &b_TILE_cell_partition);
   fChain->SetBranchAddress("TILE_cell_section", &TILE_cell_section, &b_TILE_cell_section);
   fChain->SetBranchAddress("TILE_cell_side", &TILE_cell_side, &b_TILE_cell_side);
   fChain->SetBranchAddress("TILE_cell_module", &TILE_cell_module, &b_TILE_cell_module);
   fChain->SetBranchAddress("TILE_cell_tower", &TILE_cell_tower, &b_TILE_cell_tower);
   fChain->SetBranchAddress("TILE_cell_sample", &TILE_cell_sample, &b_TILE_cell_sample);
   fChain->SetBranchAddress("mc_n",&mc_n,&b_mc_n);
   fChain->SetBranchAddress("mc_pt",&mc_pt,&b_mc_pt);
   fChain->SetBranchAddress("mc_eta",&mc_eta,&b_mc_eta);
   fChain->SetBranchAddress("mc_phi",&mc_phi,&b_mc_phi);
   fChain->SetBranchAddress("mc_m",&mc_m,&b_mc_m);
   fChain->SetBranchAddress("mc_charge",&mc_charge,&b_mc_charge);
   */
   Notify();
}

Bool_t physics::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void physics::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t physics::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef physics_cxx
