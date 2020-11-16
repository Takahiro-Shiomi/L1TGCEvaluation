#define physics_cxx
#include "physics.h"

#include "TGCRPhiCoincidenceMap.h"
#include "HistogramManager.h"
#include "RoIObj.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TLorentzVector.h>

void physics::Loop(std::string BWCW_path)
{

    //////////////////////
    //// INITIALIZE
    ////////////////////
    std::string hsName="CW_MC_Tau_Data_1105";
    HistogramManager hs(hsName);

    enum {
        TYPE_TRACKLET,
        TYPE_HIPT,
        TYPE_SL,
        TYPE_UNKNOWN,
        TYPE_TRACKLET_EIF};

    enum {NumberOfSide=2, NumberOfOctant=8};
    TGCRPhiCoincidenceMap *bwMap[NumberOfSide][NumberOfOctant];

    if(BWCW_path!=""){
        for(int i_side=0;i_side!=NumberOfSide;i_side++){
            for(int j_octant=0;j_octant!=NumberOfOctant;j_octant++){
                bwMap[i_side][j_octant] = new TGCRPhiCoincidenceMap(BWCW_path,i_side,j_octant);
            }
        }
    }
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;

    //////////////////////
    ////  EVENT LOOP
    /////////////////////
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //for (Long64_t jentry=0; jentry<2000000;jentry++) {
    //for (Long64_t jentry=2000000; jentry<4000000;jentry++) {
    //for (Long64_t jentry=4000000; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        if((int)jentry % 100000 == 0){std::cout<<"JOBS....."<<(int)jentry / 10000<<"%....COMPLETE"<<std::endl;}

        hs.Clear();

        ///Fill RoI Info
        std::vector<RoIObj> RoIs;
        if(TGC_coin_n>=1){
            bool RoIFlag = RoISL(RoIs);
            ///Fill DR/DPhi Info
            if(RoIFlag){ RoIDRDPhi(RoIs); }

            //2 to 4 candidates algorism
            if(RoIFlag){ RoITwoToFour(RoIs); }
            
            //Strip 2-station algorism
            RoIHPTLPT(RoIs);

            ///calc pT
            if(RoIs.size()!=0){
            for(auto &roi : RoIs){
                int pT = bwMap[roi.getSide()][roi.getOctant()]->test_Run3(roi.getOctant(), roi.getModuleId(),roi.getRoI(),roi.getCoinType(),roi.getDR(),roi.getDPhi());
                roi.setPt(std::abs(pT));
                if(pT<0){roi.setCharge(0);} //Charge = -
                else if(pT>0){roi.setCharge(1);} //Charge = +
                else {roi.setCharge(2);} //Charge = unknown
            }
            }
        }
        ///Fill Ntuple
        FillNtuple(hs,RoIs);
    }
}

void physics::FillNtuple(HistogramManager &hs, std::vector<RoIObj>& RoIs)
{
    //EventInformation
    hs.eventnumber = EventNumber;
    hs.bcid = bcid;

    //Offline
    hs.muon_n = mu_n;
    for(int i=0;i!=mu_n;i++){
        hs.muon_pt->push_back((*mu_pt)[i]);
        hs.muon_eta->push_back((*mu_eta)[i]);
        hs.muon_phi->push_back((*mu_phi)[i]);
        hs.muon_m->push_back((*mu_m)[i]);
        hs.muon_charge->push_back((*mu_charge)[i]);
        hs.muon_author->push_back((*mu_author)[i]);
        hs.muon_Type->push_back((*mu_muonType)[i]);
    }

    //Museg
    hs.museg_n = museg_n;
    for(int i=0;i!=museg_n;i++){
        hs.museg_x->push_back((*museg_x)[i]);
        hs.museg_y->push_back((*museg_y)[i]);
        hs.museg_z->push_back((*museg_z)[i]);
        hs.museg_px->push_back((*museg_px)[i]);
        hs.museg_py->push_back((*museg_py)[i]);
        hs.museg_pz->push_back((*museg_pz)[i]);
        hs.museg_t0->push_back((*museg_t0)[i]);
        hs.museg_t0error->push_back((*museg_t0error)[i]);
        hs.museg_chi2->push_back((*museg_chi2)[i]);
        hs.museg_ndof->push_back((*museg_ndof)[i]);
        hs.museg_sector->push_back((*museg_sector)[i]);
        hs.museg_stationName->push_back((*museg_stationName)[i]);
        hs.museg_stationEta->push_back((*museg_stationEta)[i]);
        hs.museg_author->push_back((*museg_author)[i]);
    }

    //Extrapolate
    hs.ext_mu_n = ext_mu_bias_n;
    for(int i=0;i!=ext_mu_bias_n;i++){
        hs.ext_mu_type->push_back((*ext_mu_bias_type)[i]);
        hs.ext_mu_index->push_back((*ext_mu_bias_index)[i]);
        hs.ext_mu_size->push_back((*ext_mu_bias_size)[i]);
        std::vector<int> targetVec;
        std::vector<float> DistanceVec;
        std::vector<float> EtaVec;
        std::vector<float> PhiVec;
        std::vector<float> DeltaEtaVec;
        std::vector<float> DeltaPhiVec;
        std::vector<float> PxVec;
        std::vector<float> PyVec;
        std::vector<float> PzVec;
        for(int j=0;j!=ext_mu_bias_targetVec->at(i).size();j++){
            targetVec.push_back(ext_mu_bias_targetVec->at(i).at(j));
            DistanceVec.push_back(ext_mu_bias_targetDistanceVec->at(i).at(j));
            EtaVec.push_back(ext_mu_bias_targetEtaVec->at(i).at(j));
            PhiVec.push_back(ext_mu_bias_targetPhiVec->at(i).at(j));
            DeltaEtaVec.push_back(ext_mu_bias_targetDeltaEtaVec->at(i).at(j));
            DeltaPhiVec.push_back(ext_mu_bias_targetDeltaPhiVec->at(i).at(j));
            PxVec.push_back(ext_mu_bias_targetPxVec->at(i).at(j));
            PyVec.push_back(ext_mu_bias_targetPyVec->at(i).at(j));
            PzVec.push_back(ext_mu_bias_targetPzVec->at(i).at(j));
        }
        hs.ext_mu_targetVec->push_back(targetVec);
        hs.ext_mu_targetDistanceVec->push_back(DistanceVec);
        hs.ext_mu_targetEtaVec->push_back(EtaVec);
        hs.ext_mu_targetPhiVec->push_back(PhiVec);
        hs.ext_mu_targetDeltaPhiVec->push_back(DeltaPhiVec);
        hs.ext_mu_targetDeltaEtaVec->push_back(DeltaEtaVec);
        hs.ext_mu_targetPxVec->push_back(PxVec);
        hs.ext_mu_targetPyVec->push_back(PyVec);
        hs.ext_mu_targetPzVec->push_back(PzVec);
    }

    //TGCPRD
    hs.TGC_prd_n = TGC_prd_n;
    for(int i=0;i!=TGC_prd_n;i++){
        hs.TGC_prd_x->push_back((*TGC_prd_x)[i]);
        hs.TGC_prd_y->push_back((*TGC_prd_y)[i]);
        hs.TGC_prd_z->push_back((*TGC_prd_z)[i]);
        hs.TGC_prd_shortWidth->push_back((*TGC_prd_shortWidth)[i]);
        hs.TGC_prd_longWidth->push_back((*TGC_prd_longWidth)[i]);
        hs.TGC_prd_length->push_back((*TGC_prd_length)[i]);
        hs.TGC_prd_isStrip->push_back((*TGC_prd_isStrip)[i]);
        hs.TGC_prd_gasGap->push_back((*TGC_prd_gasGap)[i]);
        hs.TGC_prd_channel->push_back((*TGC_prd_channel)[i]);
        hs.TGC_prd_eta->push_back((*TGC_prd_eta)[i]);
        hs.TGC_prd_phi->push_back((*TGC_prd_phi)[i]);
        hs.TGC_prd_station->push_back((*TGC_prd_station)[i]);
        hs.TGC_prd_bunch->push_back((*TGC_prd_bunch)[i]);
    }

    //RPCPRD
    hs.RPC_prd_n = RPC_prd_n;
    for(int i=0;i!=RPC_prd_n;i++){
        hs.RPC_prd_x->push_back((*RPC_prd_x)[i]);
        hs.RPC_prd_y->push_back((*RPC_prd_y)[i]);
        hs.RPC_prd_z->push_back((*RPC_prd_z)[i]);
        hs.RPC_prd_x2->push_back((*RPC_prd_x2)[i]);
        hs.RPC_prd_y2->push_back((*RPC_prd_y2)[i]);
        hs.RPC_prd_z2->push_back((*RPC_prd_z2)[i]);
        hs.RPC_prd_triggerInfo->push_back((*RPC_prd_triggerInfo)[i]);
        hs.RPC_prd_ambiguityFlag->push_back((*RPC_prd_ambiguityFlag)[i]);
        hs.RPC_prd_measuresPhi->push_back((*RPC_prd_measuresPhi)[i]);
        hs.RPC_prd_inRibs->push_back((*RPC_prd_inRibs)[i]);
        hs.RPC_prd_station->push_back((*RPC_prd_station)[i]);
        hs.RPC_prd_stationEta->push_back((*RPC_prd_stationEta)[i]);
        hs.RPC_prd_stationPhi->push_back((*RPC_prd_stationPhi)[i]);
        hs.RPC_prd_doubletR->push_back((*RPC_prd_doubletR)[i]);
        hs.RPC_prd_doubletZ->push_back((*RPC_prd_doubletZ)[i]);
        hs.RPC_prd_stripWidth->push_back((*RPC_prd_stripWidth)[i]);
        hs.RPC_prd_stripLength->push_back((*RPC_prd_stripLength)[i]);
        hs.RPC_prd_gasGap->push_back((*RPC_prd_gasGap)[i]);
        hs.RPC_prd_channel->push_back((*RPC_prd_channel)[i]);
    }

    //Tile
    hs.TILE_murcv_trig_n = TILE_murcv_trig_n;
    for(int i=0;i!=TILE_murcv_trig_n;i++){
        hs.TILE_murcv_trig_mod->push_back((*TILE_murcv_trig_mod)[i]);
        hs.TILE_murcv_trig_part->push_back((*TILE_murcv_trig_part)[i]);
        hs.TILE_murcv_trig_bit0->push_back((*TILE_murcv_trig_bit0)[i]);
        hs.TILE_murcv_trig_bit1->push_back((*TILE_murcv_trig_bit1)[i]);
        hs.TILE_murcv_trig_bit2->push_back((*TILE_murcv_trig_bit2)[i]);
        hs.TILE_murcv_trig_bit3->push_back((*TILE_murcv_trig_bit3)[i]);
    }
    hs.TILE_murcv_raw_n = TILE_murcv_raw_n;
    for(int i=0;i!=TILE_murcv_raw_n;i++){
        hs.TILE_murcv_raw_count->push_back((*TILE_murcv_raw_count)[i]);
        hs.TILE_murcv_raw_energy->push_back((*TILE_murcv_raw_energy)[i]);
        hs.TILE_murcv_raw_ros->push_back((*TILE_murcv_raw_ros)[i]);
        hs.TILE_murcv_raw_drawer->push_back((*TILE_murcv_raw_drawer)[i]);
        hs.TILE_murcv_raw_channel->push_back((*TILE_murcv_raw_channel)[i]);
    }
    hs.TILE_murcv_digit_n = TILE_murcv_digit_n;
    for(int i=0;i!=TILE_murcv_digit_n;i++){
        hs.TILE_murcv_digit_nSamples->push_back((*TILE_murcv_digit_nSamples)[i]);
        hs.TILE_murcv_digit_ros->push_back((*TILE_murcv_digit_ros)[i]);
        hs.TILE_murcv_digit_drawer->push_back((*TILE_murcv_digit_drawer)[i]);
        hs.TILE_murcv_digit_channel->push_back((*TILE_murcv_digit_channel)[i]);
    }

    //HLT
    hs.HLT_info_n = trigger_info_n;
    for(int i=0;i!=trigger_info_n;i++){
        hs.HLT_info_chain->push_back((*trigger_info_chain)[i]);
        hs.HLT_info_isPassed->push_back((*trigger_info_isPassed)[i]);
        std::vector<int> typeVec;
        std::vector<float> ptVec;
        std::vector<float> phiVec;
        std::vector<float> etaVec;
        for(int j=0;j!=trigger_info_typeVec->at(i).size();j++){
            typeVec.push_back(trigger_info_typeVec->at(i).at(j));
            ptVec.push_back(trigger_info_ptVec->at(i).at(j));
            phiVec.push_back(trigger_info_phiVec->at(i).at(j));
            etaVec.push_back(trigger_info_etaVec->at(i).at(j));
        }
        hs.HLT_info_typeVec->push_back(typeVec);
        hs.HLT_info_ptVec->push_back(ptVec);
        hs.HLT_info_etaVec->push_back(etaVec);
        hs.HLT_info_phiVec->push_back(phiVec);
    }

    //RoI
    hs.muctpi_ndatawords = muctpi_nDataWords;
    for(int i=0;i!=muctpi_nDataWords;i++){
        hs.muctpi_eta->push_back((*muctpi_dw_eta)[i]);
        hs.muctpi_phi->push_back((*muctpi_dw_phi)[i]);
        hs.muctpi_source->push_back((*muctpi_dw_source)[i]);
        hs.muctpi_hemisphere->push_back((*muctpi_dw_hemisphere)[i]);
        hs.muctpi_bcid->push_back((*muctpi_dw_bcid)[i]);
        hs.muctpi_sectorID->push_back((*muctpi_dw_sectorID)[i]);
        hs.muctpi_thrNumber->push_back((*muctpi_dw_thrNumber)[i]);
        hs.muctpi_roi->push_back((*muctpi_dw_roi)[i]);
        hs.muctpi_veto->push_back((*muctpi_dw_veto)[i]);
        hs.muctpi_charge->push_back((*muctpi_dw_charge)[i]);
        hs.muctpi_candidateVetoed->push_back((*muctpi_dw_candidateVetoed)[i]);
    }

    //TGC Coin Data
    hs.tgc_coin_n = TGC_coin_n;
    for(int i=0;i!=TGC_coin_n;i++){
        hs.tgc_coin_x_In->push_back((*TGC_coin_x_In)[i]);
        hs.tgc_coin_y_In->push_back((*TGC_coin_y_In)[i]);
        hs.tgc_coin_z_In->push_back((*TGC_coin_z_In)[i]);
        hs.tgc_coin_x_Out->push_back((*TGC_coin_x_Out)[i]);
        hs.tgc_coin_y_Out->push_back((*TGC_coin_y_Out)[i]);
        hs.tgc_coin_z_Out->push_back((*TGC_coin_z_Out)[i]);
        hs.tgc_coin_width_In->push_back((*TGC_coin_width_In)[i]);
        hs.tgc_coin_width_Out->push_back((*TGC_coin_width_Out)[i]);
        hs.tgc_coin_width_R->push_back((*TGC_coin_width_R)[i]);
        hs.tgc_coin_width_Phi->push_back((*TGC_coin_width_Phi)[i]);
        hs.tgc_coin_isAside->push_back((*TGC_coin_isAside)[i]);
        hs.tgc_coin_isStrip->push_back((*TGC_coin_isStrip)[i]);
        hs.tgc_coin_isForward->push_back((*TGC_coin_isForward)[i]);
        hs.tgc_coin_isInner->push_back((*TGC_coin_isInner)[i]);
        hs.tgc_coin_type->push_back((*TGC_coin_type)[i]);
        hs.tgc_coin_trackletId->push_back((*TGC_coin_trackletId)[i]);
        hs.tgc_coin_trackletIdStrip->push_back((*TGC_coin_trackletIdStrip)[i]);
        hs.tgc_coin_phi->push_back((*TGC_coin_phi)[i]);
        hs.tgc_coin_roi->push_back((*TGC_coin_roi)[i]);
        hs.tgc_coin_pt->push_back((*TGC_coin_pt)[i]);
        hs.tgc_coin_delta->push_back((*TGC_coin_delta)[i]);
        hs.tgc_coin_sub->push_back((*TGC_coin_sub)[i]);
        hs.tgc_coin_veto->push_back((*TGC_coin_veto)[i]);
        hs.tgc_coin_bunch->push_back((*TGC_coin_bunch)[i]);
        hs.tgc_coin_inner->push_back((*TGC_coin_inner)[i]);
    }

    //Run-3
    hs.TGC_Run3_n = RoIs.size(); 
    for(auto &roi : RoIs){
        hs.TGC_Run3_type->push_back(roi.getType());
        hs.TGC_Run3_pt->push_back(roi.getPt());
        hs.TGC_Run3_station->push_back(roi.getCoinType());
        hs.TGC_Run3_DR->push_back(roi.getDR());
        hs.TGC_Run3_DPhi->push_back(roi.getDPhi());
        hs.TGC_Run3_TypeDPhi->push_back(roi.getTypeDPhi());
        hs.TGC_Run3_TypeDR->push_back(roi.getTypeDR());
        hs.TGC_Run3_Side->push_back(roi.getSide());
        hs.TGC_Run3_RoI->push_back(roi.getRoI());
        hs.TGC_Run3_PhiSector->push_back(roi.getPhiSector());
        hs.TGC_Run3_IsEndcap->push_back(roi.getIsEndcap());
        hs.TGC_Run3_TrackletIdWire->push_back(roi.getTrackletIdWire());
        hs.TGC_Run3_TrackletIdStrip->push_back(roi.getTrackletIdStrip());
        hs.TGC_Run3_x->push_back(roi.getX());
        hs.TGC_Run3_y->push_back(roi.getY());
        hs.TGC_Run3_z->push_back(roi.getZ());
        hs.TGC_Run3_R->push_back(roi.getR());
        hs.TGC_Run3_Phi->push_back(roi.getPhi());
        hs.TGC_Run3_Charge->push_back(roi.getCharge());
    }

    //Truth
    /*
       hs.mc_n = mc_n;
       for(int i=0;i!=mc_n;i++){
       hs.mc_pt->push_back((*mc_pt)[i]);
       hs.mc_eta->push_back((*mc_eta)[i]);
       hs.mc_phi->push_back((*mc_phi)[i]);
       hs.mc_m->push_back((*mc_m)[i]);
       hs.mc_charge->push_back((*mc_charge)[i]);
       }
       */ 
    hs.m_tree->Fill();
}

