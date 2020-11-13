#include "physics.h"
#include "RoIObj.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <TVector3.h>

Bool_t physics::RoISL(std::vector<RoIObj>& RoIs)
{
    int Run3_n=0;
    for(int i=0;i!=TGC_coin_n;i++){
        if((*TGC_coin_type)[i]!=2){continue;}
        if((*TGC_coin_bunch)[i]!=2){continue;}//1 = previous, 2 = current, 3 = next 
        RoIObj roi_tmp;
        roi_tmp.setRoI((*TGC_coin_roi)[i]);
        roi_tmp.setType((*TGC_coin_type)[i]);
        roi_tmp.setIsEndcap(!(*TGC_coin_isForward)[i]);
        roi_tmp.setSide((*TGC_coin_isAside)[i]);
        roi_tmp.setTrackletIdWire((*TGC_coin_trackletId)[i]);
        roi_tmp.setTrackletIdStrip((*TGC_coin_trackletIdStrip)[i]);
        roi_tmp.setX((*TGC_coin_x_Out)[i]);
        roi_tmp.setY((*TGC_coin_y_Out)[i]);
        roi_tmp.setZ((*TGC_coin_z_Out)[i]);
        roi_tmp.setR((*TGC_coin_width_R)[i]);
        roi_tmp.setPhi((*TGC_coin_width_Phi)[i]);
        roi_tmp.setPhiSector((*TGC_coin_phi)[i]);
        RoIs.push_back(roi_tmp);
    }
    Run3_n=RoIs.size();
    
    if(Run3_n>0){return true;}
    else{return false;}
}

void physics::RoIDRDPhi(std::vector<RoIObj>& RoIs)
{
    for(auto &roi : RoIs){
        for(int i=0;i!=TGC_coin_n;i++){
            if((*TGC_coin_type)[i]==2 || (*TGC_coin_type)[i]==3 || (*TGC_coin_type)[i]==4){continue;}
            if((*TGC_coin_bunch)[i]!=2){continue;}//1 = previous, 2 = current, 3 = next 
            if((*TGC_coin_phi)[i] != roi.getPhiSector()){continue;}

            if((*TGC_coin_trackletId)[i] == roi.getTrackletIdStrip() && (*TGC_coin_isStrip)[i]){
                if(roi.getTypeDPhi() <= (*TGC_coin_type)[i]){
                    if((*TGC_coin_type)[i]==1){if((*TGC_coin_width_Out)[i]!=roi.getPhi()){continue;}}
                    roi.setTypeDPhi((*TGC_coin_type)[i]);
                    roi.setDPhi((*TGC_coin_delta)[i]);
                }
            } // dPhi
            if((*TGC_coin_trackletId)[i] == roi.getTrackletIdWire()  && !(*TGC_coin_isStrip)[i]){
                if(roi.getTypeDR() <= (*TGC_coin_type)[i]){
                    if((*TGC_coin_type)[i]==1){if((*TGC_coin_width_Out)[i]!=roi.getR()){continue;}}
                    roi.setTypeDR((*TGC_coin_type)[i]);
                    roi.setDR((*TGC_coin_delta)[i]);
                }
            } //DR 
        }
    }
}

void physics::RoIHPTLPT(std::vector<RoIObj>& RoIs)
{
    std::vector<int> WireID; WireID.clear();
    std::vector<int> StripID; StripID.clear();
    std::vector<int> Wtype; Wtype.clear();
    std::vector<int> Stype; Stype.clear();
    std::vector<int> Wsec; Wsec.clear();
    std::vector<int> Ssec; Ssec.clear();
    std::vector<bool> Wsource; Wsource.clear();
    std::vector<bool> Ssource; Ssource.clear();
    std::vector<bool> Wside; Wside.clear();
    std::vector<bool> Sside; Sside.clear();
    std::vector<int> RoI; RoI.clear();

    for(int i=0;i!=TGC_coin_n;i++){
        if((*TGC_coin_bunch)[i]!=2){continue;} 
        if((*TGC_coin_type)[i]==2){break;}
        if((*TGC_coin_type)[i]==0 || (*TGC_coin_type)[i]==1){
            if((*TGC_coin_width_Out)[i]==0){continue;}
            if(!(*TGC_coin_isStrip)[i]){
                int wsec = (*TGC_coin_phi)[i];
                int wtype;
                bool wsource;
                bool wside;
                wtype = (*TGC_coin_type)[i];
                if((*TGC_coin_isForward)[i]){wsource=false;}
                else{wsource=true;}
                if((*TGC_coin_isAside)[i]){wside=true;}
                else{wside=false;}

                bool samesec = false;
                if(RoIs.size()>=1){
                    for(auto &roi : RoIs){
                        if((wsec == roi.getPhiSector())&& (wsource == roi.getIsEndcap())&& (wside == roi.getSide())){ samesec = true; } 
                    }
                }
                if(!samesec){
        
    for(int j=0;j!=i && j!=TGC_coin_n;j++){
        if((*TGC_coin_bunch)[j]!=2){continue;} 
        if((*TGC_coin_type)[j]==2 || (*TGC_coin_type)[j]==1){break;}
        if((*TGC_coin_type)[j]==0){
            if((*TGC_coin_width_Out)[j]==0){continue;}
            if((*TGC_coin_isStrip)[j]){ 
                int ssec = (*TGC_coin_phi)[j];
                int stype;
                bool ssource;
                bool sside;
                stype = (*TGC_coin_type)[j];
                if((*TGC_coin_isForward)[j]){ssource=false;}
                else{ssource=true;}
                if((*TGC_coin_isAside)[j]){sside=true;}
                else{sside=false;}

                bool samesec2 = false;
                if(RoIs.size()>=1){
                    for(auto &roi : RoIs){ 
                        if((ssec == roi.getPhiSector())&& (ssource == roi.getIsEndcap())&& (sside == roi.getSide())){ samesec2 = true; } 
                    } 
                }
                if(!samesec2){

                if((wsource == ssource)&&(wside == sside)&&(wsec == ssec)){
                    bool samecan = false;
                    if(WireID.size()>=1){
                        for(int k=0;k!=WireID.size();k++){
                            if((wsource == Wsource.at(k)) && (wside == Wside.at(k)) && (wsec == Wsec.at(k))){
                                if(wtype>=Wtype.at(k) && stype>=Stype.at(k)){
                                    int roi = getRoI(i,j);
                                    if(roi!=-1){
                                        WireID.at(k) = i;
                                        StripID.at(k) = j;
                                        Wtype.at(k) = wtype;
                                        Stype.at(k) = stype;
                                        RoI.at(k) = roi;
                                    }
                                }
                                samecan = true;
                                break;
                            }
                        }
                    }
                    if(WireID.size()==0 || !samecan){
                        int roi = getRoI(i,j);
                        if(roi!=-1){
                            WireID.push_back(i);
                            StripID.push_back(j);
                            Wtype.push_back(wtype);
                            Stype.push_back(stype);
                            Wsec.push_back(wsec);
                            Ssec.push_back(ssec);
                            Wside.push_back(wside);
                            Sside.push_back(sside);
                            Wsource.push_back(wsource);
                            Ssource.push_back(ssource);
                            RoI.push_back(roi);
                        }
                    }
                }
                }
            }
        }
    }}}}}

    if(WireID.size()!=0){
    for(int i=0;i!=WireID.size();i++){
        RoIObj roi_tmp;
        roi_tmp.setRoI(RoI.at(i));
        roi_tmp.setType(2);
        roi_tmp.setIsEndcap(!(*TGC_coin_isForward)[WireID.at(i)]);
        roi_tmp.setSide((*TGC_coin_isAside)[WireID.at(i)]);
        roi_tmp.setTrackletIdWire((*TGC_coin_trackletId)[WireID.at(i)]);
        roi_tmp.setTrackletIdStrip((*TGC_coin_trackletId)[StripID.at(i)]);
        roi_tmp.setX( ((*TGC_coin_x_Out)[WireID.at(i)] + (*TGC_coin_x_Out)[StripID.at(i)])/2 );
        roi_tmp.setY( ((*TGC_coin_y_Out)[WireID.at(i)] + (*TGC_coin_y_Out)[StripID.at(i)])/2 );
        roi_tmp.setZ( ((*TGC_coin_z_Out)[WireID.at(i)] + (*TGC_coin_z_Out)[StripID.at(i)])/2 );
        roi_tmp.setR((*TGC_coin_width_Out)[WireID.at(i)]);
        roi_tmp.setPhi((*TGC_coin_width_Out)[StripID.at(i)]);
        roi_tmp.setPhiSector((*TGC_coin_phi)[WireID.at(i)]);
        roi_tmp.setTypeDR((*TGC_coin_type)[WireID.at(i)]);
        roi_tmp.setDR((*TGC_coin_delta)[WireID.at(i)]);
        roi_tmp.setTypeDPhi((*TGC_coin_type)[StripID.at(i)]);
        roi_tmp.setDPhi((*TGC_coin_delta)[StripID.at(i)]);
        RoIs.push_back(roi_tmp);
    }
    }
}

Int_t physics::getRoI(int wid, int sid)
{
    int RoI = -1;

    float x = ((*TGC_coin_x_Out)[wid] + (*TGC_coin_x_Out)[sid])/2;
    float y = ((*TGC_coin_y_Out)[wid] + (*TGC_coin_y_Out)[sid])/2;
    float z = ((*TGC_coin_z_Out)[wid] + (*TGC_coin_z_Out)[sid])/2;

    TVector3 v1;
    v1.SetXYZ(x,y,z);
    float eta=v1.PseudoRapidity();
    float phi=v1.Phi();
    //std::cout<<"EventNumber="<<EventNumber<<"___eta="<<eta<<"___phi="<<phi<<std::endl;

    if(abs(eta)<=2.75 && abs(eta)>=1.05){
    if(abs(phi)<=3.14159265){
        std::ifstream fin("/gpfs/fs7001/shiomi/ATLAS/Residual/Ntuple/L1TGCNtupleRun3InStrip2Station/L1tgcevaluation/share/TGCRoI_Mapping.dat");

        int region;
        int isCside;
        int trigSector;
        int RoINumber;
        double etamin, etamax, phimin, phimax;
        while(fin >> isCside >> region >> trigSector >> RoINumber >> etamin >> etamax >> phimin >> phimax) {
            if(isCside != !(*TGC_coin_isAside)[wid])continue;
            if(region != (*TGC_coin_isForward)[wid])continue;

            if(eta>=etamin && eta<=etamax){
                if(phi>=phimin && phi<=phimax){
                //std::cout<<"eta="<<etamin<<","<<etamax<<"___phi="<<phimin<<","<<phimax<<"___RoI="<<RoINumber<<std::endl;
                    RoI = RoINumber;
                    break;
                }
            }
        }
        fin.close();
    }
    }
    return RoI;
}

void physics::RoITwoToFour(std::vector<RoIObj>& RoIs){

    if(RoIs.size()>=2){
        bool R_same = false;
        std::vector<int> CanNo; CanNo.clear();
        for(int can1=0;can1!=(int)RoIs.size();can1++){
            for(int can2=0;can2!=(int)RoIs.size();can2++){
                if(can1 == can2){continue;}
                auto roi = RoIs.at(can1);
                auto roi2 = RoIs.at(can2);
                if( roi.getSide() != roi2.getSide() ){continue;}
                if( roi.getIsEndcap() != roi2.getIsEndcap() ){continue;}
                if( roi.getPhiSector() != roi2.getPhiSector() ){continue;}
                R_same = true;
                CanNo.push_back(can1);
            }
        }
        if(R_same){
            std::vector<int> S_ID;S_ID.clear();
            std::vector<int> W_ID;W_ID.clear();
            for(int SID=0;SID!=TGC_coin_n;SID++){
                if((*TGC_coin_type)[SID]!=0){continue;}
                if((*TGC_coin_bunch)[SID]!=2){continue;}
                bool T_same=false;
                for(int j=0;j!=(int)CanNo.size();j++){
                    auto roi = RoIs.at(CanNo.at(j));
                    if((*TGC_coin_isAside)[SID]==roi.getSide()){if((*TGC_coin_isForward)[SID]!=roi.getIsEndcap()){if((*TGC_coin_phi)[SID]==roi.getPhiSector()){T_same = true;}}}
                }
                if(T_same){
                    bool isStrip = (*TGC_coin_isStrip)[SID];
                    if(!isStrip){continue;}
                    bool S_same = false;
                        for(int j=0;j!=(int)CanNo.size();j++){
                            auto roi = RoIs.at(CanNo.at(j));
                            if((*TGC_coin_isAside)[SID]==roi.getSide()){ if((*TGC_coin_isForward)[SID]!=roi.getIsEndcap()){ if((*TGC_coin_phi)[SID]==roi.getPhiSector()){ if((*TGC_coin_trackletId)[SID]==roi.getTrackletIdStrip()){ S_same=true; }}}}
                        }
                    if(S_same){continue;}

                    for(int WID=0;WID!=TGC_coin_n;WID++){
                        if(WID == SID){continue;}
                        if((*TGC_coin_type)[WID]!=0){continue;}
                        if((*TGC_coin_bunch)[WID]!=2){continue;}
                        bool C_same=false;
                        for(int j=0;j!=(int)CanNo.size();j++){
                            auto roi = RoIs.at(CanNo.at(j));
                            if((*TGC_coin_isAside)[WID]==roi.getSide()){if((*TGC_coin_isForward)[WID]!=roi.getIsEndcap()){if((*TGC_coin_phi)[WID]==roi.getPhiSector()){C_same = true;}}}
                        }
                        if(C_same){
                            bool isWire = !(*TGC_coin_isStrip)[WID];
                            if(!isWire){continue;}
                            bool W_same = false;
                            for(int j=0;j!=(int)CanNo.size();j++){
                                auto roi = RoIs.at(CanNo.at(j));
                                if((*TGC_coin_isAside)[WID]==roi.getSide()){ if((*TGC_coin_isForward)[WID]!=roi.getIsEndcap()){ if((*TGC_coin_phi)[WID]==roi.getPhiSector()){ if((*TGC_coin_trackletId)[WID]==roi.getTrackletIdWire()){ W_same=true; }}}}
                            }
                            if(W_same){continue;}

                            if((*TGC_coin_isAside)[SID]==(*TGC_coin_isAside)[WID] && (*TGC_coin_isForward)[SID]==(*TGC_coin_isForward)[WID] && (*TGC_coin_phi)[SID]==(*TGC_coin_phi)[WID]){
                                if(S_ID.size()==0){
                                    S_ID.push_back(SID);
                                    W_ID.push_back(WID);
                                }
                                else{
                                    for(int j=0;j!=(int)S_ID.size();j++){
                                        if((*TGC_coin_isAside)[SID]==(*TGC_coin_isAside)[S_ID.at(j)] && (*TGC_coin_isForward)[SID]==(*TGC_coin_isForward)[S_ID.at(j)] && (*TGC_coin_phi)[SID]==(*TGC_coin_phi)[S_ID.at(j)] && (*TGC_coin_trackletId)[SID]==(*TGC_coin_trackletId)[S_ID.at(j)]){continue;}
                                        if((*TGC_coin_isAside)[WID]==(*TGC_coin_isAside)[W_ID.at(j)] && (*TGC_coin_isForward)[WID]==(*TGC_coin_isForward)[W_ID.at(j)] && (*TGC_coin_phi)[WID]==(*TGC_coin_phi)[W_ID.at(j)] && (*TGC_coin_trackletId)[WID]==(*TGC_coin_trackletId)[W_ID.at(j)]){continue;}
                                            S_ID.push_back(SID);
                                            W_ID.push_back(WID);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if(W_ID.size()!=0){
                std::vector<int> samesec; samesec.clear();
                for(int i=0;i!=(int)W_ID.size();i++){

                    int count = 0;
                    if(samesec.size()!=0){
                        for(int j=0;j!=(int)samesec.size();j++){
                            if((*TGC_coin_phi)[W_ID.at(i)] == samesec.at(j)){count = count+1;}
                        }
                    }
                    if(count>=2){continue;}

                    RoIObj roi_tmp;
                    roi_tmp.setRoI(getRoI(W_ID.at(i),S_ID.at(i)));
                    roi_tmp.setType(2);
                    roi_tmp.setIsEndcap(!(*TGC_coin_isForward)[W_ID.at(i)]);
                    roi_tmp.setSide((*TGC_coin_isAside)[W_ID.at(i)]);
                    roi_tmp.setTrackletIdWire((*TGC_coin_trackletId)[W_ID.at(i)]);
                    roi_tmp.setTrackletIdStrip((*TGC_coin_trackletId)[S_ID.at(i)]);
                    roi_tmp.setX( ((*TGC_coin_x_Out)[W_ID.at(i)] + (*TGC_coin_x_Out)[S_ID.at(i)])/2 );
                    roi_tmp.setY( ((*TGC_coin_y_Out)[W_ID.at(i)] + (*TGC_coin_y_Out)[S_ID.at(i)])/2 );
                    roi_tmp.setZ( ((*TGC_coin_z_Out)[W_ID.at(i)] + (*TGC_coin_z_Out)[S_ID.at(i)])/2 );
                    roi_tmp.setR((*TGC_coin_width_Out)[W_ID.at(i)]);
                    roi_tmp.setPhi((*TGC_coin_width_Out)[S_ID.at(i)]);
                    roi_tmp.setPhiSector((*TGC_coin_phi)[W_ID.at(i)]);
                    roi_tmp.setTypeDR((*TGC_coin_type)[W_ID.at(i)]);
                    roi_tmp.setDR((*TGC_coin_delta)[W_ID.at(i)]);
                    roi_tmp.setTypeDPhi((*TGC_coin_type)[S_ID.at(i)]);
                    roi_tmp.setDPhi((*TGC_coin_delta)[S_ID.at(i)]);
                    RoIs.push_back(roi_tmp);

                    samesec.push_back((*TGC_coin_phi)[W_ID.at(i)]);
                }
            } 
        }
    }
}
