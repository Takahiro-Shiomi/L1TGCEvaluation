#include "RoIObj.h"

#include <iostream>
#include <iomanip> 
#include <fstream>

RoIObj::RoIObj():
  m_pT(-99),
  m_type(-99),
  m_side(-99),
  m_roi(-1),
  m_DR(999),
  m_DPhi(999),
  m_typeDR(-1),
  m_typeDPhi(-1),
  m_phiSector(-1),
  m_isEndcap(false),
  m_isSetDR(false),
  m_isSetDPhi(false),
  m_octant(-99),
  m_moduleId(-99),
  m_cointype(-99){

}

RoIObj& RoIObj::operator+=(const RoIObj& right){
  if(this != &right){
    if(!m_isSetDR){m_DR=right.getDR();}
    else{std::cout<<__FILE__<<"  WARNING :: DR"<<std::endl;}
    if(!m_isSetDPhi){m_DPhi=right.getDPhi();}
    else{std::cout<<__FILE__<<"  WARNING :: DPhi"<<std::endl;}
  }
}


void RoIObj::setDR(int inDR){m_DR=inDR; m_isSetDR=true;}
void RoIObj::setDPhi(int inDPhi){m_DPhi=inDPhi; m_isSetDPhi=true;}

bool RoIObj::isSameRoI(RoIObj roi)const{

  if(this->getRoI()==roi.getRoI() && 
     this->getSide()==roi.getSide() &&
     this->getPhiSector()==roi.getPhiSector() &&
     this->getIsEndcap()==roi.getIsEndcap()
     ){
    return true;
  }

  return false;

};



int RoIObj::getOctant()const{

    if(m_isEndcap==true){
        if(m_phiSector>=4 && m_phiSector<=45){
        m_octant = (m_phiSector+2)/6;
        }
        if((m_phiSector>=1 && m_phiSector<=3)||(m_phiSector>=46 && m_phiSector<=48)){
        m_octant = 0;
        }
    }
    if(m_isEndcap==false){
        if(m_phiSector>=3 && m_phiSector<=23){
        m_octant = (m_phiSector)/3;
        }
        if((m_phiSector>=1 && m_phiSector<=2)||(m_phiSector==24)){
        m_octant = 0;
        }
    }

  return m_octant;

};
int RoIObj::getModuleId()const{

    if(m_isEndcap==true){

        if((m_phiSector+2) % 6==1){m_moduleId=0;}
        if((m_phiSector+2) % 6==2){m_moduleId=1;}
        if((m_phiSector+2) % 6==3){m_moduleId=3;}
        if((m_phiSector+2) % 6==4){m_moduleId=4;}
        if((m_phiSector+2) % 6==5){m_moduleId=6;}
        if((m_phiSector+2) % 6==0){m_moduleId=7;}
        
    }
    if(m_isEndcap==false){
        
        if(m_phiSector % 3==1){m_moduleId=5;}
        if(m_phiSector % 3==2){m_moduleId=8;}
        if(m_phiSector % 3==0){m_moduleId=2;}
        
    }

  return m_moduleId;

};
int RoIObj::getCoinType()const{

    if(m_typeDR==1 && m_typeDPhi==1){m_cointype=0;}
    if(m_typeDR==1 && m_typeDPhi==0){m_cointype=1;}
    if(m_typeDR==0 && m_typeDPhi==1){m_cointype=2;}
    if(m_typeDR==0 && m_typeDPhi==0){m_cointype=3;}
  
  return m_cointype;

};

int RoIObj::getModifedPhiSector()const{

    if(this->getIsEndcap()){
        if(this->getPhiSector()==48 || this->getPhiSector()==47){
            return this->getPhiSector()-47;
        }
        return this->getPhiSector()+1;
    }
    else{

        if(this->getPhiSector()==24){
            return 0;
        }
        return this->getPhiSector();

    }

};
