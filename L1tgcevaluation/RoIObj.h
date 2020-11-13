#ifndef ROIOBJ_h
#define ROIOBJ_h


#include <iostream>


class RoIObj{

 public:
  RoIObj();
  //~RoIObj();

  RoIObj& operator+(const RoIObj& right);
  RoIObj& operator+=(const RoIObj& right);



  void setType(int inType){m_type=inType;};
  void setDR(int inDR);
  void setDPhi(int inDPhi);
  void setTypeDR(int inTypeDR){m_typeDR=inTypeDR;};
  void setTypeDPhi(int inTypeDPhi){m_typeDPhi=inTypeDPhi;};
  void setSide(int inSide){m_side=inSide;};
  void setRoI(int inRoI){m_roi=inRoI;};
  void setPhiSector(int inPhiSector){m_phiSector=inPhiSector;};
  void setIsEndcap(bool inIsEndcap){m_isEndcap=inIsEndcap;};
  void setTrackletIdWire(int inTrackletIdWire){m_trackletIdWire=inTrackletIdWire;};
  void setTrackletIdStrip(int inTrackletIdStrip){m_trackletIdStrip=inTrackletIdStrip;};
  void setPt(int InPt){m_pT=InPt;};
  void setX(float Inx){m_x=Inx;};
  void setY(float Iny){m_y=Iny;};
  void setZ(float Inz){m_z=Inz;};
  void setR(float InR){m_R=InR;};
  void setPhi(float InPhi){m_Phi=InPhi;}
  void setCharge(float InCharge){m_Charge=InCharge;}

  int getType()const{return m_type;};
  int getDR()const{return m_DR;};
  int getDPhi()const{return m_DPhi;};
  int getTypeDPhi()const{return m_typeDPhi;};
  int getTypeDR()const{return m_typeDR;};
  int getSide()const{return m_side;};
  int getRoI()const{return m_roi;};
  int getPhiSector()const{return m_phiSector;};
  bool getIsEndcap()const{return m_isEndcap;};
  int getTrackletIdWire()const{return m_trackletIdWire;};
  int getTrackletIdStrip()const{return m_trackletIdStrip;};
  int getPt(){return m_pT;};
  float getX()const{return m_x;};
  float getY()const{return m_y;};
  float getZ()const{return m_z;};
  float getR()const{return m_R;};
  float getPhi()const{return m_Phi;}
  int getCharge()const{return m_Charge;}

  int getOctant()const;
  int getModuleId()const;
  int getCoinType()const;
  int getModifedPhiSector()const;

  bool isSameRoI(RoIObj roi)const;

 private:

  int m_side;//0=A 1=C
  int m_roi;
  int m_DR;
  int m_DPhi;
  int m_typeDR;
  int m_typeDPhi;
  int m_phiSector;
  bool m_isEndcap;

  bool m_isSetDR;
  bool m_isSetDPhi;

  int m_trackletIdWire;
  int m_trackletIdStrip;

  int m_pT;
  int m_type;
  int m_Charge;

  float m_x;
  float m_y;
  float m_z;
  float m_R;
  float m_Phi;

  mutable int m_octant;
  mutable int m_moduleId;
  mutable int m_cointype;


};

#endif
