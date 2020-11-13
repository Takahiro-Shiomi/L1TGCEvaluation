/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "TGCRPhiCoincidenceMap.h"
#include "TGCNumbering.h"
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <iostream>


int TGCRPhiCoincidenceMap::test_Run3(int octantId, int moduleId, int subsector, 
                                     int type, int dr, int dphi) const
{
//    std::cout<<"test_Run3__START"<<std::endl;  
  // check pt range
  if (type<0 || type>=N_TMap ) return false; 

  int sector=(moduleId-2)/3+octantId*3;


  int phimod2 = (moduleId==2||moduleId==5||moduleId==8)&&(sector%2==1) ? 1 : 0;
  int addr=SUBSECTORADD(subsector, moduleId, phimod2,type);

  std::map<int, std::map<int, std::map<int, char> > >::const_iterator it=mapDB_Run3.find(addr);
  if (it==mapDB_Run3.end()) return false;

  std::map<int, std::map<int, char> > mapR = it->second;
  std::map<int, std::map<int, char> >::const_iterator itR=mapR.find(dr);
  if (itR==mapR.end()) return false;

  std::map<int, char> mapPhi = itR->second;
  std::map<int, char>::const_iterator itPhi=mapPhi.find(dphi);

  if (itPhi==mapPhi.end()) return false;

  char pT_char=itPhi->second;
  int pT_int=pTdef.find(pT_char)->second;

  //std::cout<<pT_int<<std::endl;
  return  pT_int;
}




TGCRPhiCoincidenceMap::TGCRPhiCoincidenceMap()
 : numberOfDR(0), numberOfDPhi(0),
   m_side(),
   m_octant(0),
   m_fullCW(false)
{
  // set message label

  this->readMap_Run3();  // read Coincidence Map 
}

TGCRPhiCoincidenceMap::TGCRPhiCoincidenceMap(const std::string& cw_path,
                                             int sideId, int octantId)
 : numberOfDR(0), numberOfDPhi(0),
   m_cw_path(cw_path),
   m_side(sideId),
   m_octant(octantId),
   m_fullCW(false)
{
  // set message label

  this->readMap_Run3();  // read Coincidence Map 
}


TGCRPhiCoincidenceMap::~TGCRPhiCoincidenceMap()
{}

bool TGCRPhiCoincidenceMap::readMap_Run3() 
{
  const int NumberOfModuleType=12;
  const int ModuleNumber[NumberOfModuleType]  =
      {  0,  1,   2,   2,  3,  4,   5,   5,  6,  7,   8,  8 };
  const std::string ModuleName[NumberOfModuleType]=
      {"0a","1a","2a","2b","3a","4a","5a","5b","6a","7a","8a","8b"};
  const std::string SideName[kTotalNumTGCSide] = {"a","c"};
  const std::string OctantName[NumberOfOctant] =
      {  "0", "1", "2", "3", "4", "5", "6", "7"};

  // initialize
  std::string buf;
  std::string fn, fullName, tag;
  int ssId;
  char delimiter = '\n';

  // loop over all files...
  for(int iModule=0; iModule<NumberOfModuleType; iModule+=1) {
    int phimod2=ModuleName[iModule].find("b")!=std::string::npos ? 1 : 0;
    std::ostringstream modName;
    std::string fn =  m_cw_path+"cm_"+ SideName[m_side]+OctantName[m_octant]+ModuleName[iModule]+"_v0001.db";
    bool Forward_type1=(ModuleName[iModule]=="2b"||ModuleName[iModule]=="5a"||ModuleName[iModule]=="8b");
    bool Forward_type2=(ModuleName[iModule]=="2a"||ModuleName[iModule]=="5b"||ModuleName[iModule]=="8a");
    if(m_octant%2==0 && Forward_type1){continue;}
    if(m_octant%2==1 && Forward_type2){continue;}
    int type = -1;
    int lDR, hDR, lDPhi, hDPhi;

    std::ifstream file(fn.c_str(),std::ios::in);    
    while(getline(file,buf,delimiter)){
      std::istringstream header(buf);
      header>>tag;
      if(tag=="#") { // read header part.     
        header>>ssId>>lDR>>hDR>>lDPhi>>hDPhi;
        type = getTYPE( lDR, hDR, lDPhi, hDPhi );
        if(hDR==15){DR_offset=-15;}
        else if(hDR==7){DR_offset=-7;}
        if(hDPhi==7){DPhi_offset=-7;}
        else if(hDPhi==3){DPhi_offset=-3;}
        // check moduleNumber and ptLevel
        if( type<0 ) {
          break;
        }

        // get window data
        std::map<int, std::map<int, char> >  bWindow;//<R,<~>>
        char pT;
        for(int ir=0; ir<=hDR-DR_offset; ir++) {
          getline(file,buf,delimiter);
          std::map<int, char> aWindow;//<Phi,pT>
          for(int iphi=0; iphi<=hDPhi-DPhi_offset; iphi++){
            pT = buf[iphi];
            if (pT=='X') continue; // none of window is opened in this dR
            aWindow[iphi+DPhi_offset] = pT;
          }
          // Warning : no window 
          if (aWindow.size()==0) {
          } else {
            bWindow[ir+DR_offset]=aWindow;
          }
        }
        int addr = SUBSECTORADD(ssId,ModuleNumber[iModule],phimod2,type);
        if (mapDB_Run3.find(addr)!=mapDB_Run3.end()) {
        } else {
          mapDB_Run3[addr]=bWindow;
        }
      }
    }
      
  }

  return true;
}


bool TGCRPhiCoincidenceMap::readHotRoIList()
{
  // initialize
  std::string buf;
  std::string fullName, tag;    char delimiter = '\n';
  std::ostringstream modName;
  std::string fn = "HotRoIList.db";


  int RoI,module,module_previous=-1;
  std::ifstream file(fn.c_str(),std::ios::in); 
    
  std::map<int,bool> mapRoI;
  while(getline(file,buf,delimiter)) {
    std::istringstream header(buf);
    header>>RoI>>module;
    if(module!=module_previous){
      if(module_previous!=-1) MapIsHotRoI[module_previous] = mapRoI;
      module_previous=module;
      mapRoI.clear();
    }
    mapRoI[RoI] = true;
  }
  MapIsHotRoI[module_previous] = mapRoI;

  return true;
}


bool TGCRPhiCoincidenceMap::test_HotRoI(int moduleId,int subsector) const
{
  std::map<int, std::map<int, bool> >::const_iterator itModule=MapIsHotRoI.find(moduleId);
  if (itModule==MapIsHotRoI.end()) return true;

  std::map<int, bool> mapRoI = itModule->second;
  std::map<int, bool>::const_iterator itRoI=mapRoI.find(subsector);
  if (itRoI==mapRoI.end()) return true;

  return false;
}


