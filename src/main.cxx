#include <iostream>
#include <math.h>
#include <stdlib.h>
#include "TGCRPhiCoincidenceMap.h"
#include "physics.h"


int main(int arg, char *argv[]){

  std::string BWCW_path="";
  std::string file="";

  for(int para=0;para!=arg-1;para++){
    std::string config=argv[para+1];
    if(config=="-BWCW"){BWCW_path=argv[para+2];}
    if(config=="-h"){std::cout<<"-BWCW"<<std::endl; return 0;}
    if(config=="-Input"){file = argv[para+2];}
  }


  if(BWCW_path==""){std::cout<<"Please set BWCW path (-BWCW)"<<std::endl;}
  else{std::cout<<"BWCW path :: "<<BWCW_path.c_str()<<std::endl;}
  if(file==""){std::cout<<"Please set Input path (-Input)"<<std::endl;}
  else{  std::cout<<"Input File :: "<<file.c_str()<<std::endl;}


  TChain *tc = new TChain("physics");
  tc->Add(file.c_str());
  TTree *tt = (TTree*)tc;

  physics* p= new physics(tt);


  p->Loop(BWCW_path);

  return 0;
}
