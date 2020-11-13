# L1TGCEvaluation


asetup 21.3.5,Athena


L1tgcevaluation  build  run

//In build
cmake ../L1tgcevaluation
cd ../run
cp ../build/src/L1TGCEva .


./L1TGCEva -BWCW ../../TurnOnCurve/athena/Trigger/TrigT1/TrigT1TGC/share/BW_CW_Run3Format/ -Input group.det-muon/group.det-muon.19069326.L1TGCNtuple._000105.root


-BWCW : path of CW for BW.
-Input : file path and name.(L1TGCNtuple)


