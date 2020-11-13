################################################################################ 
# Package:  MAKE_BWCW 
################################################################################ 
# Set the minimum required CMake version: 
cmake_minimum_required( VERSION 3.2 FATAL_ERROR )

project(MAKE_BWCW CXX)
set(CMAKE_CXX_FLAGS "-pthread -std=c++1y -m64 -I/cvmfs/atlas.cern.ch/repo/sw/software/21.3/sw/lcg/releases/LCG_88/ROOT/6.08.06/x86_64-slc6-gcc62-opt/include -L/cvmfs/atlas.cern.ch/repo/sw/software/21.3/sw/lcg/releases/LCG_88/ROOT/6.08.06/x86_64-slc6-gcc62-opt/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic")

add_subdirectory(src)