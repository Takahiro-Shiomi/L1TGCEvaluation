add_executable(L1TGCEva
	       main.cxx
	       TGCRPhiCoincidenceMap.cxx
	       TGCCoinData.cxx
	       TrigT1TGCNtuple.C
	       physics.C
	       HistogramManager.C
	       RoIObj.C
	      ) 
# Include directories for this project 
set(INCLUDE_DIR
    ${PROJECT_SOURCE_DIR}/L1tgcevaluation
    )
 # Add a include files 
include_directories("${INCLUDE_DIR}") 
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

