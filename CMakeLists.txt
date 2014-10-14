PROJECT(Matlab_LRSpline)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

# Required packages
FIND_PACKAGE(Matlab REQUIRED)
FIND_PACKAGE(LRSpline REQUIRED)

# Use same definitions as was compiled with LRspline 
LIST(APPEND DEFINITIONS ${LRSpline_DEFINITIONS})

# Set folder where to build everything
SET(BUILD_DIR "${PROJECT_SOURCE_DIR}/lib")
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${BUILD_DIR})

# mex-compile c++ and put .m-files in the right place
ADD_CUSTOM_TARGET(matlab-lrspline ALL
  COMMAND ${MATLAB_MEX} CC=${MATLAB_CC} -I${PROJECT_SOURCE_DIR}/include ${DEFINITIONS} ${PROJECT_SOURCE_DIR}/src/lrsplinesurface_interface.cpp -lLRSpline -outdir ${BUILD_DIR}
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/src/matlab ${BUILD_DIR})

