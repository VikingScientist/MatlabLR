PROJECT(Matlab_LRSpline)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

# Required packages
FIND_PACKAGE(Matlab )
FIND_PACKAGE(LRSpline )

MESSAGE("matlab mex is ${MATLAB_MEX} and ${MATLAB_CC} for crying ouut loud")

file(COPY ${PROJECT_SOURCE_DIR}/src/LRSplineSurface.m DESTINATION ${CMAKE_BINARY_DIR}/lib)
ADD_CUSTOM_TARGET(matlab-lrspline ALL COMMAND ${MATLAB_MEX} CC=${MATLAB_CC} -I${PROJECT_SOURCE_DIR}/include -I${CMAKE_BINARY_DIR}/include ${PROJECT_SOURCE_DIR}/src/lrsplinesurface_interface.cpp -lLRSpline -outdir ${CMAKE_BINARY_DIR}/lib)

