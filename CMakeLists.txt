PROJECT(Matlab_LRSpline)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

# Required packages
FIND_PACKAGE(Matlab )
FIND_PACKAGE(LRSpline )

file(COPY ${PROJECT_SOURCE_DIR}/src/LRSplineSurface.m DESTINATION ${CMAKE_BINARY_DIR}/lib)
file(COPY ${PROJECT_SOURCE_DIR}/src/getBSplineBasisAndDerivative.m DESTINATION ${CMAKE_BINARY_DIR}/lib)
ADD_CUSTOM_TARGET(matlab-lrspline ALL COMMAND ${MATLAB_MEX} CC=${MATLAB_CC} -I${PROJECT_SOURCE_DIR}/include -I${CMAKE_BINARY_DIR}/include ${PROJECT_SOURCE_DIR}/src/lrsplinesurface_interface.cpp -lLRSpline -outdir ${CMAKE_BINARY_DIR}/lib)

