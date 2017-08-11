PROJECT(Matlab_LRSpline)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

# Required packages
FIND_PACKAGE(Matlab REQUIRED)
FIND_PACKAGE(LRSpline REQUIRED)
FIND_PACKAGE(GoTools)

# Use same definitions as was compiled with LRspline 
LIST(APPEND DEFINITIONS ${LRSpline_DEFINITIONS})
# Add in all include directories and libraries used by LRSplines
LIST(APPEND INCLUDES    -I${LRSpline_INCLUDE_DIRS})
LIST(APPEND LIBS        ${LRSpline_LIBRARIES})

IF(GoTools_FOUND)
  # Add in all include directories and libraries used by GoTools
  LIST(APPEND INCLUDES -I${GoTools_INCLUDE_DIRS})
  LIST(APPEND LIBS     ${GoTools_LIBRARIES})
ENDIF(GoTools_FOUND)

# Set folder where to build everything
SET(BUILD_DIR "${PROJECT_SOURCE_DIR}/lib")
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${BUILD_DIR})

# mex-compile c++ and put .m-files in the right place
ADD_CUSTOM_TARGET(matlab-lrspline ALL
  COMMAND ${MATLAB_MEX} CC=${MATLAB_CC} -I${PROJECT_SOURCE_DIR}/include ${INCLUDES} ${DEFINITIONS} ${PROJECT_SOURCE_DIR}/src/lrsplinesurface_interface.cpp ${LIBS} -outdir ${BUILD_DIR} -O
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/src/matlab ${BUILD_DIR})

