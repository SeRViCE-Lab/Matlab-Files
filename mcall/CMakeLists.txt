# example project shows how to call matlab in cpp

cmake_minimum_required (VERSION 2.8)

project (engdemo)


# matlab
#set(MATLAB_INCLUDE_DIR "/usr/local/MATLAB/R2015a/extern/include")
set(MATLAB_INCLUDE_DIR "/usr/local/MATLAB/R2015a/extern/include")
set(MATLAB_LIBRARY_DIR "/usr/local/MATLAB/R2015a/bin/glnxa64")
set(MATLAB_LIBRARIES
  ${MATLAB_LIBRARY_DIR}/libeng.so
  ${MATLAB_LIBRARY_DIR}/libmx.so
  ${MATLAB_LIBRARY_DIR}/libmex.so
  )

include_directories (${MATLAB_INCLUDE_DIR})

add_executable(engdemo engdemo.cpp )
target_link_libraries(engdemo ${MATLAB_LIBRARIES})