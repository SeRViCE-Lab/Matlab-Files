# cmake needs this line
cmake_minimum_required(VERSION 2.8.12)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBT_USE_DOUBLE_PRECISION -Wall")
# Unused warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wuninitialized -Winit-self -Wunused-function -Wunused-label -Wunused-variable -Wunused-but-set-variable -Wunused-but-set-parameter")
# Additional warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Warray-bounds -Wtype-limits -Wreturn-type -Wsequence-point -Wparentheses -Wmissing-braces -Wchar-subscripts -Wswitch -Wno-write-strings -Wenum-compare -Wempty-body -Wlogical-op")
SET(CMAKE_BUILD_TYPE RelWithDebInfo)

#Project Name
PROJECT(engdemo)

find_package(PkgConfig)
find_package(Matlab REQUIRED)

INCLUDE_DIRECTORIES(${Matlab_INCLUDE_DIRS})

SET(LIBRARIES ${Matlab_LIBRARIES}
        ${Matlab_ENG_LIBRARY})

OPTION(ENABLE_CXX11 "Enable C++11 Suport" ON)

IF(ENABLE_CXX11)
  INCLUDE(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
  CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
  IF(COMPILER_SUPPORTS_CXX11)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  ELSEIF(COMPILER_SUPPORTS_CXX0X)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
  ELSE()
    MESSAGE(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
  ENDIF()
ENDIF(ENABLE_CXX11)

#Show found directories
message(STATUS "\nMatlab Library Status")
message(STATUS "     Version: ${matlab_get_version_from_release_name()}"
message(STATUS "     Directory: ${Matlab_ROOT_DIR}")
message(STATUS "     Include DIRS: ${Matlab_INCLUDE_DIRS}")
message(STATUS "     Mex Library: ${Matlab_MEX_LIBRARY}")
message(STATUS "     Matlab Eng Lib: ${Matlab_ENG_LIBRARY}")
)

LINK_DIRECTORIES(${Matlab_Libraries_DIRS})

add_executable(engdemo engdemo.cpp)

target_link_libraries(engdemo ${Matlab_LIBRARIES})
