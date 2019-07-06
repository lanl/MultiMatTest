cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# --- custom targets: ---
INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)

project (MultiMatTest)
# The version number.
set (MultiMatTest_VERSION_MAJOR 2)
set (MultiMatTest_VERSION_MINOR 0.7)
set (MultiMatTest_VERSION "${MultiMatTest_VERSION_MAJOR}.${MultiMatTest_VERSION_MINOR}")

set (CMAKE_C_STANDARD 99)
 
if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_C_FLAGS "-std=c99")
#set(CMAKE_C_FLAGS "-std=c99 -fstrict-aliasing")
#set(CMAKE_CXX_FLAGS "-fstrict-aliasing")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -Wall -fno-tree-vectorize")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -Wall -fno-tree-vectorize")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -O3 -march=native -mtune=native")
#set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O3 -march=native -mtune=native")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -O2 -fno-tree-vectorize -march=native -mtune=native")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O2 -fno-tree-vectorize -march=native -mtune=native")

#find_package(OpenCL REQUIRED)
##  OPENCL_FOUND        : TRUE if found, FALSE otherwise
##  OPENCL_INCLUDE_DIRS : Include directories for OpenCL
##  OPENCL_LIBRARIES    : The libraries to link against
#if (OPENCL_FOUND)
#   include_directories(${OPENCL_INCLUDE_DIRS})
#   if (HAVE_CL_DOUBLE)
#      add_definitions(-DHAVE_CL_DOUBLE)
#   endif (HAVE_CL_DOUBLE)
#endif (OPENCL_FOUND)

## configure a header file to pass some of the CMake settings
## to the source code
#configure_file ( "${PROJECT_SOURCE_DIR}/config.h.in" "${PROJECT_BINARY_DIR}/config.h")
# 
## add the binary tree to the search path for include files
## so that we will find config.h
#include_directories("${PROJECT_BINARY_DIR}")

########### MultiMatTest target ###############

set(MultiMatTest_CXX_SRCS MultiMatTest.cpp input.cpp)
#set(MultiMatTest_C_SRCS ezcl.c genmalloc.c timer.c)
set(MultiMatTest_C_SRCS genmalloc.c timer.c)
#set(MultiMatTest_H_SRCS input.h ezcl.h genmalloc.h timer.h)
set(MultiMatTest_H_SRCS input.h genmalloc.h timer.h)
set(MultiMatTest_SRCS ${MultiMatTest_CXX_SRCS} ${MultiMatTest_C_SRCS} ${MultiMatTest_H_SRCS})

add_executable(MultiMatTest ${MultiMatTest_SRCS})

set_target_properties(MultiMatTest PROPERTIES COMPILE_DEFINITIONS HAVE_OPENCL)

target_link_libraries(MultiMatTest ${OPENCL_LIBRARIES})

install(TARGETS MultiMatTest DESTINATION bin)

########### install files ###############

########### make clean ##################
add_custom_target(clean-cmake-files
   COMMAND ${CMAKE_COMMAND} -P clean-all.cmake
)

## clean-all.cmake
set(cmake_generated ${CMAKE_BINARY_DIR}/CMakeCache.txt
                    ${CMAKE_BINARY_DIR}/cmake_install.cmake  
                    ${CMAKE_BINARY_DIR}/Makefile
                    ${CMAKE_BINARY_DIR}/CMakeFiles
)

foreach(file ${cmake_generated})

  if (EXISTS ${file})
     file(REMOVE_RECURSE ${file})
  endif()

endforeach(file)

########### make dist ##############

set (EXTRA_FILES AUTHORS CMakeLists.txt config.h.in multimat_kern.cl)
add_custom_target(dist
   rm -rf MultiMatTest-v${MultiMatTest_VERSION}
   COMMAND mkdir MultiMatTest-v${MultiMatTest_VERSION}
   COMMAND cp ${MultiMatTest_SRCS} ${EXTRA_FILES} MultiMatTest-v${MultiMatTest_VERSION}
   COMMAND cp -R cmake MultiMatTest-v${MultiMatTest_VERSION}
   COMMAND tar -czvf MultiMatTest-v${MultiMatTest_VERSION}.tgz MultiMatTest-v${MultiMatTest_VERSION}
   COMMAND rm -rf MultiMatTest-v${MultiMatTest_VERSION}
   )

