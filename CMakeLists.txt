
set(CMAKE_BUILD_TYPE "Debug")
SET(CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libstdc++")

if(UNIX)
    add_definitions( "-Wall" )
endif()

cmake_minimum_required(VERSION 2.8)

project(ConicOptimization)

####### Required Libraries ########################
set(EIGEN_INC "" CACHE PATH "Eigen Include Directory") #No default value
set(BOOST_INC "" CACHE PATH "Boost Lib Include Directory") #No default value
####### Required Libraries ########################

add_subdirectory(src)
include_directories(./src)

include_directories(SYSTEM ${EIGEN_INC})
include_directories(SYSTEM ${BOOST_INC})

#Pastix
include_directories(SYSTEM "/home/satya/ThirdParty/Pastix-5.2.2/include")

find_library(BLAS_LIBRARY blas PATHS ${LIB_INSTALL_DIR})
find_library(HWLOC_LIBRARY hwloc PATHS ${LIB_INSTALL_DIR})
find_library(PTSCOTCH ptscotch PATHS "/home/satya/ThirdParty/Scotch-6.0.0/lib")
find_library(SCOTCH scotch PATHS "/home/satya/ThirdParty/Scotch-6.0.0/lib")
find_library(SCOTCH_ERR scotcherrexit PATHS "/home/satya/ThirdParty/Scotch-6.0.0/lib")
find_library(PTSCOTCH_ERR ptscotcherrexit PATHS "/home/satya/ThirdParty/Scotch-6.0.0/lib")
find_library(PASTIX pastix PATHS "/home/satya/ThirdParty/Pastix-5.2.2/lib")
#Pastix

set(MAIN_SRC ./src/main.cpp)

add_executable(ConicOptimization ${MAIN_SRC})

target_link_libraries(${PROJECT_NAME}
                      ${PASTIX} ${PTSCOTCH} ${SCOTCH} ${SCOTCH_ERR} m rt  pthread
                      ${HWLOC_LIBRARY}
                      ${BLAS_LIBRARY} gfortran    )

install(TARGETS ConicOptimization RUNTIME DESTINATION bin)
