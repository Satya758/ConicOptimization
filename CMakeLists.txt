
set(CMAKE_BUILD_TYPE "Debug")
SET(CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libstdc++")

if(UNIX)
    add_definitions( "-Wall" )
endif()

cmake_minimum_required(VERSION 2.8)

####### Required Libraries ########################
set(EIGEN_INC "" CACHE PATH "Eigen Include Directory") #No default value
set(BOOST_INC "" CACHE PATH "Boost Lib Include Directory") #No default value
####### Required Libraries ########################

project(ConicOptimization)



add_subdirectory(src)
include_directories(.)

include_directories(SYSTEM ${EIGEN_INC})
include_directories(SYSTEM ${BOOST_INC})

add_executable(${PROJECT_NAME} ./src/main.cpp)

install(TARGETS ConicOptimization RUNTIME DESTINATION bin)


