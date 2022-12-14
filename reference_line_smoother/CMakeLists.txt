project(Reference_line_smoother)

cmake_minimum_required(VERSION 2.80)
set(CMAKE_CXX_STANDARD 17)

find_package(OsqpEigen REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${OsqpEigen_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIR})

include_directories(${CMAKE_SOURCE_DIR}/include)
add_library(optimization_base_lib STATIC ${CMAKE_SOURCE_DIR}/src/optimization_base.cpp
                                    ${CMAKE_SOURCE_DIR}/src/fem_reference_line_smoother.cpp)
target_link_libraries(optimization_base_lib OsqpEigen::OsqpEigen ${Boost_LIBRARIES})

add_executable(${PROJECT_NAME}_main ${CMAKE_SOURCE_DIR}/src/main.cpp)
target_link_libraries(${PROJECT_NAME}_main optimization_base_lib)