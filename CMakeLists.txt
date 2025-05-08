cmake_minimum_required(VERSION 3.18)
project(GEO1004A1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_PREFIX_PATH /opt/homebrew)

find_package(CGAL REQUIRED COMPONENTS Core)

find_package(Boost REQUIRED)

add_executable(GEO1004A1 hw1.cpp)

target_link_libraries(GEO1004A1 CGAL::CGAL CGAL::CGAL_Core Boost::boost)
