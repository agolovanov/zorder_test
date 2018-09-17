cmake_minimum_required (VERSION 3.10.2)
project (zorder_test)
set (CMAKE_CXX_FLAGS "-Wall -Wextra -Ofast")
set (CMAKE_CXX_STANDARD 11)
find_package( Boost )
include_directories( ${Boost_INCLUDE_DIRS} )
add_executable(zorder_test main.cpp PerformanceTest.cpp)
target_link_libraries(zorder_test Boost::boost)

