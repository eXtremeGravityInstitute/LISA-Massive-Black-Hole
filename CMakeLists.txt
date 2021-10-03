# mbh
cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

project(mbh VERSION 0.0.0 LANGUAGES C)

# Paths
include_directories(/opt/local/include)
link_directories(/opt/local/lib)

# External packages
find_package(GSL REQUIRED)
find_package(OpenMP REQUIRED)

# Set compiler
set(CMAKE_C_COMPILER mpicc)

# set cflags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffast-math -Wall -Wextra -g -ftree-vectorize -std=gnu99")

# Internal products

add_library(mbh STATIC Utils.c IMRPhenomD.c IMRPhenomD_internals.c Response.c IMRPhenomD.h IMRPhenomD_internals.h Declarations.h Constants.h)
include_directories(SYSTEM ${GSL_INCLUDE_DIRS})
target_link_libraries(mbh ${GSL_LIBRARIES} m)
set_target_properties(mbh PROPERTIES PUBLIC_HEADER "mbh.h;IMRPhenomD.h;IMRPhenomD_internals.h;Declarations.h;Constants.h")
install(TARGETS mbh EXPORT MBHConfig LIBRARY DESTINATION lib PUBLIC_HEADER DESTINATION include)
export(TARGETS mbh NAMESPACE MBH:: FILE "MBHConfig.cmake")
install(EXPORT MBHConfig DESTINATION lib/cmake/mbh NAMESPACE MBH::)
