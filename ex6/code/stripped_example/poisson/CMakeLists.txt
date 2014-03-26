project(poisson-example)
cmake_minimum_required(VERSION 2.6)

add_subdirectory(../common common)

enable_language(Fortran)

add_definitions(${COMMON_DEFINITIONS})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMMON_C_FLAGS}")
include_directories(${COMMON_INCLUDE_PATHS})

add_library(poisson_common solvers.c poissoncommon.c fst.f)

if(MPI_FOUND)
  add_executable(poisson2D-MPI poisson2D-MPI.c)
  target_link_libraries(poisson2D-MPI poisson_common ${COMMON_LIBRARIES})
endif()

add_executable(poisson2D poisson2D.c)
target_link_libraries(poisson2D poisson_common ${COMMON_LIBRARIES})
