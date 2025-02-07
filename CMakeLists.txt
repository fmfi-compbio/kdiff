cmake_minimum_required(VERSION 3.14)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(ExternalProject)
include(FetchContent)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(kdiff VERSION 1.0.0)
add_executable(kdiff main.cpp)
set_target_properties(kdiff PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
target_compile_options(kdiff PRIVATE -Wall -Wextra -Wno-unused-variable -Wno-sign-compare)

# find_package(OpenMP REQUIRED)

# spdlog
message(STATUS "spdlog will be built from source")
FetchContent_Declare(spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG 7e635fca68d014934b4af8a1cf874f63989352b7) # v1.12.0
FetchContent_MakeAvailable(spdlog)

# KMC
message(STATUS "KMC3 will be built from source")
ExternalProject_Add(kmc
    GIT_REPOSITORY https://github.com/refresh-bio/KMC.git
    GIT_TAG 2049648ecad2c2d99b827c5d4fecff158122d139
    UPDATE_COMMAND ""
    BUILD_IN_SOURCE 1
    PATCH_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make kmc
    INSTALL_COMMAND ""
)
ExternalProject_Get_Property(kmc SOURCE_DIR)
SET(KMC_SOURCE_DIR ${SOURCE_DIR})
SET(KMC_INCLUDE_DIR ${KMC_SOURCE_DIR})
add_library(KMC STATIC IMPORTED)
set_target_properties(KMC PROPERTIES IMPORTED_LOCATION ${KMC_INCLUDE_DIR}/bin/libkmc_core.a)
add_dependencies(KMC kmc)


add_dependencies(kdiff spdlog KMC)

target_include_directories(kdiff
  PRIVATE ${KMC_INCLUDE_DIR}
  PRIVATE ${PROJECT_SOURCE_DIR}
)

target_link_libraries(kdiff
  PUBLIC KMC
  PUBLIC z
  PUBLIC spdlog::spdlog
  # PUBLIC OpenMP::OpenMP_CXX
)
