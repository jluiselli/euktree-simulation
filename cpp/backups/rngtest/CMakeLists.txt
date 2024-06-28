cmake_minimum_required(VERSION 3.21...3.26)

project(randtest CXX)

set(CMAKE_CXX_FLAGS "-O3 -Wextra ")

## include necessary modules
#
include(FetchContent)
include(CheckLanguage)

## find p2rng, if not installed fetch it...
#
find_package(p2rng CONFIG
  COMPONENTS openmp
  OPTIONAL_COMPONENTS cuda oneapi rocm
)

if(NOT p2rng_FOUND)
  message(STATUS "Fetching p2rng library...")
  FetchContent_Declare(
    p2rng
    GIT_REPOSITORY https://github.com/arminms/p2rng.git
    GIT_TAG main
  )
  # setting required p2rng components
  set(P2RNG_COMPONENTS openmp cuda oneapi rocm
    CACHE STRING "Required components"
  )
  FetchContent_MakeAvailable(p2rng)
endif()


add_executable(randtest test_binom.cpp)
target_link_libraries(randtest PRIVATE p2rng::openmp)
