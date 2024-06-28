cmake_minimum_required(VERSION 3.21...3.26)

project(simchr CXX)

#requires openmp and p2rng (https://github.com/arminms/p2rng.git)




find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


set(CMAKE_CXX_FLAGS "-O2 -funroll-loops -fopenmp -std=c++17 -Wextra -march=native -lpthread")



#old attempts
#set(CMAKE_CXX_FLAGS "-O3 -ltbb -fopenmp=libomp -std=c++17 -Wextra -march=native -lpthread -DUSETBB")


#attempt at openmp with clang
#set(CMAKE_CXX_COMPILER clang++)
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp=libomp")


## include necessary modules
#
include(FetchContent)
include(CheckLanguage)

#find_package(TBB REQUIRED)

#p2rng is the fastest solution I found for rng
# find p2rng, if not installed fetch it...
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




add_executable(simchr main.cpp rando_xo.hpp rando_mp.hpp pdqsort.h)
target_link_libraries(simchr PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(simchr PRIVATE p2rng::openmp)
#target_link_libraries(simchr PRIVATE TBB::tbb)
target_compile_features(simchr PRIVATE cxx_std_11)



