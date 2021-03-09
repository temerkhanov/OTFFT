cmake_policy(SET CMP0091 NEW)
cmake_minimum_required(VERSION 3.19)

project(otfft NONE)

set(DEFAULT_BUILD_TYPE "Release")

include(CheckLanguage)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(CMAKE_VERBOSE_MAKEFILE ON)

enable_language(CXX)
enable_language(C)
check_language(Fortran)

if(CMAKE_Fortran_COMPILER)
  enable_language(Fortran)
else()
  message(STATUS "No Fortran support")
endif()

include(FindOpenMP)
if (OPENMP_FOUND)
  add_compile_options(${OpenMP_CXX_FLAGS})
  add_link_options(${OpenMP_EXE_LINKER_FLAGS})
endif()

if(APPLE)
  set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "Build architectures for Mac OS X" FORCE)
endif(APPLE)

if (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    add_compile_options(/Ofast /QxHost /Qipo /Qstd=c++17 /EHs)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    add_compile_options(-Ofast -ffast-math -std=c++17)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-Ofast -march=native -ffast-math -mprefer-vector-width=256 -std=c++17)
elseif (CMAKE_CXX_COMPILER_ID MATCHES  "GNU")
    add_compile_options(-Ofast -march=native -ffast-math -std=c++17)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    add_compile_options(/Ofast /Oa /EHs /QxHost /Qunroll /Qstd=c++17 /Qipo
                        /Qcomplex-limited-range /Qansi-alias-check-
                        /Qopt-subscript-in-range /Qopt-matmul
                        /Qopt-prefetch /Qopt-prefetch-issue-excl-hint /Qopt-dynamic-align
                        /Qalign-loops /Qopt-mem-layout-trans:3 /Qipo-jobs16
                        /Qfp-speculation:fast /Qfast-transcendentals /Qimf-use-svml:true
                        /Qfinite-math-only /Ob2 /Qfnalign /Qsfalign16 /Qalign-loops /Gv)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    add_compile_options(/Ox /favor:Intel64 /Ob2 /EHsc /fp:fast /Qspectre- /arch:AVX2
                        /openmp)
else()
    message(FATAL_ERROR "Unsupported compiler")
endif ()

add_subdirectory(otfft)
#add_subdirectory(C99)

if (CMAKE_FORTRAN_COMPILER_ID)
    add_subdirectory(Fortran2003)
endif()

find_library(FFTW_LIB fftw3)

add_executable(bstcheck bstcheck.cpp )
add_executable(fftbench fftbench.cpp)
add_executable(dctcheck dctcheck.cpp)
add_executable(fftcheck fftcheck.cpp)
add_executable(fftcheck2 fftcheck2.cpp)
add_executable(rfftcheck rfftcheck.cpp)
add_executable(rfftcheck2 rfftcheck2.cpp)

target_link_libraries(fftbench otfft_static)

if (FFTW_LIB)
    add_executable(fftbench0 fftbench0.cpp)
    add_executable(fftbench1 fftbench1.cpp)
    add_executable(fftbench2 fftbench2.cpp)
    add_executable(fftbench3 fftbench3.cpp)
    target_link_libraries(fftbench0 otfft_shared fftw3 fftw3_threads pthread)
    target_link_libraries(fftbench1 otfft_shared fftw3 fftw3_threads pthread)
    target_link_libraries(fftbench2 otfft_shared fftw3 fftw3_threads pthread)
    target_link_libraries(fftbench3 otfft_shared fftw3 fftw3_threads)
endif()

target_link_libraries(bstcheck otfft_shared)
target_link_libraries(dctcheck otfft_shared)
target_link_libraries(fftcheck otfft_shared)
target_link_libraries(fftcheck2 otfft_shared)
target_link_libraries(rfftcheck otfft_shared)
target_link_libraries(rfftcheck2 otfft_shared)
