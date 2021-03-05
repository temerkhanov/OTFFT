project(otfft NONE)
cmake_minimum_required(VERSION 3.16)

set(DEFAULT_BUILD_TYPE "Release")

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
enable_language(Fortran OPTIONAL)

include(FindOpenMP)
if (OPENMP_FOUND)
  add_compile_options(${OpenMP_CXX_FLAGS})
  add_link_options(${OpenMP_EXE_LINKER_FLAGS})
endif()


if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-Ofast -march=native -ffast-math -mprefer-vector-width=256)
elseif (CMAKE_CXX_COMPILER_ID MATCHES  "GNU")
    add_compile_options(-Ofast -march=native -ffast-math)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    add_compile_options(/Ofast /Oa /EHs /QxHost /Qipo /Qunroll-aggressive
                        /Qcomplex-limited-range /Qalias-const /Qvec-guard-write
                        /Qopt-subscript-in-range /Qopt-matmul /Qopt-mem-layout-trans:3
                        /Qopt-prefetch /Qopt-prefetch-issue-excl-hint /Qopt-dynamic-align
                        /Qalign-loops /Qopt-zmm-usage:high /Qbranches-within-32B-boundaries
                        /Qfp-speculation:fast /Qfast-transcendentals /Qimf-use-svml:true
                        /Qfinite-math-only /Ob2 /ZI /Qfnalign /Qsfalign16 /Qalign-loops)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    add_compile_options(/Ox /favor:Intel64 /Ob2 /EHsc /fp:fast /Qspectre- /arch:AVX2
                        /openmp)
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

target_link_libraries(bstcheck otfft)
target_link_libraries(fftbench otfft)

if (FFTW_LIB)
    add_executable(fftbench0 fftbench0.cpp)
    add_executable(fftbench1 fftbench1.cpp)
    add_executable(fftbench2 fftbench2.cpp)
    add_executable(fftbench3 fftbench3.cpp)
    target_link_libraries(fftbench0 otfft fftw3 fftw3_threads pthread)
    target_link_libraries(fftbench1 otfft fftw3 fftw3_threads pthread)
    target_link_libraries(fftbench2 otfft fftw3 fftw3_threads pthread)
    target_link_libraries(fftbench3 otfft fftw3 fftw3_threads)
endif()

target_link_libraries(dctcheck otfft)
target_link_libraries(fftcheck otfft)
target_link_libraries(fftcheck2 otfft)
target_link_libraries(rfftcheck otfft)
target_link_libraries(rfftcheck2 otfft)