
set(OTFFT_GEN_HEADERS otfft_gen_delete.h otfft_gen_fwd.h otfft_gen_fwd0.h
                      otfft_gen_fwdu.h otfft_gen_inv.h otfft_gen_invn.h
                      otfft_gen_invu.h otfft_gen_new.h otfft_gen_setup.h)

add_executable(ffttune ffttune.cpp)
add_executable(rewrite rewrite.cpp)

target_link_libraries(ffttune ${OpenMP_CXX_LIB_NAMES})

add_custom_command(OUTPUT ${OTFFT_GEN_HEADERS}
                   COMMAND ffttune
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                   USES_TERMINAL
                   COMMENT "Tuning the library")

add_custom_target(run_ffttune DEPENDS ffttune)

add_library(otfft SHARED otfft.cpp otfft_c.cpp  ${OTFFT_GEN_HEADERS})

if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set_source_files_properties(otfft_c.cpp PROPERTIES COMPILE_FLAGS "/Qstd=c99")
endif ()

include(GenerateExportHeader)
generate_export_header(otfft
                       BASE_NAME otfft
                       EXPORT_MACRO_NAME OTFFT_EXPORT
                       EXPORT_FILE_NAME otfft_exports.h
                       STATIC_DEFINE SHARED_EXPORTS_BUILT_AS_STATIC)
target_link_libraries(otfft ${OpenMP_CXX_LIB_NAMES})
target_compile_definitions(otfft PRIVATE OTFFT_LIBRARY_BUILD)
add_dependencies(otfft ffttune run_ffttune)
target_include_directories(otfft BEFORE
                           PRIVATE ${CMAKE_CURRENT_BINARY_DIR})