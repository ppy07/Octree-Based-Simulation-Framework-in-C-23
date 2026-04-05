cmake_minimum_required( VERSION 3.31 )
project( oktal LANGUAGES CXX )
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

set( CMAKE_CXX_STANDARD 23 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_COMPILE_WARNING_AS_ERROR ON )
add_compile_options(-Wall -Wextra -Wconversion)

# ----------------------------
# AdvPT helper libraries FIRST
# ----------------------------
option(ADVPT_ENABLE_HTGFILE "" ON)

include( FetchContent )
message( STATUS "Downloading AdvPT helper libraries..." )

FetchContent_Declare(
    advpt-helpers
    URL https://teaching.pages.i10git.cs.fau.de/advpt-exercise-ws2526/assets/advpt-helpers.tar.gz
)
FetchContent_MakeAvailable( advpt-helpers )

# ----------------------------
# Oktal library
# ----------------------------
add_library( oktal )

target_include_directories( oktal
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries( oktal
  PUBLIC
        advpt::htgfile
        advpt::mdspan
)

# ----------------------------
# Sources (compile AFTER link)
# ----------------------------
add_subdirectory( src )


# ----------------------------
# Tests
# ----------------------------
option( OKTAL_BUILD_TESTSUITE "Build the Oktal test suite" ON )

if( OKTAL_BUILD_TESTSUITE )
    include( CTest )
    add_subdirectory( tests )
endif()
