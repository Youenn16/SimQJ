cmake_minimum_required(VERSION 3.15)
project(simQJ LANGUAGES C CXX VERSION 1.0.0)
 
configure_file(
    ${PROJECT_SOURCE_DIR}/version.h.in    # input  — explicit full path
    ${PROJECT_SOURCE_DIR}/include/version.h       # output — in build directory
)

  
 

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# ── BLAS/LAPACK ───────────────────────────────────────────────────────────────
foreach(PREFIX
    /opt/homebrew/opt/openblas   # macOS Apple Silicon
    /usr/local/opt/openblas      # macOS Intel
    /usr/lib/x86_64-linux-gnu    # Debian/Ubuntu
    /usr/lib64                   # RHEL/CentOS
)
    if(EXISTS "${PREFIX}")
        set(OPENBLAS_ROOT "${PREFIX}")
        break()
    endif()
endforeach()

find_library(OPENBLAS_LIB openblas
    PATHS ${OPENBLAS_ROOT}/lib
    NO_DEFAULT_PATH
)

if(NOT OPENBLAS_LIB)
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    find_library(LAPACKE_LIB lapacke REQUIRED)
    add_library(LAPACKE::LAPACKE IMPORTED INTERFACE)
    set_target_properties(LAPACKE::LAPACKE PROPERTIES
        INTERFACE_LINK_LIBRARIES "${LAPACKE_LIB}"
    )
else()
    # OpenBLAS includes LAPACKE — same lib, just need the header path
    add_library(BLAS::BLAS     IMPORTED INTERFACE)
    add_library(LAPACK::LAPACK IMPORTED INTERFACE)
    add_library(LAPACKE::LAPACKE IMPORTED INTERFACE)
    foreach(T BLAS::BLAS LAPACK::LAPACK LAPACKE::LAPACKE)
        set_target_properties(${T} PROPERTIES
            INTERFACE_LINK_LIBRARIES      "${OPENBLAS_LIB}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENBLAS_ROOT}/include"
        )
    endforeach()
endif()

# ── Core ──────────────────────────────────────────────────────────────────────
add_library(core
    src/sim_math.cpp
    src/output_data.cpp
    src/simulation.cpp
    src/models/ising_model.cpp
)
 

target_include_directories(core PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_compile_options(core PRIVATE -O3 -march=native -funroll-loops -ffast-math )
# target_compile_options(core PRIVATE -O3 -march=native -funroll-loops -ffast-math -g -fsanitize=address)
# target_link_libraries(core PRIVATE -fsanitize=address)

target_link_libraries(core PUBLIC BLAS::BLAS LAPACK::LAPACK LAPACKE::LAPACKE)

# ── Executable ────────────────────────────────────────────────────────────────
add_executable(simQJ src/main.cpp)
target_link_libraries(simQJ PRIVATE core)

# ── Tests ─────────────────────────────────────────────────────────────────────
include(FetchContent)
FetchContent_Declare(Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.5.2
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(Catch2)
include(${Catch2_SOURCE_DIR}/extras/Catch.cmake)
enable_testing()
add_subdirectory(tests)