cmake_minimum_required(VERSION 3.10)
project(u1)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#find_package(Eigen3 REQUIRED) # Décommente si besoin

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

add_executable(gauge_hb src/bin/gauge_hb.cpp
    src/field/field.cpp
    src/geometry/geometry.cpp
    src/obs/observables.cpp
    src/heatbath/heatbath.cpp
    src/io/io.cpp
)
if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    target_compile_options(gauge_hb PRIVATE
            -O3
            -xHost     
            -fp-model precise
            -Wall -Wextra
    )
else()
    target_compile_options(gauge_hb PRIVATE
            -O3
            -fopenmp
            -march=native
            -Wall -Wextra -Wpedantic
    )
endif()



add_executable(gauge_ecmc src/bin/gauge_ecmc.cpp
    src/field/field.cpp
    src/geometry/geometry.cpp
    src/obs/observables.cpp
    src/ecmc/ecmc.cpp
    src/io/io.cpp
)
if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    target_compile_options(gauge_ecmc PRIVATE
            -O3
            -xHost     
            -fp-model precise
            -Wall -Wextra
    )
else()
    target_compile_options(gauge_ecmc PRIVATE
            -O3
            -fopenmp
            -march=native
            -Wall -Wextra -Wpedantic
    )
endif()
