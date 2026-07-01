cmake_minimum_required(VERSION 3.10)
project(u1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#find_package(Eigen3 REQUIRED) # Décommente si besoin

# Cible pointant vers le dossier src/
add_executable(gauge_hb src/bin/gauge_hb.cpp
    src/field/field.cpp
    src/geometry/geometry.cpp
    src/obs/observables.cpp
    src/heatbath/heatbath.cpp
    src/io/io.cpp
)

add_executable(gauge_ecmc src/bin/gauge_ecmc.cpp
    src/field/field.cpp
    src/geometry/geometry.cpp
    src/obs/observables.cpp
    src/ecmc/ecmc.cpp
    src/io/io.cpp
)
