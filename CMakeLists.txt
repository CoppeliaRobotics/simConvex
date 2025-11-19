cmake_minimum_required(VERSION 3.15)
project(simConvex)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_MACOSX_RPATH 1)

find_package(Boost REQUIRED)

if(NOT COPPELIASIM_INCLUDE_DIR)
    if(DEFINED ENV{COPPELIASIM_ROOT_DIR})
        set(COPPELIASIM_INCLUDE_DIR $ENV{COPPELIASIM_ROOT_DIR}/programming/include)
    else()
        set(COPPELIASIM_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../include)
    endif()
endif()
list(APPEND CMAKE_MODULE_PATH ${COPPELIASIM_INCLUDE_DIR}/cmake)
find_package(CoppeliaSim 4.5.1.7 REQUIRED COMPONENTS simLib-1 simLib-scriptFunctionData-1 simStack-1 simMath-1)

coppeliasim_add_plugin(
    simConvex
    LEGACY
    SOURCES
    sourceCode/simConvex.cpp
    external/hacd/hacdGraph.cpp
    external/hacd/hacdHACD.cpp
    external/hacd/hacdICHull.cpp
    external/hacd/hacdManifoldMesh.cpp
    external/hacd/hacdMeshDecimator.cpp
    external/hacd/hacdMicroAllocator.cpp
    external/hacd/hacdRaycastMesh.cpp
    external/vhacd/src/btAlignedAllocator.cpp
    external/vhacd/src/btConvexHullComputer.cpp
    external/vhacd/src/FloatMath.cpp
    external/vhacd/src/VHACD.cpp
    external/vhacd/src/VHACD-ASYNC.cpp
    external/vhacd/src/vhacdICHull.cpp
    external/vhacd/src/vhacdManifoldMesh.cpp
    external/vhacd/src/vhacdMesh.cpp
    external/vhacd/src/vhacdRaycastMesh.cpp
    external/vhacd/src/vhacdVolume.cpp
    external/qHull/userprintf_rbox.c
    external/qHull/userprintf.c
    external/qHull/usermem.c
    external/qHull/user.c
    external/qHull/stat.c
    external/qHull/rboxlib.c
    external/qHull/random.c
    external/qHull/qset.c
    external/qHull/poly2.c
    external/qHull/poly.c
    external/qHull/merge.c
    external/qHull/mem.c
    external/qHull/libqhull.c
    external/qHull/io.c
    external/qHull/global.c
    external/qHull/geom2.c
    external/qHull/geom.c
)
target_compile_definitions(simConvex PRIVATE SIM_MATH_DOUBLE)
target_include_directories(simConvex PRIVATE external/hacd)
target_include_directories(simConvex PRIVATE external/vhacd/inc)
target_include_directories(simConvex PRIVATE external/vhacd/public)
target_include_directories(simConvex PRIVATE external/vhacd/src)
target_include_directories(simConvex PRIVATE external)
target_include_directories(simConvex PRIVATE external/qHull)
target_include_directories(simConvex PRIVATE ${COPPELIASIM_INCLUDE_DIR}/simMath)
target_include_directories(simConvex PRIVATE ${COPPELIASIM_INCLUDE_DIR}/simStack)
target_include_directories(simConvex PRIVATE ${COPPELIASIM_INCLUDE_DIR}/simLib)
coppeliasim_add_lua(lua/simConvex.lua)
coppeliasim_add_lua(lua/simConvex-ce.lua)
coppeliasim_add_lua(lua/simQHull.lua)
coppeliasim_add_addon("addOns/Convex decomposition coacd.py")
coppeliasim_add_addon("addOns/Convex decomposition hacd.lua")
coppeliasim_add_addon("addOns/Convex decomposition vhacd.lua")
coppeliasim_add_addon("addOns/Convex hull.lua")
