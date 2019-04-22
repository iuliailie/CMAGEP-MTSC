cmake_minimum_required(VERSION 3.10)
project(multiclassCMAGEP)

set(CMAKE_CXX_STANDARD 11)
add_definitions(-std=c++11)
add_definitions(-w)

add_executable(multiclassCMAGEP
        Chromosome.cpp
        Chromosome.h
        FitnessFunction.cpp
        FitnessFunction.h
        Functions.cpp
        Functions.h
        Gene.cpp
        Gene.h
        GEP.cpp
        GEP.h
        Global.cpp
        Global.h
        main.cpp
        sorted_vector.h
        ${SOURCE_FILES})

find_package(PythonLibs REQUIRED)
include_directories(/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7)

target_link_libraries(multiclassCMAGEP python2.7)