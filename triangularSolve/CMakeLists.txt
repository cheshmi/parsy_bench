cmake_minimum_required(VERSION 2.7)
project(TriangOpenMP)


message(is:"$ENV{MKLROOT}")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -std=c++11 -march=x86-64 -O3   -mavx2" )

include_directories(${MKL_INC})
link_directories(${MKL_LIB})

set(SOURCE_FILES main.cpp Triangular_BCSC.h Inspection_Level.h Triangular_CSC.h
        BLAS.h Inspection_Block.h Coarse_Inspection_Level.h)
add_executable(TriangOpenMP ${SOURCE_FILES})

target_link_libraries(TriangOpenMP

        #MKL BLAS
        mkl_intel_lp64
        mkl_sequential
        mkl_core
        mkl_intel_ilp64
        mkl_gnu_thread
        pthread
        )