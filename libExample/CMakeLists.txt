cmake_minimum_required(VERSION 2.7)



project(pardiso6_cholesky)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -fopenmp ")

set(SOURCE_FILES cholesky_pardiso6.cpp)


include_directories(${MKL_INC} ${SUITE_INC} ${METIS_INC}
        ${PASTIX_INC} "../cholesky/" "../common/" )
link_directories(${MKL_LIB} ${MMM} ${SUITE_LIB} ${PARDISO6_LIB} )
add_executable(pardiso6_cholesky ${SOURCE_FILES} cholesky_pardiso6.cpp)
target_link_libraries(pardiso6_cholesky
        pardiso600-GNU720-X86-64.so

        #mkl_intel_lp64
        mkl_gf_lp64
        mkl_sequential
        iomp5
        #gomp
        mkl_core
        #mkl_rt
        #mkl_gnu_thread
        #mkl_intel_thread

        gfortran
        pthread
        m
        dl
        amd
        cholmod
        suitesparseconfig
        metis

        )