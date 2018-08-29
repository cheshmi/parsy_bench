cmake_minimum_required(VERSION 2.7)

project(triangularTest_chol)

include_directories(${MKL_INC} ${SUITE_INC} ${METIS_INC}
        "../cholesky/" "../triangularSolve/" "../common/")
link_directories(${MKL_LIB}  ${MMM} ${SUITE_LIB} ${METIS_LIB})

set(SOURCE_FILES triangularTest02.cpp )
add_executable(triangularTest_chol ${SOURCE_FILES})

target_link_libraries(triangularTest_chol

        #MKL BLAS
        mkl_intel_lp64
       # mkl_sequential
        mkl_core
        mkl_intel_ilp64
        mkl_gnu_thread
        pthread

        #gomp
        iomp5
        #mkl_sequential
        mkl_intel_lp64
        mkl_core
        #mkl_gnu_thread
        mkl_intel_thread
        pthread
        m
        dl
        amd
        cholmod
        suitesparseconfig
        metis
        )

project(triangularTest_DAG_chol)

include_directories(${MKL_INC} ${SUITE_INC} ${METIS_INC}
        "../cholesky/" "../triangularSolve/" "../common/")
link_directories(${MKL_LIB}  ${MMM} ${SUITE_LIB} ${METIS_LIB})

set(SOURCE_FILES triangularTest_DAG.cpp )
add_executable(triangularTest_DAG_chol ${SOURCE_FILES})

target_link_libraries(triangularTest_DAG_chol

        #MKL BLAS
        mkl_intel_lp64
       # mkl_sequential
        mkl_core
        mkl_intel_ilp64
        mkl_gnu_thread
        pthread

        #gomp
        iomp5
        #mkl_sequential
        mkl_intel_lp64
        mkl_core
        #mkl_gnu_thread
        mkl_intel_thread
        pthread
        m
        dl
        amd
        cholmod
        suitesparseconfig
        metis
        )


project(triangularTest_DAG_nonChordal)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMKL   -fopenmp ")

include_directories(${MKL_INC} ${SUITE_INC} ${METIS_INC}
        "../cholesky/" "../triangularSolve/" "../common/")
link_directories(${MKL_LIB}  ${MMM} ${SUITE_LIB} ${METIS_LIB})

set(SOURCE_FILES triangularTest_DAG_nonChordal.cpp )
add_executable(triangularTest_DAG_nonChordal ${SOURCE_FILES})

target_link_libraries(triangularTest_DAG_nonChordal

        iomp5
        #gomp
        #mkl_sequential
        mkl_intel_lp64
        mkl_core
        #mkl_gnu_thread
        mkl_intel_thread
        pthread

        m
        dl
        amd
        cholmod
        suitesparseconfig
        metis
        )


#####
project(choleskyTest)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMKL   -fopenmp ")

set(SOURCE_FILES choleskyTest01.cpp )

include_directories(${MKL_INC} ${SUITE_INC} ${METIS_INC} ${SCOTCH_INC}
        ${PASTIX_INC} "../cholesky/" "../common/" )
link_directories(${MKL_LIB} ${MMM} ${SUITE_LIB} ${METIS_LIB} ${SCOTCH_LIB})
add_executable(choleskyTest ${SOURCE_FILES})
target_link_libraries(choleskyTest
        #gomp
        iomp5
        #mkl_sequential
        mkl_intel_lp64
        mkl_core
        #mkl_gnu_thread
        mkl_intel_thread
        pthread
        m
        dl
        amd
        cholmod
        suitesparseconfig
        metis

 #       scotch
 #       scotcherrexit
 #       scotchmetis
        )


project(choleskyTest03)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMKL    -fopenmp ")

set(SOURCE_FILES choleskyTest03.cpp)

include_directories(${MKL_INC} ${SUITE_INC} ${METIS_INC}
        "../cholesky/" "../common/")
link_directories(${MKL_LIB} ${MMM} ${SUITE_LIB} ${METIS_LIB})
add_executable(choleskyTest03 ${SOURCE_FILES} choleskyTest03.cpp)
target_link_libraries(choleskyTest03
        #gomp
        iomp5
        #mkl_sequential
        mkl_intel_lp64
        mkl_core
        #mkl_gnu_thread
        mkl_intel_thread
        pthread
        m
        dl
        amd
        cholmod
        suitesparseconfig
        metis
        )





project(makingLowerHalf)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMKL  -fopenmp ")

set(SOURCE_FILES MakingLowerHalf.cpp)

include_directories(${MKL_INC} ${SUITE_INC} ${METIS_INC}
        "../cholesky/" "../common/")
#link_directories(${MKL_LIB} ${MMM} ${SUITE_LIB} ${METIS_LIB})
add_executable(makingLowerHalf ${SOURCE_FILES} MakingLowerHalf.cpp)
target_link_libraries(makingLowerHalf

        )


