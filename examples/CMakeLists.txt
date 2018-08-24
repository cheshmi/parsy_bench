cmake_minimum_required(VERSION 2.7)



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


