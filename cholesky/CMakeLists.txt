cmake_minimum_required(VERSION 3.1)
project(CholOpenMP)
#TODO use env instead
#set(MKL_INC "/home/kazem/programs/intel/compilers_and_libraries_2017.4.196/linux/mkl/include")
#set(MKL_LIB "/home/kazem/programs/intel/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64")
set(MMM "/home/kazem/programs/intel/compilers_and_libraries_2017.4.196/linux/compiler/lib/intel64/")

#set(SUITE_INC "/home/kazem/Downloads/SuiteSparse/include")
#set(SUITE_LIB "/home/kazem/Downloads/SuiteSparse/lib")

#set(METIS_INC "/home/kazem/programs/metis-5.1.0/include")
#set(METIS_LIB "/home/kazem/programs/metis-5.1.0/build/Linux-x86_64/libmetis")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMKL -DMETIS   -std=gnu++11  -fopenmp ")

set(SOURCE_FILES main.cpp PB_Cholesky.h MyBLAS.h Ordering.h Inspection_Block.h
        Etree.h def.h PostOrder.h DFS.h Reach.h Inspection_Prune.h ColumnCount.h
        Transpose.h Util.h LSparsity.h SparseUtils.h colCounts.h
        Inspection_BlockC.h AMD.h InspectionLevel.h parallel_PB_Cholesky.h
        parallel_PB_Cholesky_02.h parallel_PB_Cholesky_03.h InspectionLevel_02.h InspectionLevel_01.h InspectionLevel_03.h
        InspectionLevel_04.h InspectionLevel_05.h Partitioning.h Merge.h TreeUtils.h parallel_PB_Cholesky_04.h MergeWithTightSlack.h MergewithBalancedSlack.h MergewithBalancedSlackGivenLevels.h parallel_PB_Cholesky_05.h InspectionLevel_06.h PostOrderSpliting.h performanceModel.h Cholesky_Perf_Modeling.h Cholesky_Perf_Comm_modeling.h)

include_directories(${MKL_INC} ${SUITE_INC} ${METIS_INC} "../common/")
link_directories(${MKL_LIB} ${MMM} ${SUITE_LIB} ${METIS_LIB})
add_executable(CholOpenMP ${SOURCE_FILES})
target_link_libraries(CholOpenMP
        #mkl_blacs_intelmpi_ilp64
        #mkl_blacs_intelmpi_lp64
        #mkl_blacs_openmpi_ilp64
    #    mkl_blacs_openmpi_lp64
#        mkl_blacs_sgimpt_ilp64
   #     mkl_blacs_sgimpt_lp64
  #      mkl_def
        #mkl_gf_lp64

 #       mkl_blas95_ilp64
 #       mkl_cdft_core

        #mkl_intel_ilp64
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
        #mkl_blacs_openmpi_lp64
        #iomp
#[[        mkl_blacs_openmpi_lp64
        mkl_blacs_openmpi_ilp64
        mkl_pgi_thread]]
        #mkl_tbb_thread
        #Libs from CHOLMOD for ordering
        amd
        cholmod
        suitesparseconfig

        metis


        #openblas
        )
