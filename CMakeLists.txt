cmake_minimum_required(VERSION 2.7)

set(MKL_INC "$ENV{MKLROOT}/include")
set(MKL_LIB "$ENV{MKLROOT}/lib/intel64")
set(MMM "$ENV{MKLROOT}/../../../lib/intel64")
#set(MKL_INC "/home/kazem/programs/intel/compilers_and_libraries_2018.2.199/linux/mkl/include")
#set(MKL_LIB "/home/kazem/programs/intel/compilers_and_libraries_2018.2.199/linux/mkl/lib/intel64")
#set(MMM "/home/kazem/programs/intel/compilers_and_libraries_2018.2.199/linux/compiler/lib/intel64/")

#set(MKL_INC "/home/kazem/programs/intel/compilers_and_libraries_2017.4.196/linux/mkl/include")
#set(MKL_LIB "/home/kazem/programs/intel/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64")
#set(MMM "/home/kazem/programs/intel/compilers_and_libraries_2017.4.196/linux/compiler/lib/intel64/")



set(SUITE_INC "$ENV{SUITEROOT}/include")
set(SUITE_LIB "$ENV{SUITEROOT}/lib")
#set(SUITE_INC "/home/kazem/programs/SuiteSparse/include")
#set(SUITE_LIB "/home/kazem/programs/SuiteSparse/lib")
#set(SUITE_INC "/work/04218/tg835174/stampede2/programs/SuiteSparse/include")
#set(SUITE_LIB "/work/04218/tg835174/stampede2/programs/SuiteSparse/lib")


set(METIS_INC "$ENV{METISROOT}/include")
set(METIS_LIB "$ENV{METISROOT}/lib")
#set(METIS_INC "/home/kazem/programs/metis-5.1.0/include")
#set(METIS_LIB "/home/kazem/programs/metis-5.1.0/build/Linux-x86_64/libmetis")
#set(METIS_INC "/work/04218/tg835174/stampede2/programs/metis-5.1.0/include")
#set(METIS_LIB "/work/04218/tg835174/stampede2/programs/metis-5.1.0/build/Linux-x86_64/libmetis")
message("Here we go ${SUITE_LIB}")

set(SCOTCH_INC "$ENV{SCOTCHROOT}/include")
set(SCOTCH_LIB "$ENV{SCOTCHROOT}/lib/")
#set(SCOTCH_INC "/home/kazem/programs/scotch_6.0.4/include")
#set(SCOTCH_LIB "/home/kazem/programs/scotch_6.0.4/lib")

set(PASTIX_INC "$ENV{PASTIXROOT}/install")
set(PASTIX_LIB "$ENV{PASTIXROOT}/install/")
#set(PASTIX_INC "/home/kazem/programs/pastix_5.2.3/install/include")
#set(PASTIX_LIB "/home/kazem/programs/pastix_5.2.3/install/lib")

set(PARDISO6_LIB "/home/kazem/Documents/pardiso6/")



set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -DCHOL -DMKL -DMETIS -DVERIFY -DCOLDCACHE -std=c++11 -m64  ")

#add_subdirectory(triangularSolve)
#add_subdirectory(cholesky)
add_subdirectory(examples)
add_subdirectory(libExample)
