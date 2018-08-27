cmake_minimum_required(VERSION 2.7)

set(MKL_INC "$ENV{MKLROOT}/include")
set(MKL_LIB "$ENV{MKLROOT}/lib/intel64")
set(MMM "$ENV{MKLROOT}/../../../lib/intel64")

set(SUITE_INC "$ENV{SUITEROOT}/include")
set(SUITE_LIB "$ENV{SUITEROOT}/lib")

set(METIS_INC "$ENV{METISROOT}/include")
set(METIS_LIB "$ENV{METISROOT}/lib")

set(SCOTCH_INC "$ENV{SCOTCHROOT}/include")
set(SCOTCH_LIB "$ENV{SCOTCHROOT}/lib/")

set(PASTIX_INC "$ENV{PASTIXROOT}/install")
set(PASTIX_LIB "$ENV{PASTIXROOT}/install/")

set(PARDISO6_LIB "/home/kazem/Documents/pardiso6/")



set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -DCHOL -DMKL -DMETIS -DVERIFY -DCOLDCACHE -std=c++11 -m64  ")

add_subdirectory(examples)
add_subdirectory(libExample)
