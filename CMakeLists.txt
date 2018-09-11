cmake_minimum_required(VERSION 3.6)


MESSAGE("Starting configure")

#set TRILINOS_PATH of the platform
#set(TRILINOS_PATH /Users/xwu/phd_xwu/trilinos-12.12.1-Source/installation)
set(TRILINOS_PATH /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2)
#set(TRILINOS_PATH /home/xinzhewu/workstations/trilinos/install)
#set(TRILINOS_PATH /home/xinzhewu/workstations/trilinos/install-openmpi-1.10.2)
#set(TRILINOS_PATH /home/xinzhewu/workstations/trilinos/install-openmpi-1.10.2-gcc5.1.0)

SET(CMAKE_PREFIX_PATH ${TRILINOS_PATH} ${CMAKE_PREFIX_PATH})

#find trilinos
FIND_PACKAGE(Trilinos REQUIRED PATHS ${TRILINOS_PATH}/include ${TRILINOS_PATH})

#Output some info about the Trilinos installation found
IF(Trilinos_FOUND)
   MESSAGE("\nFound Trilinos!  Here are the details: ")
   MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
   MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
   MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
   MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
   MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
   MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
   MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
   MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
   # Will Tramonto work with a Shared version of Trilinos?
   MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
   #These last outputs are temporary
   MESSAGE("   Trilinos_CXX_COMPILER = ${Trilinos_CXX_COMPILER}")
   MESSAGE("   Trilinos_C_COMPILER = ${Trilinos_C_COMPILER}")
   MESSAGE("   Trilinos_Fortran_COMPILER = ${Trilinos_Fortran_COMPILER}")
   MESSAGE("   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
   MESSAGE("   Trilinos_C_COMPILER_FLAGS = ${Trilinos_C_COMPILER_FLAGS}")
   MESSAGE("   Trilinos_Fortran_COMPILER_FLAGS =
     ${Trilinos_Fortran_COMPILER_FLAGS}")
   MESSAGE("   Trilinos_LINKER = ${Trilinos_LINKER}")
   MESSAGE("   Trilinos_EXTRA_LD_FLAGS = ${Trilinos_EXTRA_LD_FLAGS}")
   MESSAGE("   Trilinos_AR = ${Trilinos_AR}")
   # Add MPI_EXEC, MPI_EXEC_MAX_NUMPROCS, MPI_EXEC_NUMPROC_FLAGS
   # These are not available in Trilinos 10.2.0, but should be in 10.4.0
   MESSAGE("End of Trilinos details\n")
ELSE()
  # We should never reach this because FIND_PACKAGE(Trilinos now uses REQUIRED
  MESSAGE(FATAL_ERROR "Could not find Trilinos 10.2 or greater!")
ENDIF()

IF(NOT CMAKE_CXX_COMPILER)
  MESSAGE("Setting CXX compiler.")
  SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})
ENDIF()

IF(NOT CMAKE_C_COMPILER)
  MESSAGE("Setting C compiler.")
  SET(CMAKE_C_COMPILER ${Trilinos_C_COMPILER})
ENDIF()

IF(NOT CMAKE_LINKER)
  MESSAGE("Setting linker.")
  SET(CMAKE_LINKER ${Trilinos_LINKER})
ENDIF()

IF(NOT CMAKE_EXE_LINKER_FLAGS)
  MESSAGE("Setting linker flags.")
  SET(CMAKE_EXE_LINKER_FLAGS ${Trilinos_EXTRA_LD_FLAGS})
ENDIF()

IF(NOT CMAKE_AR)
  MESSAGE("Setting ar.")
  SET(CMAKE_AR ${CMAKE_AR})
ENDIF()

project(trilinos_tutorial)
SET(CPACK_PACKAGE_VERSION "1.0")
set(TRILINOS_PATH "" CACHE FILEPATH "Path to Trilinos Directory")

IF(NOT CMAKE_CXX_FLAGS)
  MESSAGE("Setting CXX compiler flags.")
  SET(CMAKE_CXX_FLAGS ${Trilinos_CXX_COMPILER_FLAGS})
ENDIF()

IF(NOT CMAKE_C_FLAGS)
  MESSAGE("Setting C compiler flags.")
  SET(CMAKE_C_FLAGS ${Trilinos_C_COMPILER_FLAGS})
ENDIF()

ENABLE_TESTING()

#ADD_SUBDIRECTORY(Solvers)
#ADD_SUBDIRECTORY(Solvers/GMRES)

#ADD_SUBDIRECTORY(benchmark)

#Add Trilinos information to the include and link lines
include_directories(${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} )
link_directories(${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS} )
# /Library/Frameworks/QtCore.framework /Library/Frameworks/QtGui.framework)

#set trilinos libraries to link (LINK_LIBRARIES)
set(LINK_LIBRARIES ${Teuchos_LIBRARIES})

#set Anasazi libraries to link (LINK_LIBRARIES_ARNOLDI)
set(LINK_LIBRARIES_ARNOLDI ${Epetra_LIBRARIES} ${Anasazi_LIBRARIES} ${Amesos_LIBRARIES})

#set Belos libraries to link (LINK_LIBRARIES_GMRES)
set(LINK_LIBRARIES_GMRES ${Tpetra_LIBRARIES} ${Teuchos_LIBRARIES} ${Belos_LIBRARIES})

#add executable
add_executable(test.exe main.cpp)
add_executable(gmres.exe lsgmres.cpp)
#add_executable(gmres2.exe lsgmres.cpp)
add_executable(arnoldi.exe arnoldi.cpp)
add_executable(arnoldi2.exe arnoldi2.cpp)
add_executable(lsqr.exe lsqr.cpp)

add_executable(testlsqrBelos.exe test/testlsqrBelos.cpp)

target_link_libraries(lsqr.exe ${LINK_LIBRARIES})
target_link_libraries(arnoldi.exe ${LINK_LIBRARIES_ARNOLDI})
target_link_libraries(arnoldi2.exe ${LINK_LIBRARIES_ARNOLDI})
target_link_libraries(gmres.exe ${LINK_LIBRARIES_GMRES})
#target_link_libraries(gmres2.exe ${LINK_LIBRARIES_GMRES})
target_link_libraries(testlsqrBelos.exe ${LINK_LIBRARIES_GMRES})

list (APPEND CMAKE_MODULE_PATH "./cmake-modules")

find_package (PETSc REQUIRED)

if(PETSC_FOUND)
    message([STATUS] "PETSC is provided.")
endif()

if(PETSC_FOUND)
  find_package(SLEPc REQUIRED)
  if(SLEPC_FOUND)
     message([STATUS] "SLEPC version ${SLEPC_VERSION} is provided.")
  endif()
endif()

message ([STATUS] "CPP COMPILER ${PETSC_COMPILER}")
message ([STATUS] "PETSC INCLIUDES = ${PETSC_INCLUDES}")
message ([STATUS] "SLPEC INCLIUDES = ${SLPEC_INCLUDES}")

include_directories(${PETSC_INCLUDES} ${SLPEC_INCLUDES})

add_definitions (${PETSC_DEFINITIONS})

add_definitions (-g)

message (STATUS "PETSC_DEFINITIONS ${PETSC_DEFINITIONS}")

message (STATUS "PETSC_LIBRARIES ${PETSC_LIBRARIES}")

message (STATUS "SLPEC_LIBRARIES ${SLEPC_LIBRARIES}")

add_executable(lsqrPETSc.exe PETSc/LS/lsqrPETSc.cpp)
add_executable(arnoldiPETSc.exe PETSc/Arnoldi/arnoldiPETSc.cpp)

add_executable(testlsqrPetsc.exe test/testlsqrPetsc.cpp)

target_link_libraries(lsqrPETSc.exe ${PETSC_LIBRARIES})

target_link_libraries(testlsqrPetsc.exe ${PETSC_LIBRARIES})

target_link_libraries(arnoldiPETSc.exe ${SLEPC_LIBRARIES} ${PETSC_LIBRARIES})


include(Dart)
include(CPack)
