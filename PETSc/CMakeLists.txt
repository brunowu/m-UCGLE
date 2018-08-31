cmake_minimum_required (VERSION 3.6)

project(SLEPcBenchmark)

list (APPEND CMAKE_MODULE_PATH "${SLEPcBenchmark_SOURCE_DIR}/cmake-modules")

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

message ([STATUS] "PETSC COMPILER is found as ${PETSC_COMPILER}")

#add_definitions (-std=c++0x)

message ([STATUS] "CPP COMPILER ${PETSC_COMPILER}")
message ([STATUS] "PETSC INCLIUDES = ${PETSC_INCLUDES}")
message ([STATUS] "SLPEC INCLIUDES = ${SLPEC_INCLUDES}")

include_directories(${PETSC_INCLUDES} ${SLPEC_INCLUDES})

add_definitions (${PETSC_DEFINITIONS})

add_definitions (-g)

message (STATUS "PETSC_DEFINITIONS ${PETSC_DEFINITIONS}")

message (STATUS "PETSC_LIBRARIES ${PETSC_LIBRARIES}")

add_executable(lsqrPETSc.exe LS/lsqrPETSc.cpp)
add_executable(arnoldiPETSc.exe Arnoldi/arnoldiPETSc.cpp)


target_link_libraries(lsqrPETSc.exe ${PETSC_LIBRARIES} ${SLEPC_LIBRARIES})

target_link_libraries(arnoldiPETSc.exe ${PETSC_LIBRARIES} ${SLEPC_LIBRARIES})
