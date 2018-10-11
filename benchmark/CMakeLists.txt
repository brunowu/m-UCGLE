include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#Add Trilinos information to the include and link lines
include_directories(${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} )
link_directories(${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS} )

#set trilinos libraries to link (LINK_LIBRARIES)
set(LINK_LIBRARIES ${Epetra_LIBRARIES} ${Anasazi_LIBRARIES} ${Amesos_LIBRARIES})

#add executable
add_executable(blockKrylovSchur blockKrylovSchur.cpp)
target_link_libraries(blockKrylovSchur ${LINK_LIBRARIES})


#set Belos libraries to link (LINK_LIBRARIES_GMRES)
set(LINK_LIBRARIES_GMRES ${Tpetra_LIBRARIES} ${Teuchos_LIBRARIES} ${Belos_LIBRARIES})

add_executable(bgmres.exe blockGMRES.cpp)

target_link_libraries(bgmres.exe ${LINK_LIBRARIES_GMRES})

INCLUDE(CPack)
