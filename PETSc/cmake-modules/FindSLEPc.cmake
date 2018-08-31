# Try to find SLEPc                                         
#                                                           
# Once done this will define
#                               
#  SLEPC_FOUND                - system has SLEPc                       
#  SLEPC_DIR                  - SLEPc directory                        
#  SLPEC_INCLUDES             - SLEPc include directories                
#  SLEPC_LIBRARIES            - SLEPc library (static or dynamic)
#  SLEPC_VERSION              - SLEPc Version      
#                                                           
# Usage:                                                    
#  find_package(SLEPc)                                      
#                                                           
# Setting these changes the behavior of the search          
#  SLEPC_DIR       - SLEPc directory                        
#
# Redistribution and use is allowed according to the terms of the GNU General Public License v3.0.
#

## Try to set SLEPC_DIR
if(NOT DEFINED SLEPC_DIR)
  set(SLEPC_DIR $ENV{SLEPC_DIR})
endif()

#  Find and Set SLPEC_INCLUDES
if(EXISTS "${SLEPC_DIR}/include" AND
   EXISTS "${SLEPC_DIR}/${PETSC_ARCH}/include")
 set(SLPEC_INCLUDES "${SLEPC_DIR}/include" "${SLEPC_DIR}/${PETSC_ARCH}/include")
# message(STATUS "SLPEC_INCLUDES = ${SLPEC_INCLUDES}")
else()
  message(SEND_ERROR "SLEPc includes not found")
endif()

#  Find and Set SLEPC_LIBRARIES
if(EXISTS "${SLEPC_DIR}/${PETSC_ARCH}/lib/libslepc.so")
  set(SLEPC_LIBRARIES "${SLEPC_DIR}/${PETSC_ARCH}/lib/libslepc.so")
elseif(EXISTS "${SLEPC_DIR}/${PETSC_ARCH}/lib/libslepc.a")
  set(SLEPC_LIBRARIES "${SLEPC_DIR}/${PETSC_ARCH}/lib/libslepc.a")
else()
  message(SEND_ERROR "SLEPc library could not found")
endif()

if (EXISTS "${SLEPC_DIR}/include/slepcversion.h")
  file (STRINGS "${SLEPC_DIR}/include/slepcversion.h" vstrings REGEX "#define SLEPC_VERSION_(RELEASE|MAJOR|MINOR|SUBMINOR|PATCH) ")
  foreach (line ${vstrings})
    string (REGEX REPLACE " +" ";" fields ${line}) # break line into three fields (the first is always "#define")
    list (GET fields 1 var)
    list (GET fields 2 val)
    set (${var} ${val})
    set (${var} ${val})         # Also in local scope so we have access below
  endforeach ()
  if (SLEPC_VERSION_RELEASE)
    if ($(SLPEC_VERSION_PATCH) GREATER 0)
      set (SLEPC_VERSION "${SLEPC_VERSION_MAJOR}.${SLEPC_VERSION_MINOR}.${SLEPC_VERSION_SUBMINOR}p${SLEPC_VERSION_PATCH}" CACHE INTERNAL "SLEPC version")
    else ()
      set (SLEPC_VERSION "${SLEPC_VERSION_MAJOR}.${SLEPC_VERSION_MINOR}.${SLEPC_VERSION_SUBMINOR}" CACHE INTERNAL "SLEPC version")
    endif ()
  else ()
    # make dev version compare higher than any patch level of a released version
    set (SLEPC_VERSION "${SLEPC_VERSION_MAJOR}.${SLEPC_VERSION_MINOR}.${SLEPC_VERSION_SUBMINOR}.99" CACHE INTERNAL "SLEPC version")
  endif ()
else ()
  message (SEND_ERROR "SLEPC_DIR can not be used, ${SLEPC_DIR}/include/petscversion.h does not exist")
endif ()


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SLEPc
  "SLEPc could not be found: be sure to set SLEPC_DIR in your environment variables"
  SLEPC_LIBRARIES SLPEC_INCLUDES SLEPC_DIR)
