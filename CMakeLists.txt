# MON 1 Aug 2016 17:02 PM CET
# Didier Benoit,  didier.benoit13@gmail.com
# Thibaut Merlin, thibaut.merlin@gmail.com
# CMAKE compiling the CASToR reconstruction software

#---------
# Checking the CMAKE version
cmake_minimum_required( VERSION 2.8 FATAL_ERROR )

#---------
# Set the build type
set( CMAKE_BUILD_TYPE "Release" CACHE STRING
  "Choose the type of build, options are: Debug Release"
)

#---------
# Defining the project
project( CASToR )
set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake-modules
  ${CMAKE_MODULE_PATH}
)

#---------
# Defining a list of option
option( CASToR_64bits "Build CASToR on 64 bits architecture" YES )
option( CASToR_VERBOSE "Print lot of infos from CASToR" NO )
option( CASToR_DEBUG "Enable many additionnal checks in primitive functions" NO )
option( CASToR_BUILD_GATE_UTILITIES "Build GATE tools for CASToR (require ROOT library)" NO )
option( CASToR_BUILD_SAMPLE_UTILITIES "Build Example codes (require ROOT library)" NO )
option( CASToR_BUILD_DEVEL "Build devel tools. Must be left unchecked" NO )
option( CASToR_ROOT "Build CASToR with ROOT library" NO )
option( CASToR_SIMD "Build CASToR with SIMD OPTION" NO )
option( CASToR_OMP "Build CASToR with OMP" NO )
option( CASToR_MPI "Build CASToR with MPI" NO )
# For testing purposes
#option( CASToR_DEFORMATION_LEDESMA "Build CASToR with LEDESMA" NO )
#option( CASToR_ELASTIX "Build CASToR with ELASTIX" NO )


#---------
# Find the config file CASToR config
find_path( CASTOR_CONFIG_DIR oCASToRConfig.hh.in
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/include/algorithm
  ${PROJECT_SOURCE_DIR}/include/analytic_simulator
  ${PROJECT_SOURCE_DIR}/include/datafile
  ${PROJECT_SOURCE_DIR}/include/dynamic
  ${PROJECT_SOURCE_DIR}/include/image
  ${PROJECT_SOURCE_DIR}/include/management
  ${PROJECT_SOURCE_DIR}/include/optimizer
  ${PROJECT_SOURCE_DIR}/include/projector
  ${PROJECT_SOURCE_DIR}/include/scanner
  ${PROJECT_SOURCE_DIR}/include/utils
)

#---------
# Defining configuration files
configure_file(
  "${CASTOR_CONFIG_DIR}/oCASToRConfig.hh.in"
  "${PROJECT_BINARY_DIR}/oCASToRConfig.hh"
  @ONLY
)

#---------
# Setting the CASTOR_CONFIG variable
set( CASTOR_CONFIG ${PROJECT_SOURCE_DIR}/config CACHE PATH
  "Path to the CASTOR config repository"
)

#---------
# Setting position of the CASToR utilitaries
set( TOOLKITS_SOURCE_DIR ${PROJECT_SOURCE_DIR}/toolkits )

# To avoid conflict between the makefile from CMake and the makefile provided
# by CASToR
add_definitions( -DCASToR_USE_CMAKE )

########################################################################
# Provide add_compile_options() when not available
########################################################################
if( CMAKE_VERSION VERSION_LESS "2.8.12" )
    function( add_compile_options )
        add_definitions( ${ARGN} )
    endfunction( add_compile_options )
endif()


#---------
# Different option for Windows or Unix (Linux/MacOS)
if( WIN32 )

  #---------
  # Delete some flags by default provided by CMAKE
  #---------
  string( REPLACE "/W3" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} )
  string( REPLACE "/O2" "" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} )
  string( REPLACE "/DNDEBUG" "" CMAKE_CXX_FLAGS_RELEASE
    ${CMAKE_CXX_FLAGS_RELEASE}
  )
  string( REPLACE "/MD" "" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} )
  string( REPLACE "/O2" "" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} )

  # Visual studio flags by default /Ox
  add_compile_options( /MP )

  #---------
  # Checking for OpenMP
  if( CASToR_OMP )
    find_package( OpenMP REQUIRED )
    set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
    set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}
      ${OpenMP_EXE_LINKER_FLAGS}"
    )
    add_definitions( /DCASTOR_OMP )
  endif()

  #---------
  # Checking for MPI
  if( CASToR_MPI )
    set( MPIEXEC "C:/Program Files/Microsoft MPI/Bin/mpiexec.exe" )
    set( MPI_CXX_INCLUDE_PATH
      "C:/Program Files (x86)/Microsoft SDKs/MPI/Include"
    )
    if( "${CMAKE_SIZEOF_VOID_P}" EQUAL 8 )
      set( MS_MPI_ARCH_DIR x64 )
    else()
      set( MS_MPI_ARCH_DIR x86 )
    endif()
    set( MPI_CXX_LIBRARIES
      "C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/${MS_MPI_ARCH_DIR}/msmpi.lib"
    )
    set( MPI_C_INCLUDE_PATH "${MPI_CXX_INCLUDE_PATH}" )
    set( MPI_C_LIBRARIES "${MPI_CXX_LIBRARIES}" )
    set( MPIEXEC_NUMPROC_FLAG "-np" CACHE STRINGS
      "Flag used by MPI to specify the number of processes for MPIEXEC; the next option will be the number of processes."
    )

    include_directories( ${MPI_CXX_INCLUDE_PATH} )
    add_definitions( /DCASTOR_MPI )
  endif()

  #---------
  # Checking verbosity
  if( CASToR_VERBOSE )
    add_compile_options( /W4 )
    add_definitions( /DCASTOR_VERBOSE )
  endif()

  #---------
  # Checking debug
  if( CASToR_DEBUG )
    add_compile_options( /W4 )
    add_definitions( /DCASTOR_DEBUG )
  endif()
  
  #---------
  # Checking the mode Release/Debug
  if( CMAKE_BUILD_TYPE MATCHES Release ) # Release mode
    # Add release flags
    add_compile_options( /Ox /O2 )

  else() # Debug mode
     # Add debug flags
      add_compile_options( /Od )
      add_definitions( /DCASTOR_DEBUG )
  endif()
  
 

    #---------
  # Checking for ROOT
  if( CASToR_ROOT )
    find_package( ROOT REQUIRED )
    include_directories( ${ROOT_INCLUDE_DIR} )
    set( LINK_DIRECTORIES ${ROOT_LIBRARY_DIR} )
    add_definitions( -DCASTOR_ROOT )
  endif()


  #---------
  # GATE utilities
  if( CASToR_BUILD_GATE_UTILITIES )
    if(CASToR_ROOT)
    else()
      message( "CASToR_ROOT switched to ON (ROOT required for GATE utilities) ${CMAKE_MODULE_PATH}" )
      set( ROOT_INCLUDE_DIR $ENV{ROOT_INCLUDE_DIR} )
      set( ROOT_LIBRARY_DIR $ENV{ROOT_LIBRARY_DIR} )
      include_directories( $ENV{ROOT_INCLUDE_DIR} )
      set( LINK_DIRECTORIES ${ROOT_LIBRARY_DIR} )
      add_definitions( /DCASTOR_ROOT )
      set( CASToR_ROOT ON CACHE BOOL "Force ROOT linking" FORCE)
    endif()
  endif()
  

  #---------
  # Checking for ROOT
  if( CASToR_BUILD_SAMPLE_UTILITIES )
    if(CASToR_ROOT)
    else()
      message( "CASToR_ROOT switched to ON (ROOT required for sample code utilities)" )
      set( ROOT_INCLUDE_DIR $ENV{ROOT_INCLUDE_DIR} )
      set( ROOT_LIBRARY_DIR $ENV{ROOT_LIBRARY_DIR} )
      include_directories( $ENV{ROOT_INCLUDE_DIR} )
      set( LINK_DIRECTORIES ${ROOT_LIBRARY_DIR} )
      add_definitions( /DCASTOR_ROOT )
      set( CASToR_ROOT ON CACHE BOOL "Force ROOT linking" FORCE)
    endif()
  endif()
  
  # Checking for both ROOT & 64bits enabled
  # Force to 32 bits in this case, as ROOT can only be 32 bits on WIN
  # https://root-forum.cern.ch/t/native-windows-64-bit-build/9102/27
  # Display warning message
  
  #---------
  # Checking the architecture 64/32 bits
  if( CASToR_64bits )
    if( CASToR_ROOT )
      message( "Architecture set to 32 bits as ROOT can only be built as 32 bits !" )
      message( "The generator may have to be manually changed to a 32 bits version !" )
      set( CASToR_64bits OFF CACHE BOOL "" FORCE)
    else()
      add_definitions( /D_64 )
    endif()
  endif()
  

  
# Unix fork (Linux/MacOS)
else()

  #---------
  # Compiler-specific C++11 activation.
  if( "${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" )
    execute_process(
      COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION
    )
    if( NOT( GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7 ) )
      message( FATAL_ERROR "${PROJECT_NAME} requires g++ 4.7 or greater." )
    endif()
    add_compile_options( -std=c++11 )
  elseif( "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" )
    add_compile_options( -stdlib=libc++ )
  endif ()

  #---------
  # Checking the architecture 64/32 bits
  if( CASToR_64bits )
    add_compile_options( -m64 )
    add_definitions( -D_64 )
  else()
    add_compile_options( -m32 )
  endif()

  #---------
  # Checking verbosity
  if( CASToR_VERBOSE )
    add_compile_options( -Wall -pedantic -Wextra )
    add_definitions( -DCASTOR_VERBOSE )
    if( CASToR_SIMD )
      add_compile_options( -ftree-vectorizer-verbose=2 )
    endif()
  endif()

  #---------
  # Checking debug
  if( CASToR_DEBUG )
    add_compile_options( -Wall -pedantic -Wextra )
    add_definitions( -DCASTOR_DEBUG )
    if( CASToR_SIMD )
      add_compile_options( -ftree-vectorizer-verbose=2 )
    endif()
  endif()
  
  #---------
  # Checking for OpenMP
  if( CASToR_OMP )
    find_package( OpenMP REQUIRED )
    add_compile_options( ${OpenMP_CXX_FLAGS} )
    add_definitions( -DCASTOR_OMP )
  endif()

  #---------
  # Checking for MPI
  if( CASToR_MPI )
    find_package( MPI REQUIRED )
    include_directories( ${MPI_INCLUDE_PATH} )
    add_definitions( -DCASTOR_MPI )
  endif()

  if( CASToR_BUILD_GATE_UTILITIES )
    if(CASToR_ROOT)
    else()
      message( "CASToR_ROOT switched to ON (ROOT required for GATE utilities)" )
      set( ROOT_INCLUDE_DIR $ENV{ROOT_INCLUDE_DIR} )
      set( ROOT_LIBRARY_DIR $ENV{ROOT_LIBRARY_DIR} )
      include_directories( $ENV{ROOT_INCLUDE_DIR} )
      set( LINK_DIRECTORIES ${ROOT_LIBRARY_DIR} )
      add_definitions( /DCASTOR_ROOT )
      set( CASToR_ROOT ON CACHE BOOL "Force ROOT linking" FORCE)
    endif()
  endif()
  

  #---------
  # Checking for ROOT
  if( CASToR_BUILD_SAMPLE_UTILITIES )
    if(CASToR_ROOT)
    else()
      message( "CASToR_ROOT switched to ON (ROOT required for sample code utilities)" )
      set( ROOT_INCLUDE_DIR $ENV{ROOT_INCLUDE_DIR} )
      set( ROOT_LIBRARY_DIR $ENV{ROOT_LIBRARY_DIR} )
      include_directories( $ENV{ROOT_INCLUDE_DIR} )
      set( LINK_DIRECTORIES ${ROOT_LIBRARY_DIR} )
      add_definitions( /DCASTOR_ROOT )
      set( CASToR_ROOT ON CACHE BOOL "Force ROOT linking" FORCE)
    endif()
  endif()
  
  #---------
  # Checking for ROOT
  if( CASToR_ROOT )
    find_package( ROOT REQUIRED )
    include_directories( ${ROOT_INCLUDE_DIR} )
    set( LINK_DIRECTORIES ${ROOT_LIBRARY_DIR} )
    add_definitions( -DCASTOR_ROOT )
  endif()

  #--------- // LEDESMA deformation library. Not planned to be included in Castor, just for private testing.
  # Checking for LEDESMA 
  #if( CASToR_DEFORMATION_LEDESMA )
  #  include_directories( ${PROJECT_SOURCE_DIR}/include/image/lib )
  #  set( LINK_DIRECTORIES ${PROJECT_SOURCE_DIR}/include/image/lib -lElasticRegistration)
  #  set( DEFORMATION_LEDESMA_LIBRARIES -L${PROJECT_SOURCE_DIR}/include/image/lib\ -lElasticRegistration)
  #  add_definitions( -DCASTOR_DEFORMATION_LEDESMA )
  #endif()
  
  #---------
  # Checking for ELASTIX
  #if( CASToR_ELASTIX )
  #  set( ELASTIX_BUILD_DIR "" CACHE PATH "Path to elastix build folder" )
  #  set( ELASTIX_USE_FILE ${ELASTIX_BUILD_DIR}/UseElastix.cmake )

  #  option( USE_ITK_FOR_ELASTIX "Use ITK (required for elastix)" ON )
  #  find_package( ITK )
  #  if( ITK_FOUND )
  #    include( ${ITK_USE_FILE} )
  #  else( ITK_FOUND )
  #    message( FATAL_ERROR "Cannot build CASToR using ELASTIX without ITK. Please set ITK_DIR." )
  #  endif( ITK_FOUND )

  #  if( EXISTS ${ELASTIX_USE_FILE} )
  #    include( ${ELASTIX_USE_FILE} )
  #    link_libraries( param )
  #    link_libraries( elastix )
  #    link_libraries( transformix )
  #    add_definitions( -DCASTOR_ELASTIX )
  #  endif()

  #endif()

  #---------
  # Checking the mode Release/Debug
  if( CMAKE_BUILD_TYPE MATCHES Release ) # Release mode
    # Add release flags
    add_compile_options( -O3 )

    # Add SIMD flags
    if( CASToR_SIMD )
      add_compile_options(
        -ftree-vectorize -msse2 -fassociative-math
      )
    endif()
  else() # Debug mode
     # Add debug flags
      add_compile_options( -O0 -g3 -fno-inline -pg )
      add_definitions( -DCASTOR_DEBUG )
  endif()

endif()

#---------
# Setup include directory for this project
include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${PROJECT_SOURCE_DIR}/include/algorithm
  ${PROJECT_SOURCE_DIR}/include/analytic_simulator
  ${PROJECT_SOURCE_DIR}/include/datafile
  ${PROJECT_SOURCE_DIR}/include/dynamic
  ${PROJECT_SOURCE_DIR}/include/image
  ${PROJECT_SOURCE_DIR}/include/management
  ${PROJECT_SOURCE_DIR}/include/optimizer
  ${PROJECT_SOURCE_DIR}/include/projector
  ${PROJECT_SOURCE_DIR}/include/scanner
  ${PROJECT_SOURCE_DIR}/include/utils
)

#---------
# Locate sources and headers for CASToR
file( GLOB_RECURSE sources ${PROJECT_SOURCE_DIR}/src/*.cc )
file( GLOB_RECURSE headers ${PROJECT_SOURCE_DIR}/include/*.hh )

#---------
# Build the executables
add_executable( castor-recon castor-recon.cc ${sources} ${headers} )
if( CASToR_BUILD_DEVEL)
  add_executable( castor-proj castor-proj.cc ${sources} ${headers} )
endif()

# Build GATE tools for CASToR 
if( CASToR_BUILD_GATE_UTILITIES )
  add_executable( castor-GATERootToCastor ${TOOLKITS_SOURCE_DIR}/castor-GATERootToCastor.cc ${sources} ${headers} )
  add_executable( castor-GATEMacToGeom    ${TOOLKITS_SOURCE_DIR}/castor-GATEMacToGeom.cc ${sources}
  ${headers} )
endif()

# Build GATE tools for CASToR 
if( CASToR_BUILD_SAMPLE_UTILITIES )
  add_executable( castor-PetScannerLutEx      ${TOOLKITS_SOURCE_DIR}/castor-PetScannerLutEx.cc ${sources} ${headers} )
  add_executable( castor-datafileConversionEx ${TOOLKITS_SOURCE_DIR}/castor-datafileConversionEx.cc ${sources} ${headers} )
  # added by chengaoyu2013@gmail.com on 2018/09/03
  add_executable( cgy-GATESinoToCastor    ${TOOLKITS_SOURCE_DIR}/cgy-GATESinoToCastor.cc ${sources} ${headers})
endif()
  
if( WIN32 )
  target_link_libraries( castor-recon ${ROOT_LIBRARIES}
    ${DEFORMATION_LEDESMA_LIBRARIES} ${ITK_LIBRARIES} ${MPI_CXX_LIBRARIES}
  )
  if( CASToR_BUILD_DEVEL )
    target_link_libraries( castor-proj ${ROOT_LIBRARIES}
      ${DEFORMATION_LEDESMA_LIBRARIES} ${ITK_LIBRARIES} ${MPI_CXX_LIBRARIES}
    )
  endif()
  if( CASToR_BUILD_GATE_UTILITIES )
    target_link_libraries( castor-GATERootToCastor ${ROOT_LIBRARIES}
      ${DEFORMATION_LEDESMA_LIBRARIES} ${ITK_LIBRARIES} ${MPI_CXX_LIBRARIES}
    )
    target_link_libraries( castor-GATEMacToGeom ${ROOT_LIBRARIES}
      ${DEFORMATION_LEDESMA_LIBRARIES} ${ITK_LIBRARIES} ${MPI_CXX_LIBRARIES}
    )
  endif()
  if( CASToR_BUILD_SAMPLE_UTILITIES )
    target_link_libraries( castor-PetScannerLutEx ${ROOT_LIBRARIES}
      ${DEFORMATION_LEDESMA_LIBRARIES} ${ITK_LIBRARIES} ${MPI_CXX_LIBRARIES}
    )
    target_link_libraries( castor-datafileConversionEx ${ROOT_LIBRARIES}
      ${DEFORMATION_LEDESMA_LIBRARIES} ${ITK_LIBRARIES} ${MPI_CXX_LIBRARIES}
    )
  endif()
  
else()
  target_link_libraries( castor-recon ${ROOT_LIBRARIES} ${OpenMP_CXX_FLAGS} ${MPI_LIBRARIES}
    ${DEFORMATION_LEDESMA_LIBRARIES} ${ITK_LIBRARIES} )
  if( CASToR_BUILD_DEVEL )
    target_link_libraries( castor-proj ${ROOT_LIBRARIES} ${MPI_LIBRARIES}
      ${DEFORMATION_LEDESMA_LIBRARIES} ${ITK_LIBRARIES} ${OpenMP_CXX_FLAGS}
    )
  endif()
  if( CASToR_BUILD_GATE_UTILITIES )
    target_link_libraries( castor-GATERootToCastor ${ROOT_LIBRARIES} ${MPI_LIBRARIES}
      ${DEFORMATION_LEDESMA_LIBRARIES} ${ITK_LIBRARIES} ${OpenMP_CXX_FLAGS}
    )
    target_link_libraries( castor-GATEMacToGeom ${ROOT_LIBRARIES} ${MPI_LIBRARIES}
      ${DEFORMATION_LEDESMA_LIBRARIES} ${ITK_LIBRARIES} ${OpenMP_CXX_FLAGS}
    )
    # added by chengaoyu2013@gmail.com on 2018/09/03
    target_link_libraries( cgy-GATESinoToCastor ${ROOT_LIBRARIES} ${MPI_LIBRARIES}
      ${DEFORMATION_LEDESMA_LIBRARIES} ${ITK_LIBRARIES} ${OpenMP_CXX_FLAGS}
    )
  endif()
  if( CASToR_BUILD_SAMPLE_UTILITIES )
    target_link_libraries( castor-PetScannerLutEx ${ROOT_LIBRARIES} ${MPI_LIBRARIES}
      ${DEFORMATION_LEDESMA_LIBRARIES} ${ITK_LIBRARIES} ${OpenMP_CXX_FLAGS}
    )
    target_link_libraries( castor-datafileConversionEx ${ROOT_LIBRARIES} ${MPI_LIBRARIES}
      ${DEFORMATION_LEDESMA_LIBRARIES} ${ITK_LIBRARIES} ${OpenMP_CXX_FLAGS}
    )
  endif()
endif()

#---------
# Move the executables to the install path
install( TARGETS
  castor-recon 
  DESTINATION
  bin
)

if( CASToR_BUILD_DEVEL )
  install( TARGETS
    castor-proj
    DESTINATION
    bin
  )
endif()

if( CASToR_BUILD_GATE_UTILITIES )
  install( TARGETS
    castor-GATERootToCastor castor-GATEMacToGeom cgy-GATESinoToCastor
    DESTINATION
    bin
  )
endif()

if( CASToR_BUILD_SAMPLE_UTILITIES )
  install( TARGETS
    castor-PetScannerLutEx castor-datafileConversionEx
    DESTINATION
    bin
  )
endif()


