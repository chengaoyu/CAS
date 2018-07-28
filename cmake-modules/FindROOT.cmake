# - Finds ROOT instalation
# This module sets up ROOT information
# It defines:
# ROOT_FOUND          If the ROOT is found
# ROOT_INCLUDE_DIR    PATH to the include directory
# ROOT_LIBRARIES      Most common libraries
# ROOT_LIBRARY_DIR    PATH to the library directory





# root config is a bash script and not commonly executable under Windows
# ask to the user to manually find the root-config file
# CASToR v2.0 has been successfully associated to root v5.34 with MSVC 2015 (Win32)
IF ( WIN32 )

  # please set manually the path to root here
  # or set the path to the root-config file in the cmake gui
  # SET(ENV{ROOTSYS} "path/to/root")



  SET(ROOT_CONFIG_SEARCHPATH
    $ENV{ROOTSYS}/bin
    /usr/local/bin
    /opt/local/bin
    /root/bin
  )

  FIND_PROGRAM(ROOT_CONFIG_EXECUTABLE NAMES root-config PATHS
     ${ROOT_CONFIG_SEARCHPATH}
     NO_DEFAULT_PATH)



  IF (${ROOT_CONFIG_EXECUTABLE} MATCHES "ROOT_CONFIG_EXECUTABLE-NOTFOUND")
    MESSAGE( FATAL_ERROR "ROOT not installed in the searchpath and ROOTSYS is not set.
              Please set ROOT_CONFIG_EXECUTABLE with the root-config file
              (located on \path\to\root\bin) or add the path to your ROOT installation
              in the Macro FindROOT.cmake in the subdirectory cmake/modules.")
  ELSE (${ROOT_CONFIG_EXECUTABLE} MATCHES "ROOT_CONFIG_EXECUTABLE-NOTFOUND")
    STRING(REGEX REPLACE "(^.*)/bin/root-config" "\\1" test ${ROOT_CONFIG_EXECUTABLE})
    SET( ENV{ROOTSYS} ${test})
    set( ROOTSYS ${test})
  ENDIF (${ROOT_CONFIG_EXECUTABLE} MATCHES "ROOT_CONFIG_EXECUTABLE-NOTFOUND")

    SET(ROOT_FOUND FALSE)
    IF (ROOT_CONFIG_EXECUTABLE)
    SET(ROOT_FOUND TRUE)
    set(ROOT_INCLUDE_DIR ${ROOTSYS}/include)
    set(ROOT_LIBRARY_DIR ${ROOTSYS}/lib)
    SET(ROOT_BINARY_DIR ${ROOTSYS}/bin)
    set(ROOT_LIBRARIES -LIBPATH:${ROOT_LIBRARY_DIR} libGpad.lib libHist.lib libGraf.lib libGraf3d.lib libTree.lib libRint.lib libPostscript.lib libMatrix.lib libPhysics.lib libMathCore.lib libRIO.lib libNet.lib libThread.lib libCore.lib libCint.lib libMinuit.lib libGui.lib libSpectrum.lib)
    include(CMakeMacroParseArguments)
    FIND_PROGRAM(ROOTCINT_EXECUTABLE
      NAMES rootcint
      PATHS ${ROOT_BINARY_DIR}
      NO_DEFAULT_PATH
      )
    MESSAGE(STATUS "Found ROOT: $ENV{ROOTSYS}/bin/root (WIN32/version not identified)")
    ENDIF (ROOT_CONFIG_EXECUTABLE)

ELSE(WIN32)

  find_program(ROOT_CONFIG_EXECUTABLE root-config
  PATHS $ENV{ROOTSYS}/bin)

  if(NOT ROOT_CONFIG_EXECUTABLE)
    set(ROOT_FOUND FALSE)
  else()
    set(ROOT_FOUND TRUE)
  
    execute_process(
      COMMAND ${ROOT_CONFIG_EXECUTABLE} --prefix
      OUTPUT_VARIABLE ROOTSYS
      OUTPUT_STRIP_TRAILING_WHITESPACE)
  
    execute_process(
      COMMAND ${ROOT_CONFIG_EXECUTABLE} --version
      OUTPUT_VARIABLE ROOT_VERSION
      OUTPUT_STRIP_TRAILING_WHITESPACE)
  
    execute_process(
      COMMAND ${ROOT_CONFIG_EXECUTABLE} --incdir
      OUTPUT_VARIABLE ROOT_INCLUDE_DIR
      OUTPUT_STRIP_TRAILING_WHITESPACE)
  
    execute_process(
      COMMAND ${ROOT_CONFIG_EXECUTABLE} --libs
      OUTPUT_VARIABLE ROOT_LIBRARIES
      OUTPUT_STRIP_TRAILING_WHITESPACE)
  
    #set(ROOT_LIBRARIES ${ROOT_LIBRARIES} -lThread -lMinuit -lHtml -lVMC -lEG -lGeom -lTreePlayer -lXMLIO -lProof)
    #set(ROOT_LIBRARIES ${ROOT_LIBRARIES} -lProofPlayer -lMLP -lSpectrum -lEve -lRGL -lGed -lXMLParser -lPhysics)
    set(ROOT_LIBRARY_DIR ${ROOTSYS}/lib)
  
    # Make variables changeble to the advanced user
    mark_as_advanced(ROOT_CONFIG_EXECUTABLE)
  
    if(NOT ROOT_FIND_QUIETLY)
      message(STATUS "Found ROOT ${ROOT_VERSION} in ${ROOTSYS}")
    endif()
  endif()
ENDIF(WIN32)


#include(CMakeMacroParseArguments)
#find_program(ROOTCINT_EXECUTABLE rootcint PATHS $ENV{ROOTSYS}/bin)

#----------------------------------------------------------------------------
# function ROOT_GENERATE_DICTIONARY( dictionary
#                                    header1 header2 ...
#                                    LINKDEF linkdef1 ...
#                                    OPTIONS opt1...)
function(ROOT_GENERATE_DICTIONARY dictionary)
  CMAKE_PARSE_ARGUMENTS(ARG "" "" "LINKDEF;OPTIONS" "" ${ARGN})
  #---Get the list of header files-------------------------
  set(headerfiles)
  foreach(fp ${ARG_UNPARSED_ARGUMENTS})
    file(GLOB files ${fp})
    if(files)
      foreach(f ${files})
        if(NOT f MATCHES LinkDef)
          set(headerfiles ${headerfiles} ${f})
        endif()
      endforeach()
    else()
      set(headerfiles ${headerfiles} ${fp})
    endif()
  endforeach()
  #---Get the list of include directories------------------
  get_directory_property(incdirs INCLUDE_DIRECTORIES)
  set(includedirs)
  foreach( d ${incdirs})
   set(includedirs ${includedirs} -I${d})
  endforeach()
  #---Get LinkDef.h file------------------------------------
  set(linkdefs)
  foreach( f ${ARG_LINKDEF})
    if( IS_ABSOLUTE ${f})
      set(linkdefs ${linkdefs} ${f})
    else()
      if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/inc/${f})
        set(linkdefs ${linkdefs} ${CMAKE_CURRENT_SOURCE_DIR}/inc/${f})
      else()
        set(linkdefs ${linkdefs} ${CMAKE_CURRENT_SOURCE_DIR}/${f})
      endif()
    endif()
  endforeach()
  #---call rootcint------------------------------------------
  add_custom_command(OUTPUT ${dictionary}.cxx ${dictionary}.h
                     COMMAND ${ROOTCINT_EXECUTABLE} -cint -f  ${dictionary}.cxx
                                          -c ${ARG_OPTIONS} ${includedirs} ${headerfiles} ${linkdefs}
                     DEPENDS ${headerfiles} ${linkdefs})
endfunction()
