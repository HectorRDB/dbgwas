################################################################################
#  TOOLS
################################################################################
MACRO(SUBDIRLIST result curdir)
    FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
    SET (dirlist "")
    FOREACH(child ${children})
      IF(IS_DIRECTORY ${curdir}/${child})
          LIST(APPEND dirlist ${child})
      ENDIF()
    ENDFOREACH()
    SET(${result} ${dirlist})
ENDMACRO()

# We add the compilation options for the library
add_definitions (${gatb-core-flags})

# We add the gatb-core include directory
include_directories (${gatb-core-includes})

# We add the pstream include directory
set (PSTREAMS_SOURCE_DIR ${PROJECT_SOURCE_DIR}/pstreams/)
include_directories (${PSTREAMS_SOURCE_DIR})

# We add the whereami include directory
set (WHEREAMI_SOURCE_DIR ${PROJECT_SOURCE_DIR}/whereami/src)
include_directories (${WHEREAMI_SOURCE_DIR})

set(BOOST_ROOT ${CMAKE_BINARY_DIR}/../Dependencies/Source/boost_1_62_0)
set(Boost_USE_STATIC_LIBS   ON) #This makes boost plug the libraries we specify as static libraries (not dynamic ones - the user does not even need to have boost installed in their pc)
FIND_PACKAGE(Boost 1.62.0 COMPONENTS filesystem system regex serialization REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# We add the path for extra libraries
link_directories (${gatb-core-extra-libraries-path})
link_directories(${Boost_LIBRARY_DIR})

# We get the list of subdirectories names
SUBDIRLIST(PROGRAMS ${PROJECT_SOURCE_DIR}/tools)

# We create a target for each
FOREACH (program ${PROGRAMS})
  set (PROGRAM_SOURCE_DIR ${PROJECT_SOURCE_DIR}/tools/${program}/src)
  include_directories (${PROGRAM_SOURCE_DIR})
  file (GLOB_RECURSE  ProjectFiles  ${PROGRAM_SOURCE_DIR}/*.cpp)
  file (GLOB_RECURSE  WhereamiFiles  ${WHEREAMI_SOURCE_DIR}/*.c)
  add_executable(${program} ${ProjectFiles} ${WhereamiFiles})
  target_link_libraries(${program} ${gatb-core-libraries} ${Boost_LIBRARIES} -static-libgcc -static-libstdc++)
ENDFOREACH (program)

################################################################################
#  INSTALLATION 
################################################################################
install (TARGETS ${PROGRAMS} DESTINATION bin COMPONENT precompiled)
