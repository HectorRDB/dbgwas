################################################################################
#  TOOLS
################################################################################




################################################################################
# SeqAn related CMake code.
################################################################################
# ----------------------------------------------------------------------------
# SeqAn Dependencies
# ----------------------------------------------------------------------------

# Search for zlib as a dependency for SeqAn.
find_package (ZLIB)

# Load the SeqAn module and fail if not found.
set (CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH} ${PROJECT_SOURCE_DIR}/seqan/util/cmake")
set (SEQAN_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/seqan/include")
find_package (SeqAn REQUIRED)


# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")
################################################################################
# SeqAn related CMake code. - END
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

# We add the BOOST include directory
set(Boost_USE_STATIC_LIBS   ON) #This makes boost plug the libraries we specify as static libraries (not dynamic ones - the user does not even need to have boost installed in their pc)
FIND_PACKAGE(Boost COMPONENTS filesystem system REQUIRED)
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
  target_link_libraries(${program} ${gatb-core-libraries} ${Boost_LIBRARIES} -static-libgcc -static-libstdc++ ${SEQAN_LIBRARIES})
ENDFOREACH (program)

################################################################################
#  INSTALLATION 
################################################################################
install (TARGETS ${PROGRAMS} DESTINATION bin COMPONENT precompiled)
