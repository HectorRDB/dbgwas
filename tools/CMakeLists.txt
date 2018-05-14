################################################################################
#  TOOLS
################################################################################
include(ExternalProject)

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

# Installing BOOST using CMAKE (based on https://stackoverflow.com/questions/43692976/compile-boost-with-multi-configuration-with-cmake)
#message(STATUS "Installing Boost library.")
#set(BOOST_BOOTSTRAP_COMMAND ./bootstrap.sh)
#set(BOOST_B2_COMMAND ./b2)
#set(BOOST_SRC_DIR ${PROJECT_BINARY_DIR}/boost_1_62_0)
#set(BOOST_INSTALL_DIR ${PROJECT_BINARY_DIR}/boost_1_62_0_install)

#set(BOOST_BUILD_TYPE variant=release)
#if (${CMAKE_BUILD_TYPE} MATCHES Debug)
#  set(BOOST_BUILD_TYPE variant=debug)
#endif(${CMAKE_BUILD_TYPE} MATCHES Debug)

#ExternalProject_Add(boost
#  SOURCE_DIR ${BOOST_SRC_DIR}
#  BUILD_IN_SOURCE 1
#  GIT_REPOSITORY "https://github.com/boostorg/boost"
#  GIT_TAG "4f2bdeb93a4be13ba0dc5e9f44920a2bf67191fc"
#  CONFIGURE_COMMAND ${BOOST_BOOTSTRAP_COMMAND}
#  BUILD_COMMAND ${BOOST_B2_COMMAND} headers COMMAND ${BOOST_B2_COMMAND} install
#    link=static
#    ${BOOST_BUILD_TYPE}
#    --prefix=${BOOST_INSTALL_DIR}
#    --with-filesystem
#    --with-system
#    -j2
#  INSTALL_COMMAND ""
#)

# We add the BOOST to the project
#set(BOOST_ROOT ${BOOST_INSTALL_DIR})
set(Boost_USE_STATIC_LIBS   ON) #This makes boost plug the libraries we specify as static libraries (not dynamic ones - the user does not even need to have boost installed in their pc)
FIND_PACKAGE(Boost 1.62.0 COMPONENTS filesystem system REQUIRED)
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
