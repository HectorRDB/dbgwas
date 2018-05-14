include (ExternalProject)

set_property (DIRECTORY PROPERTY EP_BASE Dependencies)

set (DEPENDENCIES)
set (EXTRA_CMAKE_ARGS)

# Use static linking to avoid issues with system-wide installations of Boost.
list (APPEND DEPENDENCIES boost_1_62_0)
ExternalProject_Add (boost_1_62_0
  URL ${PROJECT_SOURCE_DIR}/thirdparty/boost_1_62_0.tar.gz
  URL_MD5 6f4571e7c5a66ccc3323da6c24be8f05
  CONFIGURE_COMMAND ./bootstrap.sh --with-libraries=filesystem,system
  BUILD_COMMAND ./b2 link=static
  BUILD_IN_SOURCE 1
  INSTALL_COMMAND ""
  )
list (APPEND EXTRA_CMAKE_ARGS
  -DBOOST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/boost_1_62_0
  -DBoost_NO_SYSTEM_PATHS=ON)

# FIXME add to default target "all"?
ExternalProject_Add (build_full_DBGWAS_project
  DEPENDS ${DEPENDENCIES}
  SOURCE_DIR ${PROJECT_SOURCE_DIR}
  CMAKE_ARGS -DUSE_SUPERBUILD=OFF ${EXTRA_CMAKE_ARGS}
  INSTALL_COMMAND ""
  BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/DBGWAS)
