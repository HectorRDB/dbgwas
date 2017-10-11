project(DBGWAS)

cmake_minimum_required(VERSION 2.6)

################################################################################
# The version number.
################################################################################
SET (gatb-tool_VERSION_MAJOR 0)
SET (gatb-tool_VERSION_MINOR 3)
SET (gatb-tool_VERSION_PATCH 7)

#write the version number
file(WRITE ${PROJECT_SOURCE_DIR}/tools/DBGWAS/src/version.h "\#define VERSION \"${gatb-tool_VERSION_MAJOR}.${gatb-tool_VERSION_MINOR}.${gatb-tool_VERSION_PATCH}\"")


IF (DEFINED MAJOR)
    SET (gatb-tool_VERSION_MAJOR ${MAJOR})
ENDIF()
IF (DEFINED MINOR)
    SET (gatb-tool_VERSION_MINOR ${MINOR})
ENDIF()
IF (DEFINED PATCH)
    SET (gatb-tool_VERSION_PATCH ${PATCH})
ENDIF()

set (gatb-tool-version ${gatb-tool_VERSION_MAJOR}.${gatb-tool_VERSION_MINOR}.${gatb-tool_VERSION_PATCH})

################################################################################
# Define cmake modules directory
################################################################################
SET (GATB_CORE_HOME  ${PROJECT_SOURCE_DIR}/thirdparty/gatb-core)
SET (CMAKE_MODULE_PATH ${GATB_CORE_HOME}/cmake)

################################################################################
# SUPPORTED KMER SIZES
################################################################################

# One can uncomment this line and set the wanted values
#set (KSIZE_LIST "32   64   96  128  160  192  224  256")
set (KSIZE_LIST "32   64   96  128")

################################################################################
# THIRD PARTIES
################################################################################

# We don't want to install some GATB-CORE artifacts
SET (GATB_CORE_EXCLUDE_TOOLS     1)
SET (GATB_CORE_EXCLUDE_TESTS     1)
SET (GATB_CORE_EXCLUDE_EXAMPLES  1)

# GATB CORE
include (GatbCore)

################################################################################
# TOOLS
################################################################################

ADD_SUBDIRECTORY(tools)

################################################################################
# Packaging
################################################################################
# We set the version number
SET (CPACK_PACKAGE_DESCRIPTION_SUMMARY  "gatb-tool ${PROJECT_NAME}")
SET (CPACK_PACKAGE_VENDOR               "Genscale team (INRIA)")
SET (CPACK_PACKAGE_VERSION_MAJOR        "${gatb-tool_VERSION_MAJOR}")
SET (CPACK_PACKAGE_VERSION_MINOR        "${gatb-tool_VERSION_MINOR}")
SET (CPACK_PACKAGE_VERSION_PATCH        "${gatb-tool_VERSION_PATCH}")
SET (CPACK_PACKAGE_VERSION              "${gatb-tool-version}")  

# We set the kind of archive
SET (CPACK_GENERATOR                    "TGZ")
SET (CPACK_SOURCE_GENERATOR             "TGZ")

# Packaging the source ; we ignore unwanted files 
SET (CPACK_SOURCE_IGNORE_FILES
        "^${PROJECT_SOURCE_DIR}/build/"
        "^${PROJECT_SOURCE_DIR}/builddebug/"
        "^${PROJECT_SOURCE_DIR}/buildlinux/"
        "^${PROJECT_SOURCE_DIR}/cytoscapejs/"
        "^${PROJECT_SOURCE_DIR}/TODO"
        "^${GATB_CORE_HOME}/.project"
        "^${GATB_CORE_HOME}/.gitignore"
        "^${GATB_CORE_HOME}/doc"
        "^${GATB_CORE_HOME}/DELIVERY.md"
)

# Packaging the binary ; we want to include some additional files
INSTALL (FILES   ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE                     DESTINATION .)
INSTALL (FILES   ${CMAKE_CURRENT_SOURCE_DIR}/README.md                     DESTINATION .)
#TODO: remove this, do we need it?
#INSTALL (FILES   ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/THIRDPARTIES.md  DESTINATION .)
INSTALL (FILES   ${PROJECT_SOURCE_DIR}/csjs/cytoscape_template.html  DESTINATION ./bin/)
INSTALL (PROGRAMS   ${PROJECT_SOURCE_DIR}/gemma/gemma.0.93b  DESTINATION ./bin/)
INSTALL (PROGRAMS   ${PROJECT_SOURCE_DIR}/phantomjs/phantomjs  DESTINATION ./bin/)
INSTALL (PROGRAMS   ${PROJECT_SOURCE_DIR}/blast/blastn  DESTINATION ./bin/)
INSTALL (PROGRAMS   ${PROJECT_SOURCE_DIR}/blast/blastx  DESTINATION ./bin/)
INSTALL (PROGRAMS   ${PROJECT_SOURCE_DIR}/blast/makeblastdb  DESTINATION ./bin/)
INSTALL (FILES   ${PROJECT_SOURCE_DIR}/phantomjs/render_graph.js  DESTINATION ./bin/)
INSTALL (FILES   ${PROJECT_SOURCE_DIR}/bugwas_scripts/build_bmx_data.R  DESTINATION ./bin/)
INSTALL (FILES   ${PROJECT_SOURCE_DIR}/bugwas_scripts/cdbg_all_plots.R  DESTINATION ./bin/)
INSTALL (FILES   ${PROJECT_SOURCE_DIR}/bugwas_scripts/cdbg_lin_loc.R  DESTINATION ./bin/)
INSTALL (FILES   ${PROJECT_SOURCE_DIR}/bugwas_scripts/cdbg_ridge_regression.R  DESTINATION ./bin/)
INSTALL (FILES   ${PROJECT_SOURCE_DIR}/bugwas_scripts/DBGWAS.R  DESTINATION ./bin/)
INSTALL (FILES   ${PROJECT_SOURCE_DIR}/html_templates/index_template.html  DESTINATION ./bin/)
INSTALL (FILES   ${PROJECT_SOURCE_DIR}/html_templates/index_table_template.html  DESTINATION ./bin/)
INSTALL (DIRECTORY   ${PROJECT_SOURCE_DIR}/sample_example  DESTINATION .)
INSTALL (DIRECTORY   ${PROJECT_SOURCE_DIR}/csjs/lib  DESTINATION ./bin/)
INSTALL (DIRECTORY   ${PROJECT_SOURCE_DIR}/sharedobjects  DESTINATION ./bin/)

include (CPack)
