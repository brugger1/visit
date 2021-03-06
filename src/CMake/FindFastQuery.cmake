# Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
# Project developers.  See the top-level LICENSE file for dates and other
# details.  No copyright assignment is required to contribute to VisIt.

#****************************************************************************
# Modifications:
#   Kathleen Biagas, Tues Oct 1 09:33:47 MST 2013
#   Removed VISIT_MSVC_VERSION from windows handling.
#
#****************************************************************************/

# Use the FASTQUERY_DIR hint from the config-site .cmake file 

INCLUDE(${VISIT_SOURCE_DIR}/CMake/SetUpThirdParty.cmake)

IF (WIN32)
    SET_UP_THIRD_PARTY(FASTQUERY lib include/fastquery fastquery)
ELSE (WIN32)
    IF("${VISIT_CMAKE_PLATFORM}" STREQUAL "Linux")
        # Linux requires librt to resolve "clock_gettime"
        # add this as a general dep:
        #SET(FASTQUERY_LIBDEP /usr/lib rt "${FASTQUERY_LIBDEP}")
    ENDIF("${VISIT_CMAKE_PLATFORM}" STREQUAL "Linux")
    SET_UP_THIRD_PARTY(FASTQUERY lib include/fastquery fastquery)
ENDIF (WIN32)

