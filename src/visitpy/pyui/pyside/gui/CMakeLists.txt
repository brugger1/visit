# Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
# Project developers.  See the top-level LICENSE file for dates and other
# details.  No copyright assignment is required to contribute to VisIt.

#*****************************************************************************
# Modifications:
#    Gunther H. Weber, Thu Nov 21 18:24:53 PST 2013
#    Use instead PYSIDE_INCLUDE_DIR (set by findpackage PySide) instead of
#    VISIT_PYSIDE_DIR to support system PySide installs.
#
#    Kathleen Biagas, Tue Jan 24 11:06:12 PST 2017
#    Add support for Qt5.
#
#*****************************************************************************

set(pyside_gui_include_paths
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${SHIBOKEN_INCLUDE_DIR}
    ${SHIBOKEN_PYTHON_INCLUDE_DIR}
    ${PYSIDE_INCLUDE_DIR}
    ${pyside_qt_includes}
    ${PYSIDE_INCLUDE_DIR}
    ${PYSIDE_INCLUDE_DIR}/QtCore
    ${PYSIDE_INCLUDE_DIR}/QtGui
    ${PYSIDE_INCLUDE_DIR}/QtOpenGL
    ${PYSIDE_INCLUDE_DIR}/QtWidgets
    ${VISIT_SOURCE_DIR}/gui/
    ${VISIT_SOURCE_DIR}/gui/main
    ${VISIT_SOURCE_DIR}/viewer/subjectproxy
    ${VISIT_SOURCE_DIR}/viewer/core
    ${VISIT_SOURCE_DIR}/viewer/main 
    ${VISIT_SOURCE_DIR}/viewer/main/ui
    ${VISIT_SOURCE_DIR}/viewer/proxy
    ${VISIT_SOURCE_DIR}/viewer/rpc
    ${VISIT_SOURCE_DIR}/vtkqt
    ${VISIT_SOURCE_DIR}/common/utility
    ${VISIT_SOURCE_DIR}/common/state
    ${VISIT_SOURCE_DIR}/common/misc
    ${VISIT_SOURCE_DIR}/common/comm
    ${VISIT_SOURCE_DIR}/common/plugin
    ${VISIT_SOURCE_DIR}/common/Exceptions/Pipeline
    ${VISIT_SOURCE_DIR}/avt/DBAtts/MetaData/
    ${VISIT_SOURCE_DIR}/avt/DBAtts/SIL
    ${VISIT_SOURCE_DIR}/avt/Math
    ${VISIT_SOURCE_DIR}/avt/VisWindow/VisWindow
    ${VISIT_SOURCE_DIR}/winutil
    ${VISIT_SOURCE_DIR}/visitpy/pyui/common
)

set(pyside_gui_link_libs
    #gui
    guiwrapper
    #viewersubjectproxy
    #${SHIBOKEN_PYTHON_LIBRARIES}
    #${SHIBOKEN_LIBRARY}
    #${PYSIDE_LIBRARY}
    #${QT_QTCORE_LIBRARY}
    ${QWT_LIBRARY}
    )

set(pyside_gui_sources ${CMAKE_CURRENT_SOURCE_DIR}/pysidegui.C) 
set(pyside_gui_gen_sources 
${CMAKE_CURRENT_BINARY_DIR}/pyside_gui/pysidegui_wrapper.cpp
${CMAKE_CURRENT_BINARY_DIR}/pyside_gui/guiwrapper_wrapper.cpp
${CMAKE_CURRENT_BINARY_DIR}/pyside_gui/pyside_gui_module_wrapper.cpp)

LINK_DIRECTORIES(${QWT_LIBRARY_DIR})

#
# For dev build, we may need to link to mesa on Linux.
#
IF(OSMESA_FOUND)
    set(pyside_gui_link_libs ${pyside_gui_link_libs} ${OSMESA_LIBRARIES})
ENDIF()

PYSIDE_ADD_MODULE(pyside_gui
                  site-packages/visit
                  pyside_gui_sources
                  pyside_gui_gen_sources
                  pyside_gui_include_paths
                  pyside_gui_link_libs
                  global.h typesystem.xml)

