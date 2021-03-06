// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//  File: BoundaryGUIPluginInfo.C
// ************************************************************************* //

#include <BoundaryPluginInfo.h>
#include <BoundaryAttributes.h>
#include <QApplication>
#include <QvisBoundaryPlotWindow.h>

VISIT_PLOT_PLUGIN_ENTRY(Boundary,GUI)

// ****************************************************************************
//  Method: BoundaryGUIPluginInfo::GetMenuName
//
//  Purpose:
//    Return a pointer to the name to use in the GUI menu.
//
//  Returns:    A pointer to the name to use in the GUI menu.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

QString *
BoundaryGUIPluginInfo::GetMenuName() const
{
    return new QString(qApp->translate("PlotNames", "Boundary"));
}


// ****************************************************************************
//  Method: BoundaryGUIPluginInfo::CreatePluginWindow
//
//  Purpose:
//    Return a pointer to an plot's attribute window.
//
//  Arguments:
//    type      The type of the plot.
//    attr      The attribute subject for the plot.
//    notepad   The notepad to use for posting the window.
//
//  Returns:    A pointer to the plot's attribute window.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

QvisPostableWindowObserver *
BoundaryGUIPluginInfo::CreatePluginWindow(int type, AttributeSubject *attr,
    const QString &caption, const QString &shortName, QvisNotepadArea *notepad)
{
    return new QvisBoundaryPlotWindow(type, (BoundaryAttributes *)attr,
        caption, shortName, notepad);
}

// ****************************************************************************
//  Method: BoundaryGUIPluginInfo::XPMIconData
//
//  Purpose:
//    Return a pointer to the icon data.
//
//  Returns:    A pointer to the icon data.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

#include <Boundary.xpm>
const char **
BoundaryGUIPluginInfo::XPMIconData() const
{
    return Boundary_xpm;
}

