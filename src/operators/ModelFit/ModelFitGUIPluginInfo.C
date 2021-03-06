// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//  File: ModelFitGUIPluginInfo.C
// ************************************************************************* //

#include <ModelFitPluginInfo.h>
#include <ModelFitAtts.h>
#include <QApplication>
#include <QvisModelFitWindow.h>

VISIT_OPERATOR_PLUGIN_ENTRY(ModelFit,GUI)

// ****************************************************************************
//  Method: ModelFitGUIPluginInfo::GetMenuName
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
ModelFitGUIPluginInfo::GetMenuName() const
{
    return new QString(qApp->translate("OperatorNames", "ModelFit"));
}


// ****************************************************************************
//  Method: ModelFitGUIPluginInfo::CreatePluginWindow
//
//  Purpose:
//    Return a pointer to an operator's attribute window.
//
//  Arguments:
//    type      The type of the operator.
//    attr      The attribute subject for the operator.
//    notepad   The notepad to use for posting the window.
//
//  Returns:    A pointer to the operator's attribute window.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

QvisPostableWindowObserver *
ModelFitGUIPluginInfo::CreatePluginWindow(int type, AttributeSubject *attr,
    const QString &caption, const QString &shortName, QvisNotepadArea *notepad)
{
    return new QvisModelFitWindow(type, (ModelFitAtts *)attr,
        caption, shortName, notepad);
}

