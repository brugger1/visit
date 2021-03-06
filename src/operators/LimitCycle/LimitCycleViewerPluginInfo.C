// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//  File: LimitCycleViewerPluginInfo.C
// ************************************************************************* //

#include <LimitCyclePluginInfo.h>
#include <LimitCycleAttributes.h>

VISIT_OPERATOR_PLUGIN_ENTRY_EV(LimitCycle,Viewer)


// ****************************************************************************
//  Method: LimitCycleViewerPluginInfo::XPMIconData
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

#include <LimitCycle.xpm>
const char **
LimitCycleViewerPluginInfo::XPMIconData() const
{
    return LimitCycle_xpm;
}

