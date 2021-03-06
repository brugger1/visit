// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//  File: ParallelCoordinatesPluginInfo.C
// ************************************************************************* //

#include <ParallelCoordinatesPluginInfo.h>
#include <ParallelCoordinatesAttributes.h>

#include <visit-config.h>
VISIT_PLUGIN_VERSION(ParallelCoordinates,PLOT_EXPORT)

VISIT_PLOT_PLUGIN_ENTRY(ParallelCoordinates,General)

// ****************************************************************************
//  Method: ParallelCoordinatesGeneralPluginInfo::GetName
//
//  Purpose:
//    Return the name of the plot plugin.
//
//  Returns:    A pointer to the name of the plot plugin.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

const char *
ParallelCoordinatesGeneralPluginInfo::GetName() const
{
    return "ParallelCoordinates";
}

// ****************************************************************************
//  Method: ParallelCoordinatesGeneralPluginInfo::GetVersion
//
//  Purpose:
//    Return the version of the plot plugin.
//
//  Returns:    A pointer to the version of the plot plugin.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

const char *
ParallelCoordinatesGeneralPluginInfo::GetVersion() const
{
    return "1.0";
}

// ****************************************************************************
//  Method: ParallelCoordinatesGeneralPluginInfo::GetID
//
//  Purpose:
//    Return the id of the plot plugin.
//
//  Returns:    A pointer to the id of the plot plugin.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

const char *
ParallelCoordinatesGeneralPluginInfo::GetID() const
{
    return "ParallelCoordinates_1.0";
}
// ****************************************************************************
//  Method: ParallelCoordinatesGeneralPluginInfo::EnabledByDefault
//
//  Purpose:
//    Return true if this plugin should be enabled by default; false otherwise.
//
//  Returns:    true/false
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

bool
ParallelCoordinatesGeneralPluginInfo::EnabledByDefault() const
{
    return true;
}
