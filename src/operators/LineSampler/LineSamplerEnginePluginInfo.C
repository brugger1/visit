// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//  File: LineSamplerEnginePluginInfo.C
// ************************************************************************* //

#include <LineSamplerPluginInfo.h>
#include <avtLineSamplerFilter.h>

VISIT_OPERATOR_PLUGIN_ENTRY_EV(LineSampler,Engine)

// ****************************************************************************
//  Method: LineSamplerEnginePluginInfo::AllocAvtPluginFilter
//
//  Purpose:
//    Return a pointer to a newly allocated avtPluginFilter.
//
//  Returns:    A pointer to the newly allocated avtPluginFilter.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

avtPluginFilter *
LineSamplerEnginePluginInfo::AllocAvtPluginFilter()
{
    return new avtLineSamplerFilter;
}
