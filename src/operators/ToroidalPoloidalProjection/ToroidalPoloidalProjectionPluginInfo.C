// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//  File: ToroidalPoloidalProjectionPluginInfo.C
// ************************************************************************* //

#include <ToroidalPoloidalProjectionPluginInfo.h>
#include <ToroidalPoloidalProjection.h>

#include <visit-config.h>
VISIT_PLUGIN_VERSION(ToroidalPoloidalProjection,OP_EXPORT)

VISIT_OPERATOR_PLUGIN_ENTRY(ToroidalPoloidalProjection,General)

// ****************************************************************************
//  Method: ToroidalPoloidalProjectionGeneralPluginInfo::GetName
//
//  Purpose:
//    Return the name of the operator plugin.
//
//  Returns:    A pointer to the name of the operator plugin.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

const char *
ToroidalPoloidalProjectionGeneralPluginInfo::GetName() const
{
    return "ToroidalPoloidalProjection";
}

// ****************************************************************************
//  Method: ToroidalPoloidalProjectionGeneralPluginInfo::GetVersion
//
//  Purpose:
//    Return the version of the operator plugin.
//
//  Returns:    A pointer to the version of the operator plugin.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

const char *
ToroidalPoloidalProjectionGeneralPluginInfo::GetVersion() const
{
    return "1.0";
}

// ****************************************************************************
//  Method: ToroidalPoloidalProjectionGeneralPluginInfo::GetID
//
//  Purpose:
//    Return the id of the operator plugin.
//
//  Returns:    A pointer to the id of the operator plugin.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

const char *
ToroidalPoloidalProjectionGeneralPluginInfo::GetID() const
{
    return "ToroidalPoloidalProjection_1.0";
}
// ****************************************************************************
//  Method: ToroidalPoloidalProjectionGeneralPluginInfo::EnabledByDefault
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
ToroidalPoloidalProjectionGeneralPluginInfo::EnabledByDefault() const
{
    return false;
}

// ****************************************************************************
//  Method: ToroidalPoloidalProjectionGeneralPluginInfo::GetCategoryName
//
//  Purpose:
//    Return the category name to which the operator belongs.
//
//  Returns:    Return the category name to which the operator belongs.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

const char *
ToroidalPoloidalProjectionGeneralPluginInfo::GetCategoryName() const
{
    return "Transforms";
}
