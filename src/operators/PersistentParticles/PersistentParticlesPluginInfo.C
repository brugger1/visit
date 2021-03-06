// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//  File: PersistentParticlesPluginInfo.C
// ************************************************************************* //

#include <PersistentParticlesPluginInfo.h>
#include <PersistentParticlesAttributes.h>

#include <visit-config.h>
VISIT_PLUGIN_VERSION(PersistentParticles,OP_EXPORT)

VISIT_OPERATOR_PLUGIN_ENTRY(PersistentParticles,General)

// ****************************************************************************
//  Method: PersistentParticlesGeneralPluginInfo::GetName
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
PersistentParticlesGeneralPluginInfo::GetName() const
{
    return "PersistentParticles";
}

// ****************************************************************************
//  Method: PersistentParticlesGeneralPluginInfo::GetVersion
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
PersistentParticlesGeneralPluginInfo::GetVersion() const
{
    return "2.0";
}

// ****************************************************************************
//  Method: PersistentParticlesGeneralPluginInfo::GetID
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
PersistentParticlesGeneralPluginInfo::GetID() const
{
    return "PersistentParticles_2.0";
}
// ****************************************************************************
//  Method: PersistentParticlesGeneralPluginInfo::EnabledByDefault
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
PersistentParticlesGeneralPluginInfo::EnabledByDefault() const
{
    return true;
}

// ****************************************************************************
//  Method: PersistentParticlesGeneralPluginInfo::GetCategoryName
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
PersistentParticlesGeneralPluginInfo::GetCategoryName() const
{
    return "Analysis";
}
