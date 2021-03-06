// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//  File: unvPluginInfo.C
// ************************************************************************* //

#include <unvPluginInfo.h>

#include <visit-config.h>
VISIT_PLUGIN_VERSION(unv,DBP_EXPORT)

VISIT_DATABASE_PLUGIN_ENTRY(unv,General)

// ****************************************************************************
//  Method: unvGeneralPluginInfo::GetName
//
//  Purpose:
//    Return the name of the database plugin.
//
//  Returns:    A pointer to the name of the database plugin.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

const char *
unvGeneralPluginInfo::GetName() const
{
    return "unv";
}

// ****************************************************************************
//  Method: unvGeneralPluginInfo::GetVersion
//
//  Purpose:
//    Return the version of the database plugin.
//
//  Returns:    A pointer to the version of the database plugin.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

const char *
unvGeneralPluginInfo::GetVersion() const
{
    return "";
}

// ****************************************************************************
//  Method: unvGeneralPluginInfo::GetID
//
//  Purpose:
//    Return the id of the database plugin.
//
//  Returns:    A pointer to the id of the database plugin.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

const char *
unvGeneralPluginInfo::GetID() const
{
    return "unv_";
}
// ****************************************************************************
//  Method: unvGeneralPluginInfo::EnabledByDefault
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
unvGeneralPluginInfo::EnabledByDefault() const
{
    return true;
}
// ****************************************************************************
//  Method: unvGeneralPluginInfo::HasWriter
//
//  Purpose:
//    Return true if this plugin has a database writer.
//
//  Returns:    true/false
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

bool
unvGeneralPluginInfo::HasWriter() const
{
    return false;
}
// ****************************************************************************
//  Method:  unvGeneralPluginInfo::GetDefaultFilePatterns
//
//  Purpose:
//    Returns the default patterns for a unv database.
//
//  Programmer:  generated by xml2info
//  Creation:    omitted
//
// ****************************************************************************
std::vector<std::string>
unvGeneralPluginInfo::GetDefaultFilePatterns() const
{
    std::vector<std::string> defaultPatterns;
    defaultPatterns.push_back("*.unv");
    defaultPatterns.push_back("*.unv.gz");

    return defaultPatterns;
}

// ****************************************************************************
//  Method:  unvGeneralPluginInfo::AreDefaultFilePatternsStrict
//
//  Purpose:
//    Returns if the file patterns for a unv database are
//    intended to be interpreted strictly by default.
//
//  Programmer:  generated by xml2info
//  Creation:    omitted
//
// ****************************************************************************
bool
unvGeneralPluginInfo::AreDefaultFilePatternsStrict() const
{
    return false;
}

// ****************************************************************************
//  Method:  unvGeneralPluginInfo::OpensWholeDirectory
//
//  Purpose:
//    Returns if the unv plugin opens a whole directory name
//    instead of a single file.
//
//  Programmer:  generated by xml2info
//  Creation:    omitted
//
// ****************************************************************************
bool
unvGeneralPluginInfo::OpensWholeDirectory() const
{
    return false;
}
