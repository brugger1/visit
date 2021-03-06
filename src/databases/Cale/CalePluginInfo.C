// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//  File: CalePluginInfo.C
// ************************************************************************* //

#include <CalePluginInfo.h>

#include <visit-config.h>
VISIT_PLUGIN_VERSION(Cale,DBP_EXPORT)

VISIT_DATABASE_PLUGIN_ENTRY(Cale,General)

// ****************************************************************************
//  Method: CaleGeneralPluginInfo::GetName
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
CaleGeneralPluginInfo::GetName() const
{
    return "Cale";
}

// ****************************************************************************
//  Method: CaleGeneralPluginInfo::GetVersion
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
CaleGeneralPluginInfo::GetVersion() const
{
    return "1.1";
}

// ****************************************************************************
//  Method: CaleGeneralPluginInfo::GetID
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
CaleGeneralPluginInfo::GetID() const
{
    return "Cale_1.1";
}
// ****************************************************************************
//  Method: CaleGeneralPluginInfo::EnabledByDefault
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
CaleGeneralPluginInfo::EnabledByDefault() const
{
    return true;
}
// ****************************************************************************
//  Method: CaleGeneralPluginInfo::HasWriter
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
CaleGeneralPluginInfo::HasWriter() const
{
    return false;
}
// ****************************************************************************
//  Method:  CaleGeneralPluginInfo::GetDefaultFilePatterns
//
//  Purpose:
//    Returns the default patterns for a Cale database.
//
//  Programmer:  generated by xml2info
//  Creation:    omitted
//
// ****************************************************************************
std::vector<std::string>
CaleGeneralPluginInfo::GetDefaultFilePatterns() const
{
    std::vector<std::string> defaultPatterns;
    defaultPatterns.push_back("*.pdb");
    defaultPatterns.push_back("*.cale");

    return defaultPatterns;
}

// ****************************************************************************
//  Method:  CaleGeneralPluginInfo::AreDefaultFilePatternsStrict
//
//  Purpose:
//    Returns if the file patterns for a Cale database are
//    intended to be interpreted strictly by default.
//
//  Programmer:  generated by xml2info
//  Creation:    omitted
//
// ****************************************************************************
bool
CaleGeneralPluginInfo::AreDefaultFilePatternsStrict() const
{
    return false;
}

// ****************************************************************************
//  Method:  CaleGeneralPluginInfo::OpensWholeDirectory
//
//  Purpose:
//    Returns if the Cale plugin opens a whole directory name
//    instead of a single file.
//
//  Programmer:  generated by xml2info
//  Creation:    omitted
//
// ****************************************************************************
bool
CaleGeneralPluginInfo::OpensWholeDirectory() const
{
    return false;
}
