// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//  File: PICS_TesterPluginInfo.C
// ************************************************************************* //

#include <PICS_TesterPluginInfo.h>

#include <visit-config.h>
VISIT_PLUGIN_VERSION(PICS_Tester,DBP_EXPORT)

VISIT_DATABASE_PLUGIN_ENTRY(PICS_Tester,General)

// ****************************************************************************
//  Method: PICS_TesterGeneralPluginInfo::GetName
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
PICS_TesterGeneralPluginInfo::GetName() const
{
    return "PICS_Tester";
}

// ****************************************************************************
//  Method: PICS_TesterGeneralPluginInfo::GetVersion
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
PICS_TesterGeneralPluginInfo::GetVersion() const
{
    return "1.0";
}

// ****************************************************************************
//  Method: PICS_TesterGeneralPluginInfo::GetID
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
PICS_TesterGeneralPluginInfo::GetID() const
{
    return "PICS_Tester_1.0";
}
// ****************************************************************************
//  Method: PICS_TesterGeneralPluginInfo::EnabledByDefault
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
PICS_TesterGeneralPluginInfo::EnabledByDefault() const
{
    return true;
}
// ****************************************************************************
//  Method: PICS_TesterGeneralPluginInfo::HasWriter
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
PICS_TesterGeneralPluginInfo::HasWriter() const
{
    return false;
}
// ****************************************************************************
//  Method:  PICS_TesterGeneralPluginInfo::GetDefaultFilePatterns
//
//  Purpose:
//    Returns the default patterns for a PICS_Tester database.
//
//  Programmer:  generated by xml2info
//  Creation:    omitted
//
// ****************************************************************************
std::vector<std::string>
PICS_TesterGeneralPluginInfo::GetDefaultFilePatterns() const
{
    std::vector<std::string> defaultPatterns;
    defaultPatterns.push_back("*.pics");

    return defaultPatterns;
}

// ****************************************************************************
//  Method:  PICS_TesterGeneralPluginInfo::AreDefaultFilePatternsStrict
//
//  Purpose:
//    Returns if the file patterns for a PICS_Tester database are
//    intended to be interpreted strictly by default.
//
//  Programmer:  generated by xml2info
//  Creation:    omitted
//
// ****************************************************************************
bool
PICS_TesterGeneralPluginInfo::AreDefaultFilePatternsStrict() const
{
    return false;
}

// ****************************************************************************
//  Method:  PICS_TesterGeneralPluginInfo::OpensWholeDirectory
//
//  Purpose:
//    Returns if the PICS_Tester plugin opens a whole directory name
//    instead of a single file.
//
//  Programmer:  generated by xml2info
//  Creation:    omitted
//
// ****************************************************************************
bool
PICS_TesterGeneralPluginInfo::OpensWholeDirectory() const
{
    return false;
}
