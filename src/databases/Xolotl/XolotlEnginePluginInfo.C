// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers. See the top-level LICENSE file for dates and other
// details. No copyright assignment is required to contribute to VisIt.

#include <XolotlPluginInfo.h>

// ****************************************************************************
//  Function:  GetEngineInfo
//
//  Purpose:
//    Return a new EnginePluginInfo for the Xolotl database.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************
extern "C" DBP_EXPORT EngineDatabasePluginInfo* Xolotl_GetEngineInfo()
{
    return new XolotlEnginePluginInfo;
}

// ****************************************************************************
//  Method: XolotlEnginePluginInfo::GetWriter
//
//  Purpose:
//      Sets up a Xolotl writer.
//
//  Returns:    A Xolotl writer.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************
avtDatabaseWriter *
XolotlEnginePluginInfo::GetWriter(void)
{
    return NULL;
}
