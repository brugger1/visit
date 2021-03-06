// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#include <DenovoPluginInfo.h>
#include <avtDenovoFileFormat.h>
#include <avtSTMDFileFormatInterface.h>
#include <avtGenericDatabase.h>

// ****************************************************************************
//  Method:  DenovoCommonPluginInfo::GetDatabaseType
//
//  Purpose:
//    Returns the type of a Denovo database.
//
//  Programmer:  generated by xml2info
//  Creation:    omitted
//
// ****************************************************************************
DatabaseType
DenovoCommonPluginInfo::GetDatabaseType()
{
    return DB_TYPE_STMD;
}

// ****************************************************************************
//  Method: DenovoCommonPluginInfo::SetupDatabase
//
//  Purpose:
//      Sets up a Denovo database.
//
//  Arguments:
//      list    A list of file names.
//      nList   The number of timesteps in list.
//      nBlocks The number of blocks in the list.
//
//  Returns:    A Denovo database from list.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************
avtDatabase *
DenovoCommonPluginInfo::SetupDatabase(const char *const *list,
                                   int nList, int nBlock)
{
    avtSTMDFileFormat **ffl = new avtSTMDFileFormat*[nList];
    for (int i = 0 ; i < nList ; i++)
    {
        ffl[i] = new avtDenovoFileFormat(list[i]);
    }
    avtSTMDFileFormatInterface *inter 
           = new avtSTMDFileFormatInterface(ffl, nList);
    return new avtGenericDatabase(inter);
}
// ****************************************************************************
//  Method: DenovoCommonPluginInfo::GetLicense
//
//  Purpose:
//      Gets the write options.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

std::string
DenovoCommonPluginInfo::GetLicense() const
{
    return std::string();
}
