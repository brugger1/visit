// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#include <paraDISPluginInfo.h>
#include <avtparaDISFileFormat.h>
#include <avtSTSDFileFormatInterface.h>
#include <avtGenericDatabase.h>
#include <avtparaDISOptions.h>

// ****************************************************************************
//  Method:  paraDISCommonPluginInfo::GetDatabaseType
//
//  Purpose:
//    Returns the type of a paraDIS database.
//
//  Programmer:  generated by xml2info
//  Creation:    omitted
//
// ****************************************************************************
DatabaseType
paraDISCommonPluginInfo::GetDatabaseType()
{
    return DB_TYPE_STSD;
}

// ****************************************************************************
//  Method: paraDISCommonPluginInfo::SetupDatabase
//
//  Purpose:
//      Sets up a paraDIS database.
//
//  Arguments:
//      list    A list of file names.
//      nList   The number of timesteps in list.
//      nBlocks The number of blocks in the list.
//
//  Returns:    A paraDIS database from list.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************
avtDatabase *
paraDISCommonPluginInfo::SetupDatabase(const char *const *list,
                                   int nList, int nBlock)
{
    int nTimestep = nList / nBlock;
    avtSTSDFileFormat ***ffl = new avtSTSDFileFormat**[nTimestep];
    for (int i = 0 ; i < nTimestep ; i++)
    {
        ffl[i] = new avtSTSDFileFormat*[nBlock];
        for (int j = 0 ; j < nBlock ; j++)
        {
            ffl[i][j] = new avtparaDISFileFormat(list[i*nBlock + j], readOptions);
        }
    }
    avtSTSDFileFormatInterface *inter 
           = new avtSTSDFileFormatInterface(ffl, nTimestep, nBlock);
    return new avtGenericDatabase(inter);
}

// ****************************************************************************
//  Method: paraDISCommonPluginInfo::GetReadOptions
//
//  Purpose:
//      Gets the read options.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

DBOptionsAttributes *
paraDISCommonPluginInfo::GetReadOptions() const
{
    return GetparaDISReadOptions();
}
// ****************************************************************************
//  Method: paraDISCommonPluginInfo::GetWriteOptions
//
//  Purpose:
//      Gets the write options.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

DBOptionsAttributes *
paraDISCommonPluginInfo::GetWriteOptions() const
{
    return GetparaDISWriteOptions();
}
