// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#include <MFIXCDFPluginInfo.h>
#include <avtMFIXCDFFileFormat.h>
#include <avtSTMDFileFormatInterface.h>
#include <avtGenericDatabase.h>
#include <avtMFIXCDFOptions.h>

// ****************************************************************************
//  Method:  MFIXCDFCommonPluginInfo::GetDatabaseType
//
//  Purpose:
//    Returns the type of a MFIXCDF database.
//
//  Programmer:  generated by xml2info
//  Creation:    omitted
//
// ****************************************************************************
DatabaseType
MFIXCDFCommonPluginInfo::GetDatabaseType()
{
    return DB_TYPE_STMD;
}

// ****************************************************************************
//  Method: MFIXCDFCommonPluginInfo::SetupDatabase
//
//  Purpose:
//      Sets up a MFIXCDF database.
//
//  Arguments:
//      list    A list of file names.
//      nList   The number of timesteps in list.
//      nBlocks The number of blocks in the list.
//
//  Returns:    A MFIXCDF database from list.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************
avtDatabase *
MFIXCDFCommonPluginInfo::SetupDatabase(const char *const *list,
                                   int nList, int nBlock)
{
    avtSTMDFileFormat **ffl = new avtSTMDFileFormat*[nList];
    for (int i = 0 ; i < nList ; i++)
    {
        ffl[i] = new avtMFIXCDFFileFormat(list[i], readOptions);
    }
    avtSTMDFileFormatInterface *inter 
           = new avtSTMDFileFormatInterface(ffl, nList);
    return new avtGenericDatabase(inter);
}

// ****************************************************************************
//  Method: MFIXCDFCommonPluginInfo::GetReadOptions
//
//  Purpose:
//      Gets the read options.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

DBOptionsAttributes *
MFIXCDFCommonPluginInfo::GetReadOptions() const
{
    return GetMFIXCDFReadOptions();
}
// ****************************************************************************
//  Method: MFIXCDFCommonPluginInfo::GetWriteOptions
//
//  Purpose:
//      Gets the write options.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

DBOptionsAttributes *
MFIXCDFCommonPluginInfo::GetWriteOptions() const
{
    return GetMFIXCDFWriteOptions();
}
