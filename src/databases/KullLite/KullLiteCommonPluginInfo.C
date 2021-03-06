// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#include <KullLitePluginInfo.h>
#include <avtKullLiteFileFormat.h>
#include <avtSTMDFileFormatInterface.h>
#include <avtGenericDatabase.h>

// ****************************************************************************
//  Method:  KullLiteCommonPluginInfo::GetDatabaseType
//
//  Purpose:
//    Returns the type of a KullLite database.
//
//  Programmer:  generated by xml2info
//  Creation:    omitted
//
// ****************************************************************************
DatabaseType
KullLiteCommonPluginInfo::GetDatabaseType()
{
    return DB_TYPE_STMD;
}

// ****************************************************************************
//  Method: KullLiteCommonPluginInfo::SetupDatabase
//
//  Purpose:
//      Sets up a KullLite database.
//
//  Arguments:
//      list    A list of file names.
//      nList   The number of timesteps in list.
//      nBlocks The number of blocks in the list.
//
//  Returns:    A KullLite database from list.
//
//  Programmer: childs -- generated by xml2info
//  Creation:   Mon Mar 22 09:17:37 PDT 2004
//
//  Modifications:
//    Brad Whitlock, Thu Oct 9 14:25:58 PST 2003
//    Fixed a memory leak that can happen if the file format decides not to
//    actually open the file.
//
// ****************************************************************************

avtDatabase *
KullLiteCommonPluginInfo::SetupDatabase(const char *const *list,
                                   int nList, int nBlock)
{
    avtSTMDFileFormat **ffl = new avtSTMDFileFormat*[nList];
    for (int i = 0 ; i < nList ; i++)
        ffl[i] = 0;

    TRY
    {
        for (int i = 0 ; i < nList ; i++)
            ffl[i] = new avtKullLiteFileFormat(list[i]);
    }
    CATCH(VisItException)
    {
        for (int i = 0 ; i < nList ; i++)
            delete ffl[i];
        delete [] ffl;

        RETHROW;
    }
    ENDTRY

    avtSTMDFileFormatInterface *inter 
           = new avtSTMDFileFormatInterface(ffl, nList);
    return new avtGenericDatabase(inter);
}

