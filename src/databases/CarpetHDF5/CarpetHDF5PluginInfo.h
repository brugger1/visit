// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ****************************************************************************
//                               CarpetHDF5PluginInfo.h
// ****************************************************************************

#ifndef CARPETHDF5_PLUGIN_INFO_H
#define CARPETHDF5_PLUGIN_INFO_H
#include <DatabasePluginInfo.h>
#include <database_plugin_exports.h>

class avtDatabase;
class avtDatabaseWriter;

// ****************************************************************************
//  Class: CarpetHDF5DatabasePluginInfo
//
//  Purpose:
//    Classes that provide all the information about the CarpetHDF5 plugin.
//    Portions are separated into pieces relevant to the appropriate
//    components of VisIt.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
//  Modifications:
//
// ****************************************************************************

class CarpetHDF5GeneralPluginInfo : public virtual GeneralDatabasePluginInfo
{
  public:
    virtual const char *GetName() const;
    virtual const char *GetVersion() const;
    virtual const char *GetID() const;
    virtual bool  EnabledByDefault() const;
    virtual bool  HasWriter() const;
    virtual std::vector<std::string> GetDefaultFilePatterns() const;
    virtual bool  AreDefaultFilePatternsStrict() const;
    virtual bool  OpensWholeDirectory() const;
};

class CarpetHDF5CommonPluginInfo : public virtual CommonDatabasePluginInfo, public virtual CarpetHDF5GeneralPluginInfo
{
  public:
    virtual DatabaseType              GetDatabaseType();
    virtual avtDatabase              *SetupDatabase(const char * const *list,
                                                    int nList, int nBlock);
};

class CarpetHDF5MDServerPluginInfo : public virtual MDServerDatabasePluginInfo, public virtual CarpetHDF5CommonPluginInfo
{
  public:
    // this makes compilers happy... remove if we ever have functions here
    virtual void dummy();
};

class CarpetHDF5EnginePluginInfo : public virtual EngineDatabasePluginInfo, public virtual CarpetHDF5CommonPluginInfo
{
  public:
    virtual avtDatabaseWriter        *GetWriter(void);
};

#endif
