// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//  File: ConnCompReducePluginInfo.h
// ************************************************************************* //

#ifndef CONNCOMPREDUCE_PLUGIN_INFO_H
#define CONNCOMPREDUCE_PLUGIN_INFO_H
#include <OperatorPluginInfo.h>
#include <operator_plugin_exports.h>

class ConnCompReduceAttributes;

// ****************************************************************************
//  Class: ConnCompReducePluginInfo
//
//  Purpose:
//    Five classes that provide all the information about an ConnCompReduce operator
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
//  Modifications:
//
// ****************************************************************************

class ConnCompReduceGeneralPluginInfo : public virtual GeneralOperatorPluginInfo
{
  public:
    virtual const char *GetName() const;
    virtual const char *GetVersion() const;
    virtual const char *GetID() const;
    virtual bool  EnabledByDefault() const;
    virtual const char *GetCategoryName() const;
};

class ConnCompReduceCommonPluginInfo : public virtual CommonOperatorPluginInfo, public virtual ConnCompReduceGeneralPluginInfo
{
  public:
    virtual AttributeSubject *AllocAttributes();
    virtual void CopyAttributes(AttributeSubject *to, AttributeSubject *from);
};

class ConnCompReduceGUIPluginInfo : public virtual GUIOperatorPluginInfo, public virtual ConnCompReduceCommonPluginInfo
{
  public:
    virtual QString *GetMenuName() const;
    virtual QvisPostableWindowObserver *CreatePluginWindow(int type,
        AttributeSubject *attr, const QString &caption, const QString &shortName,
        QvisNotepadArea *notepad);
};

class ConnCompReduceViewerEnginePluginInfo : public virtual ViewerEngineOperatorPluginInfo, public virtual ConnCompReduceCommonPluginInfo
{
  public:
    virtual AttributeSubject *GetClientAtts();
    virtual AttributeSubject *GetDefaultAtts();
    virtual void SetClientAtts(AttributeSubject *atts);
    virtual void GetClientAtts(AttributeSubject *atts);

    virtual void InitializeOperatorAtts(AttributeSubject *atts,
                                        const avtPlotMetaData &plot,
                                        const bool fromDefault);
    virtual void UpdateOperatorAtts(AttributeSubject *atts,
                                    const avtPlotMetaData &plot);
    virtual const char *GetMenuName() const;

    static void InitializeGlobalObjects();
  private:
    static ConnCompReduceAttributes *defaultAtts;
    static ConnCompReduceAttributes *clientAtts;
};

class ConnCompReduceViewerPluginInfo : public virtual ViewerOperatorPluginInfo, public virtual ConnCompReduceViewerEnginePluginInfo
{
  public:
};

class ConnCompReduceEnginePluginInfo : public virtual EngineOperatorPluginInfo, public virtual ConnCompReduceViewerEnginePluginInfo
{
  public:
    virtual avtPluginFilter *AllocAvtPluginFilter();
};

class ConnCompReduceScriptingPluginInfo : public virtual ScriptingOperatorPluginInfo, public virtual ConnCompReduceCommonPluginInfo
{
  public:
    virtual void InitializePlugin(AttributeSubject *subj, void *data);
    virtual void *GetMethodTable(int *nMethods);
    virtual bool TypesMatch(void *pyobject);
    virtual char *GetLogString();
    virtual void SetDefaults(const AttributeSubject *atts);
};

#endif
