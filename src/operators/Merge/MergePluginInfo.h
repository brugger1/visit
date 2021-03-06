// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//  File: MergePluginInfo.h
// ************************************************************************* //

#ifndef MERGE_PLUGIN_INFO_H
#define MERGE_PLUGIN_INFO_H
#include <OperatorPluginInfo.h>
#include <operator_plugin_exports.h>

class MergeOperatorAttributes;

// ****************************************************************************
//  Class: MergePluginInfo
//
//  Purpose:
//    Five classes that provide all the information about an Merge operator
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
//  Modifications:
//
// ****************************************************************************

class MergeGeneralPluginInfo : public virtual GeneralOperatorPluginInfo
{
  public:
    virtual const char *GetName() const;
    virtual const char *GetVersion() const;
    virtual const char *GetID() const;
    virtual bool  EnabledByDefault() const;
    virtual const char *GetCategoryName() const;
};

class MergeCommonPluginInfo : public virtual CommonOperatorPluginInfo, public virtual MergeGeneralPluginInfo
{
  public:
    virtual AttributeSubject *AllocAttributes();
    virtual void CopyAttributes(AttributeSubject *to, AttributeSubject *from);
};

class MergeGUIPluginInfo : public virtual GUIOperatorPluginInfo, public virtual MergeCommonPluginInfo
{
  public:
    virtual QString *GetMenuName() const;
    virtual QvisPostableWindowObserver *CreatePluginWindow(int type,
        AttributeSubject *attr, const QString &caption, const QString &shortName,
        QvisNotepadArea *notepad);
};

class MergeViewerEnginePluginInfo : public virtual ViewerEngineOperatorPluginInfo, public virtual MergeCommonPluginInfo
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
    static MergeOperatorAttributes *defaultAtts;
    static MergeOperatorAttributes *clientAtts;
};

class MergeViewerPluginInfo : public virtual ViewerOperatorPluginInfo, public virtual MergeViewerEnginePluginInfo
{
  public:
};

class MergeEnginePluginInfo : public virtual EngineOperatorPluginInfo, public virtual MergeViewerEnginePluginInfo
{
  public:
    virtual avtPluginFilter *AllocAvtPluginFilter();
};

class MergeScriptingPluginInfo : public virtual ScriptingOperatorPluginInfo, public virtual MergeCommonPluginInfo
{
  public:
    virtual void InitializePlugin(AttributeSubject *subj, void *data);
    virtual void *GetMethodTable(int *nMethods);
    virtual bool TypesMatch(void *pyobject);
    virtual char *GetLogString();
    virtual void SetDefaults(const AttributeSubject *atts);
};

#endif
