// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//                             VectorPluginInfo.h                            //
// ************************************************************************* //

#ifndef VECTOR_PLUGIN_INFO_H
#define VECTOR_PLUGIN_INFO_H
#include <PlotPluginInfo.h>
#include <plot_plugin_exports.h>

class VectorAttributes;

// ****************************************************************************
//  Class: VectorPluginInfo
//
//  Purpose:
//    Five classes that provide all the information about a Vector
//    plot plugin.  The information is broken up into five classes since
//    portions of it are only relevant to particular components within
//    visit.  There is the general information which all the components
//    are interested in, the gui information which the gui is interested in,
//    the viewer information which the viewer is interested in, the
//    engine information which the engine is interested in, and finally a.
//    scripting portion that enables the Python VisIt extension to use the
//    plugin.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

class VectorGeneralPluginInfo: public virtual GeneralPlotPluginInfo
{
  public:
    virtual const char *GetName() const;
    virtual const char *GetVersion() const;
    virtual const char *GetID() const;
    virtual bool  EnabledByDefault() const;
};

class VectorCommonPluginInfo : public virtual CommonPlotPluginInfo, public virtual VectorGeneralPluginInfo
{
  public:
    virtual AttributeSubject *AllocAttributes();
    virtual void CopyAttributes(AttributeSubject *to, AttributeSubject *from);
    virtual int GetVariableTypes() const;
};

class VectorGUIPluginInfo : public virtual GUIPlotPluginInfo, public virtual VectorCommonPluginInfo
{
  public:
    virtual QString *GetMenuName() const;
    virtual QvisPostableWindowObserver *CreatePluginWindow(int type,
        AttributeSubject *attr, const QString &caption, const QString &shortName,
        QvisNotepadArea *notepad);
    virtual const char **XPMIconData() const;
};

class VectorViewerEnginePluginInfo : public virtual ViewerEnginePlotPluginInfo, public virtual VectorCommonPluginInfo
{
  public:
    virtual AttributeSubject *GetClientAtts();
    virtual AttributeSubject *GetDefaultAtts();
    virtual void SetClientAtts(AttributeSubject *atts);
    virtual void GetClientAtts(AttributeSubject *atts);

    virtual avtPlot *AllocAvtPlot();

    virtual void InitializePlotAtts(AttributeSubject *atts, const avtPlotMetaData &);
    virtual bool SupportsAnimation() const;
    virtual bool AnimationReset(AttributeSubject *atts, const avtPlotMetaData &plot);
    virtual bool AnimationStep(AttributeSubject *atts, const avtPlotMetaData  &plot);
    virtual const char *GetMenuName() const;
    static void InitializeGlobalObjects();
  private:
    static VectorAttributes *defaultAtts;
    static VectorAttributes *clientAtts;
};

class VectorViewerPluginInfo : public virtual ViewerPlotPluginInfo, public virtual VectorViewerEnginePluginInfo
{
  public:
    virtual const char **XPMIconData() const;

};

class VectorEnginePluginInfo : public virtual EnginePlotPluginInfo, public virtual VectorViewerEnginePluginInfo
{
  public:
};

class VectorScriptingPluginInfo : public virtual ScriptingPlotPluginInfo, public virtual VectorCommonPluginInfo
{
  public:
    virtual void InitializePlugin(AttributeSubject *subj, void *data);
    virtual void *GetMethodTable(int *nMethods);
    virtual bool TypesMatch(void *pyobject);
    virtual char *GetLogString();
    virtual void SetDefaults(const AttributeSubject *atts);
};

#endif
