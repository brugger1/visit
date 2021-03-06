// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//                        SurfaceNormalScriptingPluginInfo.C
// ************************************************************************* //
#include <PySurfaceNormalAttributes.h>
#include <SurfaceNormalPluginInfo.h>

VISIT_OPERATOR_PLUGIN_ENTRY(SurfaceNormal,Scripting)

// ****************************************************************************
// Method: SurfaceNormalScriptingPluginInfo::InitializePlugin
//
// Purpose: 
//   Calls the initialization function for the plugin.
//
// Arguments:
//   subj    : A pointer to the plugin's state object.
//   data    : A pointer to data to be used by the observer function.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

void
SurfaceNormalScriptingPluginInfo::InitializePlugin(AttributeSubject *subj,
    void *data)
{
    PySurfaceNormalAttributes_StartUp((SurfaceNormalAttributes *)subj, data);
}

// ****************************************************************************
// Method: SurfaceNormalScriptingPluginInfo::GetMethodTable
//
// Purpose: 
//   Returns a pointer to the plugin's Python method table. These methods are
//   added to the top-level visit module's methods.
//
// Arguments:
//   nMethods : Returns the number of methods in the method table.
//
// Returns:    A pointer to the method table.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

void *
SurfaceNormalScriptingPluginInfo::GetMethodTable(int *nMethods)
{
    return PySurfaceNormalAttributes_GetMethodTable(nMethods);
}

// ****************************************************************************
// Method: SurfaceNormalScriptingPluginInfo::TypesMatch
//
// Purpose: 
//   Returns whether or not the input PyObject is SurfaceNormal plot attributes.
//
// Arguments:
//   pyobject : A PyObject cast to void*.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

bool
SurfaceNormalScriptingPluginInfo::TypesMatch(void *pyobject)
{
    return PySurfaceNormalAttributes_Check((PyObject *)pyobject);
}

// ****************************************************************************
// Method: SurfaceNormalScriptingPluginInfo::GetLogString
//
// Purpose: 
//   Gets a string representation of the current attributes.
//
// Arguments:
//   val : Whether or not to log state information.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

char *
SurfaceNormalScriptingPluginInfo::GetLogString()
{
    std::string s(PySurfaceNormalAttributes_GetLogString());
    char *v = new char[s.size() + 1];
    strcpy(v, s.c_str());
    return v;
}

// ****************************************************************************
// Method: SurfaceNormalScriptingPluginInfo::SetDefaults
//
// Purpose: 
//   Used to set the default values for a plugin's state object.
//
// Arguments:
//   atts : The new state.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

void
SurfaceNormalScriptingPluginInfo::SetDefaults(const AttributeSubject *atts)
{
    PySurfaceNormalAttributes_SetDefaults((const SurfaceNormalAttributes *)atts);
}
