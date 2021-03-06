// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//  File: TensorViewerEnginePluginInfo.C
// ************************************************************************* //

#include <TensorPluginInfo.h>
#include <avtTensorPlot.h>
#include <TensorAttributes.h>

//
// Storage for static data elements.
//
TensorAttributes *TensorViewerEnginePluginInfo::clientAtts = NULL;
TensorAttributes *TensorViewerEnginePluginInfo::defaultAtts = NULL;

// ****************************************************************************
//  Method:  TensorViewerEnginePluginInfo::InitializeGlobalObjects
//
//  Purpose:
//    Initialize the plot atts.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************
void
TensorViewerEnginePluginInfo::InitializeGlobalObjects()
{
    TensorViewerEnginePluginInfo::clientAtts  = new TensorAttributes;
    TensorViewerEnginePluginInfo::defaultAtts = new TensorAttributes;
}

// ****************************************************************************
//  Method: TensorViewerEnginePluginInfo::GetClientAtts
//
//  Purpose:
//    Return a pointer to the viewer client attributes.
//
//  Returns:    A pointer to the viewer client attributes.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

AttributeSubject *
TensorViewerEnginePluginInfo::GetClientAtts()
{
    return clientAtts;
}

// ****************************************************************************
//  Method: TensorViewerEnginePluginInfo::GetDefaultAtts
//
//  Purpose:
//    Return a pointer to the viewer default attributes.
//
//  Returns:    A pointer to the viewer default attributes.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

AttributeSubject *
TensorViewerEnginePluginInfo::GetDefaultAtts()
{
    return defaultAtts;
}

// ****************************************************************************
//  Method: TensorViewerEnginePluginInfo::SetClientAtts
//
//  Purpose:
//    Set the viewer client attributes.
//
//  Arguments:
//    atts      A pointer to the new client attributes.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

void
TensorViewerEnginePluginInfo::SetClientAtts(AttributeSubject *atts)
{
    *clientAtts = *(TensorAttributes *)atts;
    clientAtts->Notify();
}

// ****************************************************************************
//  Method: TensorViewerEnginePluginInfo::GetClientAtts
//
//  Purpose:
//    Get the viewer client attributes.
//
//  Arguments:
//    atts      A pointer to return the client default attributes in.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

void
TensorViewerEnginePluginInfo::GetClientAtts(AttributeSubject *atts)
{
    *(TensorAttributes *)atts = *clientAtts;
}

// ****************************************************************************
//  Method: TensorViewerEnginePluginInfo::AllocAvtPlot
//
//  Purpose:
//    Return a pointer to a newly allocated avt plot.
//
//  Returns:    A pointer to the newly allocated avt plot.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

avtPlot *
TensorViewerEnginePluginInfo::AllocAvtPlot()
{
    return new avtTensorPlot;
}

// ****************************************************************************
//  Method: TensorViewerEnginePluginInfo::InitializePlotAtts
//
//  Purpose:
//    Initialize the plot attributes to the default attributes.
//
//  Arguments:
//    atts      The attribute subject to initialize.
//    plot      The viewer plot whose attributes are getting initialized.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

void
TensorViewerEnginePluginInfo::InitializePlotAtts(AttributeSubject *atts,
    const avtPlotMetaData &)
{
    *(TensorAttributes*)atts = *defaultAtts;
}
// ****************************************************************************
//  Method: TensorViewerEnginePluginInfo::GetMenuName
//
//  Purpose:
//    Return a pointer to the name to use in the viewer menus.
//
//  Returns:    A pointer to the name to use in the viewer menus.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

const char *
TensorViewerEnginePluginInfo::GetMenuName() const
{
    return "Tensor";
}

