// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//                       avtPseudocolorMapper.C                              //
// ************************************************************************* //

#include <avtPseudocolorMapper.h>

#include <vtkActor.h>
#include <vtkDataSet.h>
#include <vtkMultiRepMapper.h>
#include <vtkPointGlyphMapper.h>
#include <vtkProperty.h>

using std::string;
using std::vector;

// ****************************************************************************
//  Method: avtPseudocolorMapper constructor
//
//  Programmer: Kathleen Biagas
//  Creation:   Autust 24, 2016
//
//  Modifications:
//
// ****************************************************************************

avtPseudocolorMapper::avtPseudocolorMapper() : avtVariableMapper()
{
    drawSurface   = true;
    drawWireframe = false;
    drawPoints    = false;
    pointSize = 2;
    wireframeColor[0] = wireframeColor[1] = wireframeColor[2] = 0.;
    pointsColor[0] = pointsColor[1] = pointsColor[2] = 0.;
}


// ****************************************************************************
//  Method: avtPseudocolorMapper destructor
//
//  Programmer: Kathleen Biagas
//  Creation:   August 24, 2016
//
// ****************************************************************************

avtPseudocolorMapper::~avtPseudocolorMapper()
{
}


// ****************************************************************************
//  Method: avtPseudocolorMapper::CreateMapper
//
//  Purpose:
//    Creates a vtkMultiRepMapper or vtkPointGlyphMapper, depending on
//    input index.
//
//  Arguments:
//    index     The index for the mapper to be created.  Allows choosing type
//              of  mapper based on label assigned to the index'd dataset.
//
//  Programmer: Kathleen Biagas
//  Creation:   August 24, 2016
//
// ****************************************************************************

vtkDataSetMapper *
avtPseudocolorMapper::CreateMapper(int index)
{
    if (labels.size() == nMappers &&  index >= 0 && index < labels.size())
    {
        if (labels[index].compare(0, 10, string("pc_points_")) == 0)
            return ((vtkDataSetMapper*)vtkPointGlyphMapper::New());
    }
    return vtkMultiRepMapper::New();
}


// ****************************************************************************
//  Method: avtPseudocolorMapper::CustomizeMappers
//
//  Purpose:
//    Adds our flags to the vtk mapper.
//
//  Programmer: Kathleen Biagas
//  Creation:   August 24, 2016
//
//  Modifications:
//    Kathleen Biagas, Tue Aug 27 09:19:30 PDT 2019
//    Use label to determine which underlying mapper and settings to use.
//
// ****************************************************************************

void
avtPseudocolorMapper::CustomizeMappers()
{
    avtVariableMapper::CustomizeMappers();

    for (int i = 0; i < nMappers; ++i)
    {
        if (mappers[i] == NULL)
            continue;
        if(labels.empty() || labels[i].compare(0, 9, string("pc_polys_")) == 0)
        {
            vtkMultiRepMapper *mrm = (vtkMultiRepMapper*)mappers[i];
            mrm->SetDrawSurface(drawSurface);
            mrm->SetDrawWireframe(drawWireframe);
            mrm->SetDrawPoints(drawPoints);
            mrm->SetWireframeColor(wireframeColor);
            mrm->SetPointsColor(pointsColor);
        }
        else if (!labels.empty() && labels[i].compare(0, 10, string("pc_points_")) == 0)
        {
            vtkPointGlyphMapper *pm = (vtkPointGlyphMapper*)mappers[i];
            pm->SetSpatialDimension(
                GetInput()->GetInfo().GetAttributes().GetSpatialDimension());
            pm->SetLookupTable(avtPointMapper::pmLUT);
            pm->SetGlyphType(avtPointMapper::glyphType);
            if(drawPoints)
                pm->ColorByScalarOff();
            else
                pm->ColorByScalarOn(avtPointMapper::coloringVarName);
            actors[i]->GetProperty()->SetColor(pointsColor);
        }
        actors[i]->GetProperty()->SetPointSize(pointSize);
    }
    if (dataScaling)
    {
        // need to call ScaleByVar rather than DataScalingOn, because
        // scalingVarDim may not yet have been set correctly.
        ScaleByVar(avtPointMapper::scalingVarName);
    }
    else
        DataScalingOff();

    SetScale(avtPointMapper::scale);
}


// ****************************************************************************
//  Method: avtPseudocolorMapper::SetDrawSurface
//
//  Purpose:
//     Toggles the surface representation mode
//
//  Programmer: Kathleen Biagas
//  Creation:   August 24, 2016
//
//  Modifications:
//    Kathleen Biagas, Tue Aug 20 16:59:20 PDT 2019
//    Only apply this setting to the non-glyphed datasets.
//
// ****************************************************************************

void
avtPseudocolorMapper::SetDrawSurface(bool val)
{
    if (drawSurface != val)
    {
        drawSurface = val;
        for (int i = 0; i < nMappers; ++i)
        {
            if (mappers[i] == NULL)
                continue;

            if (labels.empty() || labels[i].compare(0, 9, string("pc_polys_")) == 0)
                ((vtkMultiRepMapper *)mappers[i])->SetDrawSurface(drawSurface);
        }
    }
}


// ****************************************************************************
//  Method: avtPseudocolorMapper::SetDrawWireframe
//
//  Purpose:
//     Toggles the Wireframe representation mode
//
//  Programmer: Kathleen Biagas
//  Creation:   August 24, 2016
//
//  Modifications:
//    Kathleen Biagas, Tue Aug 20 16:59:20 PDT 2019
//    Only apply this setting to the non-glyphed datasets.
//
// ****************************************************************************

void
avtPseudocolorMapper::SetDrawWireframe(bool val)
{
    if (drawWireframe != val)
    {
        drawWireframe = val;
        for (int i = 0; i < nMappers; ++i)
        {
            if (mappers[i] == NULL)
                continue;

             if (labels.empty() || labels[i].compare(0, 9, string("pc_polys_")) == 0)
                ((vtkMultiRepMapper *)mappers[i])->SetDrawWireframe(drawWireframe);
        }
    }
}


// ****************************************************************************
//  Method: avtPseudocolorMapper::SetDrawPoints
//
//  Purpose:
//     Toggles the Points representation mode
//
//  Programmer: Kathleen Biagas
//  Creation:   August 24, 2016
//
//  Modifications:
//    Kathleen Biagas, Tue Aug 20 16:59:20 PDT 2019
//    Only apply this setting to the non-glyphed datasets.
//
// ****************************************************************************

void
avtPseudocolorMapper::SetDrawPoints(bool val)
{
    if (drawPoints != val)
    {
        drawPoints = val;
        for (int i = 0; i < nMappers; ++i)
        {
            if (mappers[i] == NULL)
                continue;

            if (labels.empty() || labels[i].compare(0, 9, string("pc_polys_")) == 0)
                ((vtkMultiRepMapper *)mappers[i])->SetDrawPoints(drawPoints);
        }
    }
}


// ****************************************************************************
//  Method: avtPseudocolorMapper::SetPointSize
//
//  Purpose:
//     Sets the point size.
//
//  Programmer: Kathleen Biagas
//  Creation:   April 10, 2019
//
// ****************************************************************************

void
avtPseudocolorMapper::SetPointSize(int ps)
{
    if (pointSize != ps)
    {
        pointSize = ps;
        for (int i = 0; i < nMappers; ++i)
        {
            if (actors[i] != NULL)
                actors[i]->GetProperty()->SetPointSize(pointSize);
        }
    }
}


// ****************************************************************************
//  Method: ColorsAreDifferent
//
//  Purpose:
//     Helper method for comparing rgb colors.
//
//  Programmer: Kathleen Biagas
//  Creation:   June 30, 2016
//
// ****************************************************************************

bool
ColorsAreDifferent(double a[3], double b[3])
{
   return ((a[0] != b[0]) ||
           (a[1] != b[1]) ||
           (a[2] != b[2]));
}


// ****************************************************************************
//  Method: avtPseudocolorMapper::SetWireframeColor
//
//  Purpose:
//     Sets color to be used for the wirefame mode
//
//  Programmer: Kathleen Biagas
//  Creation:   August 24, 2016
//
//  Modifications:
//    Kathleen Biagas, Tue Aug 27 09:21:10 PDT 2019
//    Check type of mapper.
//
// ****************************************************************************

void
avtPseudocolorMapper::SetWireframeColor(double rgb[3])
{
    if (ColorsAreDifferent(wireframeColor, rgb))
    {
        wireframeColor[0] = rgb[0];
        wireframeColor[1] = rgb[1];
        wireframeColor[2] = rgb[2];
        for (int i = 0; i < nMappers; ++i)
        {
            if (mappers[i] != NULL && mappers[i]->IsA("vtkMultiRepMapper"))
                ((vtkMultiRepMapper *)mappers[i])->SetWireframeColor(wireframeColor);
        }
    }
}


// ****************************************************************************
//  Method: avtPseudocolorMapper::SetPointsColor
//
//  Purpose:
//     Sets color to be used for the wirefame mode
//
//  Programmer: Kathleen Biagas
//  Creation:   August 24, 2016
//
//  Modifications:
//    Kathleen Biagas, Tue Aug 27 09:21:10 PDT 2019
//    Check type of mapper.
//
// ****************************************************************************

void
avtPseudocolorMapper::SetPointsColor(double rgb[3])
{
    if (ColorsAreDifferent(pointsColor, rgb))
    {
        pointsColor[0] = rgb[0];
        pointsColor[1] = rgb[1];
        pointsColor[2] = rgb[2];
        for (int i = 0; i < nMappers; ++i)
        {
            if (mappers[i] != NULL)
            {
                 if (mappers[i]->IsA("vtkMultiRepMapper"))
                 {
                     ((vtkMultiRepMapper *)mappers[i])->SetPointsColor(pointsColor);
                 }
                 else //if (mappers[i]->IsA("vtkPointGlyphMapper"))
                 {
                     actors[i]->GetProperty()->SetColor(pointsColor);
                 }
            }
        }
    }
}


// ****************************************************************************
//  Method: avtPseudocolorMapper::SetLabels
//
//  Purpose:
//     Saves the labels coming from the input data tree.
//
//  Programmer: Kathleen Biagas
//  Creation:   August 20, 2019
//
// ****************************************************************************

void
avtPseudocolorMapper::SetLabels(vector<string> &l, bool fromTree)
{
    if (!fromTree)
        return;

    labels = l;
}


// ****************************************************************************
//  Method: avtPseudocolorMapper::SetLookupTable
//
//  Purpose:
//     Sends the lookuptable to both parent classes.
//
//  Programmer: Kathleen Biagas
//  Creation:   August 26, 2019
//
// ****************************************************************************

void
avtPseudocolorMapper::SetLookupTable(vtkLookupTable *LUT)
{
    avtVariableMapper::SetLookupTable(LUT);
    avtPointMapper::SetLUT(LUT);
}

