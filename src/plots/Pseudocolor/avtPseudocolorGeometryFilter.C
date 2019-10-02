// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//                   avtPseudocolorGeometryFilter.C                          //
// ************************************************************************* //

#include <avtPseudocolorGeometryFilter.h>

#include <avtCallback.h>
#include <DebugStream.h>

#include <vtkAppendPolyData.h>
#include <vtkCellData.h>
#include <vtkConeSource.h>
#include <vtkExtractCellsByType.h>
#include <vtkGeometryFilter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkRibbonFilter.h>
#include <vtkSphereSource.h>
#include <vtkTubeFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVertexFilter.h>

#include <string>

using std::string;

double GetBBoxSize( double *bbox )
{
    double vol = 1;
    int    numDims = 0;
    if (bbox[1] > bbox[0])
    {
        vol *= (bbox[1]-bbox[0]);
        numDims++;
    }
    if (bbox[3] > bbox[2])
    {
        vol *= (bbox[3]-bbox[2]);
        numDims++;
    }
    if (bbox[5] > bbox[4])
    {
        vol *= (bbox[5]-bbox[4]);
        numDims++;
    }

    double length = pow(vol, 1.0/numDims);
    return length;
}


// ****************************************************************************
//  Method: avtPseudocolorGeometryFilter constructor
//
//  Programmer: Kathleen Biagas
//  Creation:   August 20, 2019
//
//  Modifications:
//
// ****************************************************************************

avtPseudocolorGeometryFilter::avtPseudocolorGeometryFilter()
{
    bboxSize = 0;
}


// ****************************************************************************
//  Method: avtPseudocolorGeometryFilter destructor
//
//  Programmer: Kathleen Biagas
//  Creation:   August 20, 2019
//
//  Modifications:
//
// ****************************************************************************

avtPseudocolorGeometryFilter::~avtPseudocolorGeometryFilter()
{
}


// ****************************************************************************
//  Method: avtPseudocolorGeometryFilter::SetPlotAtts
//
//  Purpose:    Sets the PseudcolorAttributes needed for this filter.
//
//  Programmer: Kathleen Biagas
//  Creation:   August 20, 2019
//
// ****************************************************************************

void
avtPseudocolorGeometryFilter::SetPlotAtts(const PseudocolorAttributes *atts)
{
    plotAtts = *atts;
}


// ****************************************************************************
//  Method: avtPseudocolorGeometryFilter::ExecuteData
//
//  Purpose:
//      Applies tubes or ribbons to lines and and glyphs to line
//      end points (depending on input and plotAtts). Separates vertex
//      data into separate dataset for processing by the glyph mapper.
//
//  Arguments:
//      inDR    The input data representation.
//
//  Returns:    The output data tree.
//
//  Programmer: Kathleen Biagas
//  Creation:   August 20, 2019
//
// ****************************************************************************

avtDataTree_p
avtPseudocolorGeometryFilter::ExecuteDataTree(avtDataRepresentation *inDR)
{
    //
    // Get the VTK data set.
    //
    vtkDataSet *inDS = inDR->GetDataVTK();
    int domain       = inDR->GetDomain();
    string domString = std::to_string(domain);
    string label     = inDR->GetLabel();

    //
    // Verify input is of the right type
    //
    if ((inDS->GetDataObjectType() != VTK_POLY_DATA) &&
        (inDS->GetDataObjectType() != VTK_UNSTRUCTURED_GRID ))
    {
        return new avtDataTree(inDS, domain, label);
    }


    //
    // Find out whether we can/should apply glyphs.
    //
    bool canGlyphPoints = plotAtts.GetRenderSurfaces() ||
                          plotAtts.GetRenderPoints();

    bool glyphingLines = plotAtts.GetRenderSurfaces() &&
                         (plotAtts.GetLineType() != PseudocolorAttributes::Line);
    bool glyphingEnds  = plotAtts.GetRenderSurfaces() &&
                        !(plotAtts.GetTailStyle() == PseudocolorAttributes::None &&
                          plotAtts.GetHeadStyle() == PseudocolorAttributes::None);
    bool glyphingVerts = canGlyphPoints && plotAtts.GetPointType() != Point;
    int topoDim = GetInput()->GetInfo().GetAttributes().GetTopologicalDimension();


    if (!glyphingLines && !glyphingEnds && !glyphingVerts)
    {
        // Even if glyphing isn't being applied here, if the
        // input is all points, the glyph mapper should be used
        // instead of the regular surface mapper.  The labels associated
        // with the dataset help the PC mapper determine which vtk mapper
        // to utilize, so create the points label if appropriate

        if (topoDim == 0) // indicates points data
        {
            label = string("pc_points_") + label + domString;
        }
        else
        {
            if (inDS->GetDataObjectType() == VTK_POLY_DATA)
            {
                // Check for verts only, create label if so
                vtkPolyData *pd  = vtkPolyData::SafeDownCast(inDS);
                if (pd->GetNumberOfVerts() == pd->GetNumberOfCells())
                {
                    label = string("pc_points_") + label + domString;
                }
            }
            else  // unstructured grid
            {
                // Check for verts only, create label if so
                vtkUnstructuredGrid *ug  = vtkUnstructuredGrid::SafeDownCast(inDS);
                if (ug->IsHomogeneous() && (ug->GetCellType(0) == VTK_VERTEX ||
                    ug->GetCellType(0) == VTK_POLY_VERTEX))
                {
                    label = string("pc_points_") + label + domString;
                }
            }
        }
        return new avtDataTree(inDS, domain, label);
    }


    vtkPolyData *inPolys = NULL;
    if (inDS->GetDataObjectType() == VTK_POLY_DATA)
    {
        inPolys = vtkPolyData::SafeDownCast(inDS);
    }
    else  // already checked for polydata/ugrid, so this is the ugrid case
    {
        // convert it to polydata
        vtkNew<vtkGeometryFilter> geo;
        geo->SetInputData(inDS);
        geo->Update();
        inPolys = geo->GetOutput();
        inPolys->Register(NULL);
    }


    //
    // Verify input has right type of cells
    //
    if ((glyphingLines || glyphingEnds) && inPolys->GetNumberOfLines() < 1)
    {
        avtCallback::IssueWarning("Applying tubes/ribbons/endpoints only works"
            " on poly data or unstructured grids containing LINE cells.");
        glyphingLines = false;
        glyphingEnds = false;
    }
    if (!plotAtts.GetRenderPoints() &&
        (glyphingVerts && inPolys->GetNumberOfVerts() < 1))
    {
        avtCallback::IssueWarning("Applying glyphing to points only works"
            " on poly data or unstructured grids containing VERTEX cells."
            " If you want to glyph all the points, turn on "
            "'Draw objects as Points'");
        glyphingVerts = false;
    }
    if (!glyphingLines && !glyphingEnds && !glyphingVerts)
    {
        return new avtDataTree(inDS, domain, label);
    }

    //
    // If input contains more than lines, want to extract lines for processing
    //
    vtkPolyData *linesData = inPolys;
    bool removeLinesFromInput = false;
    if ((glyphingLines || glyphingEnds) &&
        inPolys->GetNumberOfLines() < inPolys->GetNumberOfCells())
    {
        vtkNew<vtkExtractCellsByType> extractLines;
        extractLines->AddCellType(VTK_LINE);
        extractLines->AddCellType(VTK_POLY_LINE);
        extractLines->SetInputData(inPolys);
        extractLines->Update();
        linesData = vtkPolyData::SafeDownCast(extractLines->GetOutput());
        linesData->Register(NULL);
        removeLinesFromInput = glyphingLines;
    }

    // Get bounding box size, used for scaling tubes, ribbons, line-end glyphs
    double bbox[6] = {0.,1.,0.,1.,0.,1.};
    GetInput()->GetInfo().GetAttributes().GetOriginalSpatialExtents()->CopyTo(bbox);
    double bboxSize = GetBBoxSize(bbox);

    //
    // Handle tubes or ribbons applied to lines
    //
    vtkNew<vtkPolyData> processedLines;
    if (plotAtts.GetLineType() == PseudocolorAttributes::Tube)
    {
        AddTubes(linesData, processedLines.GetPointer(), bboxSize);
    }
    else if (plotAtts.GetLineType() == PseudocolorAttributes::Ribbon)
    {
        AddRibbons(linesData, processedLines.GetPointer(), bboxSize);
    }

    //
    // Handle glyphs applied to line endpoints
    //
    vtkNew<vtkPolyData> endPoints;
    if (glyphingEnds)
    {
        AddEndPoints(linesData, endPoints, bboxSize);
    }

    //
    // If input contains Verts, we want to extract them into a separate
    // dataset for handling by the glyph mapper.
    // In a mixed-topology case (verts + other cell types), send the
    // input through an extractor.
    // In the simple case of all Verts, simply copy the input.
    //
    vtkNew<vtkPolyData> vertsData;

    bool removeVertsFromInput = false;
    if (glyphingVerts)
    {
        if (plotAtts.GetRenderPoints())
        {
            // Need to create a new dataset with vertex cells
            // for all the points.
            vtkNew<vtkVertexFilter> createVerts;
            createVerts->VertexAtPointsOn();
            createVerts->SetInputData(inPolys);
            createVerts->Update();
            vertsData->ShallowCopy(createVerts->GetOutput());
            removeVertsFromInput = true;
        }
        else if (inPolys->GetNumberOfVerts() < inPolys->GetNumberOfCells())
        {
            // Need to extract only vertex cells to a new dataset
            vtkNew<vtkExtractCellsByType> extractVerts;
            extractVerts->AddCellType(VTK_VERTEX);
            extractVerts->AddCellType(VTK_POLY_VERTEX);
            extractVerts->SetInputData(inPolys);
            extractVerts->Update();
            vertsData->ShallowCopy(extractVerts->GetOutput());
            removeVertsFromInput = true;
        }
    }
    else if (glyphingVerts && inPolys->GetNumberOfVerts() == inPolys->GetNumberOfCells())
    {
        vertsData->ShallowCopy(inPolys);
    }


    //
    //  Process the inputPolys into an output possibly devoid
    //  of lines and vertices.
    //
    vtkNew<vtkPolyData> outPolys;

    if (removeLinesFromInput || removeVertsFromInput)
    {
        vtkNew<vtkExtractCellsByType> remover;
        remover->SetInputData(inPolys);
        // tell the extractor we want all cell types ...
        remover->AddAllCellTypes();
        if(removeLinesFromInput)
        {
            // except Lines
            remover->RemoveCellType(VTK_LINE);
            remover->RemoveCellType(VTK_POLY_LINE);
        }
        if(removeVertsFromInput)
        {
            // except Verts
            remover->RemoveCellType(VTK_VERTEX);
            remover->RemoveCellType(VTK_POLY_VERTEX);
        }
        remover->Update();
        outPolys->ShallowCopy(remover->GetOutput());
    }
    else
    {
        // only want inPolys in output if it wasn't all verts
        if (!(glyphingVerts && !removeVertsFromInput))
            outPolys->ShallowCopy(inPolys);
    }

    //
    // append all glyphed output (tubes/ribbon, endpoints and) into a single dataset.
    //
    vtkNew<vtkPolyData> outLineGlyphs;
    if (glyphingLines || glyphingEnds)
    {
        vtkNew<vtkAppendPolyData>appender;
        if(glyphingLines)
            appender->AddInputData(processedLines);
        if (glyphingEnds)
            appender->AddInputData(endPoints);

        appender->Update();
        outLineGlyphs->ShallowCopy(appender->GetOutput());
    }

    //
    //  Want line and endpoint glyphs separate from polys because we don't
    //  want wireframe or points mode applied to them.
    //
    int numOutputs = outPolys->GetNumberOfCells() > 0 ? 1 : 0;
    numOutputs += outLineGlyphs->GetNumberOfCells() > 0 ? 1 : 0;
    numOutputs += vertsData->GetNumberOfCells() > 0 ? 1 : 0;

    // create the output data tree with datasets labeled
    avtDataTree_p rv = NULL;
    vtkDataSet **outs = new vtkDataSet *[numOutputs];
    stringVector l;
    int idx = 0;
    if (outPolys->GetNumberOfCells() > 0)
    {
        outs[idx++] = outPolys;
        // this label sends the dataset throught the vtkMultiRepMapper
        // and allows wireframe and points mode
        l.push_back(string("pc_polys_") + label + domString);
    }
    if (outLineGlyphs->GetNumberOfCells() > 0 )
    {
        outs[idx++] = outLineGlyphs;
        // this label sends the dataset throught the vtkMultiRepMapper
        // but does not allows wireframe and points mode
        l.push_back(string("pc_lineglyphs_") + label + domString);
    }
    if (vertsData->GetNumberOfCells() > 0)
    {
        outs[idx++] = vertsData;
        // this label sends the dataset throught the vtkPointGlyphMapper.
        l.push_back(string("pc_points_") + label + domString);
    }
    rv = new avtDataTree(numOutputs, outs, domain, l);
    delete [] outs;
    return rv;
}


// ****************************************************************************
//  Method: avtPseudocolorGeometryFilter::AddTubes
//
//  Purpose:  Applies vtkTubeFilter to input.
//
//  Notes:    Content taken from avtPolylineToTubeFilter.
//
//  Arguments:
//    input     The input. (Original input minus all but line cells).
//    ouput     A place to store the output.
//    bboxSize  The size of the bounding box.
//
//  Programmer: Kathleen Biagas
//  Creation:   August 20, 2019
//
//  Modifications:
//
// ****************************************************************************

void
avtPseudocolorGeometryFilter::AddTubes(vtkPolyData *input,
                                       vtkPolyData *output,
                                       double bboxSize)
{
    vtkNew<vtkTubeFilter> tubeFilter;
    tubeFilter->SetInputData(input);
    if (plotAtts.GetTubeRadiusSizeType() == PseudocolorAttributes::Absolute)
        tubeFilter->SetRadius(plotAtts.GetTubeRadiusAbsolute());
    else
        tubeFilter->SetRadius(plotAtts.GetTubeRadiusBBox() *bboxSize);
    tubeFilter->SetNumberOfSides(plotAtts.GetTubeResolution());
    tubeFilter->SetCapping(1);

    if (plotAtts.GetTubeRadiusVarEnabled())
    {
        string radiusVar = plotAtts.GetTubeRadiusVar();
        if (!radiusVar.empty() && radiusVar != "default")
        {
            int fieldAssociation = vtkDataObject::FIELD_ASSOCIATION_POINTS;
            if (input->GetCellData()->HasArray(radiusVar.c_str()))
            {
                fieldAssociation = vtkDataObject::FIELD_ASSOCIATION_CELLS;
            }
            tubeFilter->SetInputArrayToProcess(0, 0,0, fieldAssociation,
                                               radiusVar.c_str());
        }
        tubeFilter->SetVaryRadiusToVaryRadiusByScalar();
        tubeFilter->SetRadiusFactor(plotAtts.GetTubeRadiusVarRatio());
    }
    tubeFilter->Update();

    output->ShallowCopy(tubeFilter->GetOutput());
}


// ****************************************************************************
//  Method: avtPseudocolorGeometryFilter::AddRibbons
//
//  Purpose:  Applies vtkRibbonFilter to input, stores in output.
//
//  Notes:    Content taken from avtPolylineToRibbonFilter.
//
//  Arguments:
//  Arguments:
//    input     The input. (Original input minus all but line cells).
//    ouput     A place to store the output.
//    bboxSize  The size of the bounding box.
//
//  Programmer: Kathleen Biagas
//  Creation:   August 20, 2019
//
//  Modifications:
//
// ****************************************************************************

void
avtPseudocolorGeometryFilter::AddRibbons(vtkPolyData *input,
                                         vtkPolyData *output,
                                         double bboxSize)
{
    vtkNew<vtkRibbonFilter> ribbonFilter;
    ribbonFilter->SetInputData(input);
    if(plotAtts.GetTubeRadiusSizeType() == PseudocolorAttributes::Absolute)
        ribbonFilter->SetWidth(plotAtts.GetTubeRadiusAbsolute());
    else
        ribbonFilter->SetWidth(plotAtts.GetTubeRadiusBBox()*bboxSize);

    ribbonFilter->SetVaryWidth(plotAtts.GetTubeRadiusVarEnabled());

    if (plotAtts.GetTubeRadiusVarEnabled())
    {
        string widthVar(plotAtts.GetTubeRadiusVar());
        if (widthVar != "" && widthVar != "\0" && widthVar != "default")
        {
            int fieldAssociation = vtkDataObject::FIELD_ASSOCIATION_POINTS;
            if (input->GetCellData()->HasArray(widthVar.c_str()))
            {
                fieldAssociation = vtkDataObject::FIELD_ASSOCIATION_CELLS;
            }
            ribbonFilter->SetInputArrayToProcess(0, 0,0, fieldAssociation,
                                                 widthVar.c_str());
        }
    }
    ribbonFilter->Update();

    output->ShallowCopy(ribbonFilter->GetOutput());
}


// ****************************************************************************
//  Method: avtPseudocolorGeometryFilter::AddEndPoints
//
//  Purpose:  Applies glyphs to line endpoints.
//
//  Notes:    Content taken from avtPolylineAddEndPointsFilter.
//
//  Arguments:
//    input     The input. (Original input minus all but line cells).
//    ouput     A place to store the output.
//    bboxSize  The size of the bounding box.
//
//  Programmer: Kathleen Biagas
//  Creation:   August 20, 2019
//
//  Modifications:
//
// ****************************************************************************

void
avtPseudocolorGeometryFilter::AddEndPoints(vtkPolyData *input,
                                           vtkPolyData *output,
                                           double bboxSize)
{
    vtkNew<vtkAppendPolyData> appender;
    double radius = 0.;
    if( plotAtts.GetEndPointRadiusSizeType() == PseudocolorAttributes::Absolute )
        radius = plotAtts.GetEndPointRadiusAbsolute();
    else
        radius = plotAtts.GetEndPointRadiusBBox() * bboxSize;

    const avtDataAttributes &datts = GetInput()->GetInfo().GetAttributes();
    string activeVar = datts.GetVariableName();

    double ratio           = plotAtts.GetEndPointRatio();
    bool varyRadius        = plotAtts.GetEndPointRadiusVarEnabled();
    std::string radiusVar  = plotAtts.GetEndPointRadiusVar();
    double  radiusFactor   = plotAtts.GetEndPointRadiusVarRatio();
    int resolution         = plotAtts.GetEndPointResolution();

    vtkDataArray *radiusArray = NULL;
    double range[2] = {0,1}, scale = 1;

    if (varyRadius && !radiusVar.empty())
    {
        if (radiusVar == "default")
            radiusVar = activeVar;

        radiusArray = input->GetPointData()->GetArray(radiusVar.c_str());
        if (!radiusArray)
            radiusArray = input->GetCellData()->GetArray(radiusVar.c_str());

        radiusArray->GetRange(range, 0);

        if ((range[1] - range[0]) == 0.0)
            range[1] = range[0] + 1.0;

        scale = (radiusFactor-1) / (range[1]-range[0]);
    }

    vtkCellArray *lines  = input->GetLines();
    vtkPoints    *points = input->GetPoints();

    vtkIdType numPts;
    vtkIdType *ptIndexs;

    vtkIdType lineIndex = 0;
    lines->InitTraversal();

    vtkCellData  *inputCellData  = input->GetCellData();
    vtkPointData *inputPointData = input->GetPointData();

    while (lines->GetNextCell(numPts, ptIndexs))
    {
        vtkPolyData *outPD;

        double p0[3], p1[3];

        avtVector vec;

        // Do the two endpoints in a loop. The first iteration is for the
        // tail and the second is for the head.
        for (int i = 0; i < 2; ++i)
        {
            if ((i == 0 && plotAtts.GetTailStyle() != PseudocolorAttributes::None) ||
                (i == 1 && plotAtts.GetHeadStyle() != PseudocolorAttributes::None))
            {
                int style, tip, tail;

                if (i == 0)
                {
                    style = plotAtts.GetTailStyle();
                    tip  = 0;
                    tail = 1;
                }
                else
                {
                    style = plotAtts.GetHeadStyle();
                    tip  = numPts - 1;
                    tail = numPts - 2;
                }

                points->GetPoint(ptIndexs[tip], p0);

                double scaledRadius = radius;

                if (varyRadius && radiusArray)
                {
                    scaledRadius *=
                      (1.0 + (radiusArray->GetComponent( ptIndexs[tip], 0 ) -
                              range[0]) * scale);
                }

                if (style == PseudocolorAttributes::Spheres)
                {
                    vtkNew<vtkSphereSource> sphere;

                    sphere->SetRadius(scaledRadius);
                    sphere->SetPhiResolution(resolution);
                    sphere->SetThetaResolution(resolution);
                    sphere->SetCenter(p0);
                    sphere->Update();

                    outPD = sphere->GetOutput();
                    outPD->Register(NULL);

                }
                else if (style == PseudocolorAttributes::Cones)
                {
                    points->GetPoint(ptIndexs[tail], p1);

                    vec = avtVector(p0[0]-p1[0], p0[1]-p1[1],  p0[2]-p1[2]);
                    vec.normalize();

                    p0[0] += (vec * scaledRadius * ratio / 2.0).x;
                    p0[1] += (vec * scaledRadius * ratio / 2.0).y;
                    p0[2] += (vec * scaledRadius * ratio / 2.0).z;

                    vtkNew<vtkConeSource> cone;

                    cone->SetRadius(scaledRadius);
                    cone->SetHeight(scaledRadius * ratio);
                    cone->SetResolution(resolution);
                    cone->SetCenter(p0);
                    cone->SetDirection(vec.x, vec.y, vec.z);
                    cone->CappingOn();
                    cone->Update();

                    outPD = cone->GetOutput();
                    outPD->Register(NULL);
                }
                int npts   = outPD->GetNumberOfPoints();
                int ncells = outPD->GetNumberOfCells();

                vtkCellData  *outputCellData  = outPD->GetCellData();
                vtkPointData *outputPointData = outPD->GetPointData();

                // remove the generated normals
                outputPointData->RemoveArray("Normals");


                // Copy over all of the point data from the lines to the
                // glyph's points

                int nArrays = inputPointData->GetNumberOfArrays();
                for (int j = 0; j < nArrays; ++j)
                {
                    vtkDataArray *array = inputPointData->GetArray(j);
                    vtkDataArray *scalars = array->NewInstance();
                    scalars->Allocate(npts);
                    scalars->SetName(array->GetName());

                    outputPointData->AddArray(scalars);
                    if (array->GetName() == activeVar)
                    {
                        outputPointData->SetActiveScalars(activeVar.c_str());
                    }

                    double scalar = array->GetComponent(ptIndexs[tip], 0);

                    for (int k = 0; k < npts; ++k)
                        scalars->InsertTuple(k, array->GetTuple(ptIndexs[tip]));

                    scalars->Delete();
                }

                // Copy over all of the cell data from the lines to the
                // glyph's cells.
                nArrays = inputCellData->GetNumberOfArrays();
                for (int j = 0; j < nArrays; ++j)
                {
                    vtkDataArray *array = inputCellData->GetArray(j);
                    vtkDataArray *scalars = array->NewInstance();
                    scalars->Allocate(ncells);
                    scalars->SetName(array->GetName());
                    outputCellData->AddArray(scalars);
                    if (array->GetName() == activeVar)
                    {
                        outputCellData->SetActiveScalars(activeVar.c_str());
                    }

                    for (int k = 0; k < ncells; ++k)
                        scalars->InsertTuple(k, array->GetTuple(lineIndex));

                    scalars->Delete();
                }
                appender->AddInputData(outPD);

                outPD->Delete();
            }
        }
        ++lineIndex;
    }
    appender->Update();
    output->ShallowCopy(appender->GetOutput());
}


// ****************************************************************************
//  Method: avtPseudocolorGeometryFilter::UpdateDataObjectInfo
//
//  Purpose:  Sets flags in the pipeline.
//
//  Programmer: Kathleen Biagas
//  Creation:   August 20, 2019
//
//  Modifications:
//
// ****************************************************************************

void
avtPseudocolorGeometryFilter::UpdateDataObjectInfo(void)
{
    if (GetInput()->GetInfo().GetAttributes().GetTopologicalDimension() == 1)
        GetOutput()->GetInfo().GetValidity().InvalidateZones();
}


// ****************************************************************************
//  Method: avtPseudocolorGeometryFilter::PostExcecute
//
//  Purpose:
//    Sets the output's label attributes to reflect what is currently
//    present in the tree.
//
//  Programmer: Kathleen Biagas
//  Creation:   August 20, 2019
//
//  Modifications:
//
// ****************************************************************************

void
avtPseudocolorGeometryFilter::PostExecute(void)
{
    // Use labels to ensure lines/polys/verts aren't merged back together
    // during CompactTreeFilter
    stringVector treeLabels;
    GetDataTree()->GetAllUniqueLabels(treeLabels);
    GetOutput()->GetInfo().GetAttributes().SetLabels(treeLabels);
}

