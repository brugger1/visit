// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//                       avtPseudocolorMapper.h                              //
// ************************************************************************* //

#ifndef AVT_PSEUDOCOLORMAPPER_H
#define AVT_PSEUDOCOLORMAPPER_H

#include <avtVariablePointGlyphMapper.h>

#include <string>
#include <vector>

// ****************************************************************************
//  Class:  avtPseudocolorMapper
//
//  Purpose:
//      Pseudocolor plot specific mapper, that utilizes a specialized
//      vtkDataSetMapper allowing Multiple representations of the same dataset
//      to be rendered at the same time( eg Surface, Wireframe, and Points).
//
//  Programmer: Kathleen Biagas
//  Creation:   August 24, 2016
//
//  Modifications:
//    Kathleen Biagas, Wed Apr 10 09:06:05 PDT 2019
//    Added pointSize.
//
//    Kathleen Biagas, Tue Aug 27 09:26:01 PDT 2019
//    Derive from avtVariablePointGlyphMapper.
//    Added SetLookupTable and labels.
//
//    Kathleen Biagas, Mon Oct 14 20:26:00 PDT 2019
//    Add TurnLightingOn, to override avtVariableMapper::TurnLightingOn.
//
// ****************************************************************************

class avtPseudocolorMapper : public avtVariablePointGlyphMapper
{
  public:
                               avtPseudocolorMapper();
    virtual                   ~avtPseudocolorMapper();

    void                       SetDrawSurface(bool);
    void                       SetDrawWireframe(bool);
    void                       SetDrawPoints(bool);
    void                       SetPointSize(int);
    void                       SetWireframeColor(double rgb[3]);
    void                       SetPointsColor(double rgb[3]);
    void                       SetLookupTable(vtkLookupTable *);

    void                       TurnLightingOn(void);

  protected:
    // these are called from avtMapper
    virtual vtkDataSetMapper  *CreateMapper(int);
    virtual void               CustomizeMappers(void);
    virtual void               SetLabels(std::vector<std::string> &, bool);

  private:

    bool   drawSurface;
    bool   drawWireframe;
    bool   drawPoints;
    int    pointSize;
    double wireframeColor[3];
    double pointsColor[3];
    std::vector<std::string> labels;
};

#endif

