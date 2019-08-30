// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//                            avtPointMapper.h                               //
// ************************************************************************* //

#ifndef AVT_POINT_MAPPER_H
#define AVT_POINT_MAPPER_H

#include <avtMapper.h>
#include <plotter_exports.h>
#include <GlyphTypes.h>

#include <string>


class vtkLookupTable;

// ****************************************************************************
//  Class: avtPointMapper
//
//  Purpose:
//      A mapper for points.  This extends the functionality of a mapper by
//      mapping a glyph onto a dataset with a scalar variable.
//
//  Programmer: Kathleen Biagas
//  Creation:   August 17, 2016
//
//  Modifications:
//    Kathleen Biagas, Tue Aug 27 08:41:20 PDT 2019
//    Made ivars protected instead of private. Rename lut to pmLUT. (There are
//    mappers that derive both from this class and avtVariableMapper, which has
//    a protected member 'lut').
//
// ****************************************************************************

class PLOTTER_API  avtPointMapper : virtual public avtMapper
{
  public:
                               avtPointMapper();
    virtual                   ~avtPointMapper();


    void                       ColorByScalarOn(const std::string &);
    void                       ColorByScalarOff(void);

    void                       ScaleByVar(const std::string &);
    void                       DataScalingOn(const std::string &, int = 1);
    void                       DataScalingOff(void);

    void                       SetScale(double);
    void                       SetGlyphType(GlyphType);
    void                       SetPointSize(double s);

    virtual bool               SetFullFrameScaling(bool, const double *);

    bool                       GetColorByScalar(void)
                                  { return colorByScalar ; }
    void                       SetLUT(vtkLookupTable *);

  protected:
    virtual void               CustomizeMappers(void);
    virtual vtkDataSetMapper  *CreateMapper(void);

    bool                       colorByScalar;

    int                        spatialDim;
    double                     scale;
    std::string                scalingVarName;
    std::string                coloringVarName;
    int                        scalingVarDim;
    GlyphType                  glyphType;
    double                     pointSize;
    bool                       dataScaling;
    vtkLookupTable            *pmLUT;

};


#endif


