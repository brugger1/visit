// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#include "vtkMultiRepMapper.h"

#include <vtkActor.h>
#include <vtkObjectFactory.h>
#include <vtkProperty.h>

//-----------------------------------------------------------------------------
vtkStandardNewMacro(vtkMultiRepMapper)

//-----------------------------------------------------------------------------
vtkMultiRepMapper::vtkMultiRepMapper()
{
  this->DrawSurface = true;
  this->DrawWireframe = false;
  this->DrawPoints = false;
  this->WireframeColor[0] = 1.;
  this->WireframeColor[1] = 1.;
  this->WireframeColor[2] = 1.;
  this->PointsColor[0] = 1.;
  this->PointsColor[1] = 1.;
  this->PointsColor[2] = 1.;
  this->currentScalarVis = this->ScalarVisibility;
}

//-----------------------------------------------------------------------------
vtkMultiRepMapper::~vtkMultiRepMapper()
{
}

//----------------------------------------------------------------------------
void vtkMultiRepMapper::Render(vtkRenderer *ren, vtkActor *act)
{
  if (this->DrawSurface)
    {
    act->GetProperty()->SetRepresentationToSurface();
    this->Superclass::Render(ren,act);
    }
  if (this->DrawWireframe)
    {
    bool sv = this->ScalarVisibility;
    this->ScalarVisibilityOff();
    act->GetProperty()->SetRepresentationToWireframe();
    act->GetProperty()->SetColor(this->WireframeColor);
    this->Superclass::Render(ren,act);
    this->SetScalarVisibility(sv);
    }
  if (this->DrawPoints)
    {
    bool sv = this->ScalarVisibility;
    this->ScalarVisibilityOff();
    act->GetProperty()->SetRepresentationToPoints();
    act->GetProperty()->SetColor(this->PointsColor);
    this->Superclass::Render(ren,act);
    this->SetScalarVisibility(sv);
    }
}

//-----------------------------------------------------------------------------
void vtkMultiRepMapper::PrintSelf(ostream& os, vtkIndent indent)
{
  os << "Draw surface: "  << this->DrawSurface << endl;
  os << "Draw wireframe: " << this->DrawWireframe << endl;
  os << "Draw points: "    << this->DrawPoints << endl;
  os << "Wireframe color: " << this->WireframeColor[0]  << " "
                            << this->WireframeColor[1]  << " "
                            << this->WireframeColor[2] ;
  os << "Points color: " << this->PointsColor[0]  << " "
                         << this->PointsColor[1]  << " "
                         << this->PointsColor[2] ;
  this->Superclass::PrintSelf(os, indent);

}

//-----------------------------------------------------------------------------
void
vtkMultiRepMapper::TurnLightingOn(vtkProperty *prop)
{
    if(this->DrawSurface)
    {
        // This method can get called multiple times.  If there are
        // multiple lights, then blindly calling SetAmbient to 0.0 may
        // turn off an ambient light that augments a normal light.  So
        // only set the default lighting attributes if we know we are
        // in "unlit" mode, which would mean diffuse would be not 1.0.
        if (prop->GetDiffuse() != 1.0)
        {
            prop->SetAmbient(0.0);
            prop->SetDiffuse(1.0);
        }
    }
    else // don't want lighting if only drawing lines or points
    {
        prop->SetAmbient(1.0);
        prop->SetDiffuse(0.0);
    }
}
