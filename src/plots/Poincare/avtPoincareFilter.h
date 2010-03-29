/*****************************************************************************
*
* Copyright (c) 2000 - 2010, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-400124
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

// ************************************************************************* //
//                              avtPoincareFilter.h                              //
// ************************************************************************* //

#ifndef AVT_Poincare_FILTER_H
#define AVT_Poincare_FILTER_H


#include <avtStreamlineFilter.h>

// ****************************************************************************
//  Class: avtPoincareFilter
//
//  Purpose:
//      This operator is the implied operator associated with an Poincare plot.
//
//  Programmer: Dave Pugmire -- generated by xml2avt
//  Creation:   Tue Oct 7 09:02:52 PDT 2008
//
//  Modifications:
//    Dave Pugmire (for Allen Sanderson), Wed Feb 25 09:52:11 EST 2009
//    Add terminate by steps, add AdamsBashforth solver, Allen Sanderson's new code.
//
//    Dave Pugmire, Tue Apr 28 09:26:06 EDT 2009
//    Changed color to colorBy
//
//    Dave Pugmire, Wed May 27 15:03:42 EDT 2009
//    Removed GetStreamlinePoints().
//
//    Dave Pugmire, Tue Aug 18 09:10:49 EDT 2009
//    Add ability to restart streamline integration.
//
// ****************************************************************************

#include "StreamlineAnalyzerLib.h"


class avtPoincareFilter : public avtStreamlineFilter
{
  public:
    avtPoincareFilter();
    virtual                  ~avtPoincareFilter();

    virtual const char       *GetType(void) { return "avtPoincareFilter"; };
    virtual const char       *GetDescription(void) {
      return "Performing Poincare"; };

    void SetMaxPunctures( double punctures ) { maxPunctures = punctures; }

    void SetMaxToroidalWinding( unsigned int value ) {
      maxToroidalWinding = value; };
    void SetOverrideToroidalWinding( unsigned int value) { override = value; }

    void SetHitRate( double val ) { hitrate = val; }

    void SetAdjustPlane( int val ) { adjust_plane = val; }

    void SetOverlaps( unsigned int val ) { overlaps = val; }


    void SetShowCurves( unsigned int val ) { is_curvemesh = val; }

    void SetClipPlanes( vector< double > planeAngles ) {
      planes = planeAngles; }

    void SetDataValue( unsigned int val ) { dataValue = val; }

    void SetShowOPoints( bool val ) { showOPoints = val; }
    void SetShowIslands( bool val ) { showIslands = val; }
    void SetShowLines( bool val )   { showLines = val; }
    void SetShowPoints( bool val )  { showPoints = val; }
    void SetVerboseFlag( bool val ) { verboseFlag = val; }
    void SetShowRidgelines( bool val )   { showRidgelines = val; }

  protected:
    // Streamline overides.
    virtual void              Execute(void);
    virtual bool              ContinueExecute();
    virtual void              PreExecute(void);
    virtual void              PostExecute(void);
    virtual avtContract_p     ModifyContract(avtContract_p);
    virtual void              UpdateDataObjectInfo(void);
    virtual void              CreateStreamlineOutput( 
                                   vector<avtStreamlineWrapper *> &sls);

  virtual void findIslandCenter( avtDataTree *dt,
                                 vector< vector < vector < Point > > > &nodes,
                                 unsigned int color,
                                 double color_value );

  virtual void loadCurve( avtDataTree *dt,
                          vector< vector < vector < Point > > > &nodes,
                          unsigned int nnodes,
                          unsigned int islands,
                          unsigned int skip,
                          unsigned int color,
                          double color_value );
  
  virtual void loadCurve( avtDataTree *dt,
                          vector< vector < vector < Point > > > &nodes,
                          unsigned int nnodes,
                          unsigned int color,
                          double color_value );
  
  virtual void loadSurface( avtDataTree *dt,
                            vector< vector < vector < Point > > > &nodes,
                            unsigned int nnodes,
                            unsigned int islands,
                            unsigned int skip,
                            unsigned int color,
                            double color_value);

  virtual void loadPoints( avtDataTree *dt,
                           vector < Point  > &nodes,
                           unsigned int period,
                           unsigned int nnodes,
                           unsigned int islands,
                           unsigned int poloidalWindings,
                           unsigned int color,
                           double color_value,
                           bool ptFlag );

    // Poincare filter methods.
    bool                      ClassifyStreamlines();
    avtDataTree               *CreatePoincareOutput();


    FusionPSE::FieldlineLib FLlib;         

    double maxPunctures;
    unsigned int maxToroidalWinding;
    unsigned int override;

    double hitrate;
    int adjust_plane;
    unsigned int overlaps;

    bool is_curvemesh;
    std::vector< double > planes;

    unsigned int dataValue;
  
    bool showOPoints, showIslands;
    bool showLines, showPoints, showRidgelines, verboseFlag;

    class SLHelper
    {
      public:
        SLHelper() {}
        ~SLHelper() {}
        avtStreamlineWrapper *slSeg;
        std::vector<avtVector> streamlinePts;
    };

    std::vector<SLHelper> streamlines;
    std::vector<FieldlineInfo> poincareClassification;
};


#endif
