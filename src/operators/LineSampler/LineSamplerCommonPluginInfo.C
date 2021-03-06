// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//  File: LineSamplerCommonPluginInfo.C
// ************************************************************************* //

#include <LineSamplerPluginInfo.h>
#include <LineSamplerAttributes.h>

// ****************************************************************************
//  Method: LineSamplerCommonPluginInfo::AllocAttributes
//
//  Purpose:
//    Return a pointer to a newly allocated attribute subject.
//
//  Returns:    A pointer to the newly allocated attribute subject.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

AttributeSubject *
LineSamplerCommonPluginInfo::AllocAttributes()
{
    return new LineSamplerAttributes;
}

// ****************************************************************************
//  Method: LineSamplerCommonPluginInfo::CopyAttributes
//
//  Purpose:
//    Copy a LineSampler attribute subject.
//
//  Arguments:
//    to        The destination attribute subject.
//    from      The source attribute subject.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

void 
LineSamplerCommonPluginInfo::CopyAttributes(AttributeSubject *to,
    AttributeSubject *from)
{
    *((LineSamplerAttributes *) to) = *((LineSamplerAttributes *) from);
}
