// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//  File: RemapCommonPluginInfo.C
// ************************************************************************* //

#include <RemapPluginInfo.h>
#include <RemapAttributes.h>

// ****************************************************************************
//  Method: RemapCommonPluginInfo::AllocAttributes
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
RemapCommonPluginInfo::AllocAttributes()
{
    return new RemapAttributes;
}

// ****************************************************************************
//  Method: RemapCommonPluginInfo::CopyAttributes
//
//  Purpose:
//    Copy a Remap attribute subject.
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
RemapCommonPluginInfo::CopyAttributes(AttributeSubject *to,
    AttributeSubject *from)
{
    *((RemapAttributes *) to) = *((RemapAttributes *) from);
}
