// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//                           avtExodusOptions.C                              //
// ************************************************************************* //

#include <avtExodusOptions.h>
#include <DBOptionsAttributes.h>

#include <vector>
#include <string>

using std::vector;
using std::string;
using namespace ExodusDBOptions;

// ****************************************************************************
//  Function: GetExodusReadOptions
//
//  Purpose:
//      Creates the options for Exodus readers.
//
//  Programmer: miller -- generated by xml2avt
//  Mark C. Miller, Tue Dec  9 10:04:22 PST 2014
//
//    Mark C. Miller, Wed Feb 11 17:05:25 PST 2015
//    Elaborated on help string for defining expressions.
// ****************************************************************************

DBOptionsAttributes *
GetExodusReadOptions(void)
{
    DBOptionsAttributes *rv = new DBOptionsAttributes;

    rv->SetBool(EXODUS_DETECT_COMPOUND_VARS, true);

    vector<string> materialConvention;
    materialConvention.push_back("None");   // 0
    materialConvention.push_back("ALEGRA"); // 1
    materialConvention.push_back("CTH");    // 2
    materialConvention.push_back("Custom"); // 3
    rv->SetEnum(EXODUS_MATERIAL_CONVENTION, 0); // None
    rv->SetEnumStrings(EXODUS_MATERIAL_CONVENTION, materialConvention);
    rv->SetInt(EXODUS_MATERIAL_COUNT, -1);
    rv->SetString(EXODUS_VOLFRAC_NAMESCHEME, "");
    rv->SetString(EXODUS_MATSPEC_NAMESCHEME, "");

    char helpStr[4096];
    snprintf(helpStr, sizeof(helpStr),
        "<p><b>%s</b>: Checking this option will cause the plugin to try to guess that similarly "
        "named variables are the scalar components of an aggregate type such as a vector, "
        "tensor or array variable. The plugin will then automatically define expressions "
        "for these aggregate typed variables. Note that this is just a convenience to free "
        "users from having to define expressions manally within their VisIt session."
        "<p> "
        "<p><b>%s</b>: Ordinarily, the plugin will determine the material count from the "
        "material convention nameschemes. However, if it is having trouble getting the "
        "correct count, you can specify it manually with this option. "
        "<p><b>%s</b>: A few pre-defined conventions for handling mixed materials from Exodus files "
        "are supported. In addition, you can define your own custom conventions as well. "
        "For a custom convention, you must define the <i>namescheme</i> that will produce "
        "the names of the scalar variables holding material volume fractions. Optionally, "
        "you can specify a <i>namescheme</i> to produce the names of the scalar variables "
        "holding material-specific values given the name of a non-material-specific variable. "
        "For more information on nameschemes, please consult the description of DBMakeNamescheme "
        "in the <a href=\"https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/"
        "silo/LLNL-SM-654357.pdf#page=226\">Silo user's manual</a>. The nameschemes used here are identical "
        "to those described in the Silo user's manual with one extension. The conversion specifier %%V "
        "is used to denote the basename (non-material-specific) name of a set of scalar variables "
        "holding material specific values."
        "<p> "
        "<p>The ALEGRA nameschemes for volume fraction and material specific variables  are "
        "\"%s\" and \"%s\"."
        "<p> "
        "<p>The CTH nameschemes are \"%s\" and \"%s\"."
        "<p> "
        "<p>Finally, it is assumed materials are identified starting from one (1). The special "
        "material id of zero (0) is used to denote void.",
            EXODUS_DETECT_COMPOUND_VARS,
            EXODUS_MATERIAL_COUNT, 
            EXODUS_MATERIAL_CONVENTION,
            EXODUS_VOLFRAC_NAMESCHEME_ALEGRA,
            EXODUS_MATSPEC_NAMESCHEME_ALEGRA,
            EXODUS_VOLFRAC_NAMESCHEME_CTH,
            EXODUS_MATSPEC_NAMESCHEME_CTH);
    rv->SetHelp(helpStr);

    return rv;
}


// ****************************************************************************
//  Function: GetExodusWriteOptions
//
//  Purpose:
//      Creates the options for Exodus writers.
//
//  Programmer: miller -- generated by xml2avt
//  Mark C. Miller, Tue Dec  9 10:04:22 PST 2014
//
// ****************************************************************************

DBOptionsAttributes *
GetExodusWriteOptions(void)
{
    DBOptionsAttributes *rv = new DBOptionsAttributes;
    return rv;
}
