/*****************************************************************************
*
* Copyright (c) 2000 - 2011, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-442911
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

#include <ConnectedComponentsAttributes.h>
#include <DataNode.h>

// ****************************************************************************
// Method: ConnectedComponentsAttributes::ConnectedComponentsAttributes
//
// Purpose: 
//   Init utility for the ConnectedComponentsAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

void ConnectedComponentsAttributes::Init()
{
    EnableGhostNeighborsOptimization = true;

    ConnectedComponentsAttributes::SelectAll();
}

// ****************************************************************************
// Method: ConnectedComponentsAttributes::ConnectedComponentsAttributes
//
// Purpose: 
//   Copy utility for the ConnectedComponentsAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

void ConnectedComponentsAttributes::Copy(const ConnectedComponentsAttributes &obj)
{
    EnableGhostNeighborsOptimization = obj.EnableGhostNeighborsOptimization;

    ConnectedComponentsAttributes::SelectAll();
}

// Type map format string
const char *ConnectedComponentsAttributes::TypeMapFormatString = CONNECTEDCOMPONENTSATTRIBUTES_TMFS;
const AttributeGroup::private_tmfs_t ConnectedComponentsAttributes::TmfsStruct = {CONNECTEDCOMPONENTSATTRIBUTES_TMFS};


// ****************************************************************************
// Method: ConnectedComponentsAttributes::ConnectedComponentsAttributes
//
// Purpose: 
//   Default constructor for the ConnectedComponentsAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

ConnectedComponentsAttributes::ConnectedComponentsAttributes() : 
    AttributeSubject(ConnectedComponentsAttributes::TypeMapFormatString)
{
    ConnectedComponentsAttributes::Init();
}

// ****************************************************************************
// Method: ConnectedComponentsAttributes::ConnectedComponentsAttributes
//
// Purpose: 
//   Constructor for the derived classes of ConnectedComponentsAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

ConnectedComponentsAttributes::ConnectedComponentsAttributes(private_tmfs_t tmfs) : 
    AttributeSubject(tmfs.tmfs)
{
    ConnectedComponentsAttributes::Init();
}

// ****************************************************************************
// Method: ConnectedComponentsAttributes::ConnectedComponentsAttributes
//
// Purpose: 
//   Copy constructor for the ConnectedComponentsAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

ConnectedComponentsAttributes::ConnectedComponentsAttributes(const ConnectedComponentsAttributes &obj) : 
    AttributeSubject(ConnectedComponentsAttributes::TypeMapFormatString)
{
    ConnectedComponentsAttributes::Copy(obj);
}

// ****************************************************************************
// Method: ConnectedComponentsAttributes::ConnectedComponentsAttributes
//
// Purpose: 
//   Copy constructor for derived classes of the ConnectedComponentsAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

ConnectedComponentsAttributes::ConnectedComponentsAttributes(const ConnectedComponentsAttributes &obj, private_tmfs_t tmfs) : 
    AttributeSubject(tmfs.tmfs)
{
    ConnectedComponentsAttributes::Copy(obj);
}

// ****************************************************************************
// Method: ConnectedComponentsAttributes::~ConnectedComponentsAttributes
//
// Purpose: 
//   Destructor for the ConnectedComponentsAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

ConnectedComponentsAttributes::~ConnectedComponentsAttributes()
{
    // nothing here
}

// ****************************************************************************
// Method: ConnectedComponentsAttributes::operator = 
//
// Purpose: 
//   Assignment operator for the ConnectedComponentsAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

ConnectedComponentsAttributes& 
ConnectedComponentsAttributes::operator = (const ConnectedComponentsAttributes &obj)
{
    if (this == &obj) return *this;

    ConnectedComponentsAttributes::Copy(obj);

    return *this;
}

// ****************************************************************************
// Method: ConnectedComponentsAttributes::operator == 
//
// Purpose: 
//   Comparison operator == for the ConnectedComponentsAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

bool
ConnectedComponentsAttributes::operator == (const ConnectedComponentsAttributes &obj) const
{
    // Create the return value
    return ((EnableGhostNeighborsOptimization == obj.EnableGhostNeighborsOptimization));
}

// ****************************************************************************
// Method: ConnectedComponentsAttributes::operator != 
//
// Purpose: 
//   Comparison operator != for the ConnectedComponentsAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

bool
ConnectedComponentsAttributes::operator != (const ConnectedComponentsAttributes &obj) const
{
    return !(this->operator == (obj));
}

// ****************************************************************************
// Method: ConnectedComponentsAttributes::TypeName
//
// Purpose: 
//   Type name method for the ConnectedComponentsAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

const std::string
ConnectedComponentsAttributes::TypeName() const
{
    return "ConnectedComponentsAttributes";
}

// ****************************************************************************
// Method: ConnectedComponentsAttributes::CopyAttributes
//
// Purpose: 
//   CopyAttributes method for the ConnectedComponentsAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

bool
ConnectedComponentsAttributes::CopyAttributes(const AttributeGroup *atts)
{
    if(TypeName() != atts->TypeName())
        return false;

    // Call assignment operator.
    const ConnectedComponentsAttributes *tmp = (const ConnectedComponentsAttributes *)atts;
    *this = *tmp;

    return true;
}

// ****************************************************************************
// Method: ConnectedComponentsAttributes::CreateCompatible
//
// Purpose: 
//   CreateCompatible method for the ConnectedComponentsAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

AttributeSubject *
ConnectedComponentsAttributes::CreateCompatible(const std::string &tname) const
{
    AttributeSubject *retval = 0;
    if(TypeName() == tname)
        retval = new ConnectedComponentsAttributes(*this);
    // Other cases could go here too. 

    return retval;
}

// ****************************************************************************
// Method: ConnectedComponentsAttributes::NewInstance
//
// Purpose: 
//   NewInstance method for the ConnectedComponentsAttributes class.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

AttributeSubject *
ConnectedComponentsAttributes::NewInstance(bool copy) const
{
    AttributeSubject *retval = 0;
    if(copy)
        retval = new ConnectedComponentsAttributes(*this);
    else
        retval = new ConnectedComponentsAttributes;

    return retval;
}

// ****************************************************************************
// Method: ConnectedComponentsAttributes::SelectAll
//
// Purpose: 
//   Selects all attributes.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

void
ConnectedComponentsAttributes::SelectAll()
{
    Select(ID_EnableGhostNeighborsOptimization, (void *)&EnableGhostNeighborsOptimization);
}

///////////////////////////////////////////////////////////////////////////////
// Set property methods
///////////////////////////////////////////////////////////////////////////////

void
ConnectedComponentsAttributes::SetEnableGhostNeighborsOptimization(bool EnableGhostNeighborsOptimization_)
{
    EnableGhostNeighborsOptimization = EnableGhostNeighborsOptimization_;
    Select(ID_EnableGhostNeighborsOptimization, (void *)&EnableGhostNeighborsOptimization);
}

///////////////////////////////////////////////////////////////////////////////
// Get property methods
///////////////////////////////////////////////////////////////////////////////

bool
ConnectedComponentsAttributes::GetEnableGhostNeighborsOptimization() const
{
    return EnableGhostNeighborsOptimization;
}

///////////////////////////////////////////////////////////////////////////////
// Keyframing methods
///////////////////////////////////////////////////////////////////////////////

// ****************************************************************************
// Method: ConnectedComponentsAttributes::GetFieldName
//
// Purpose: 
//   This method returns the name of a field given its index.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

std::string
ConnectedComponentsAttributes::GetFieldName(int index) const
{
    switch (index)
    {
    case ID_EnableGhostNeighborsOptimization: return "EnableGhostNeighborsOptimization";
    default:  return "invalid index";
    }
}

// ****************************************************************************
// Method: ConnectedComponentsAttributes::GetFieldType
//
// Purpose: 
//   This method returns the type of a field given its index.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

AttributeGroup::FieldType
ConnectedComponentsAttributes::GetFieldType(int index) const
{
    switch (index)
    {
    case ID_EnableGhostNeighborsOptimization: return FieldType_bool;
    default:  return FieldType_unknown;
    }
}

// ****************************************************************************
// Method: ConnectedComponentsAttributes::GetFieldTypeName
//
// Purpose: 
//   This method returns the name of a field type given its index.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

std::string
ConnectedComponentsAttributes::GetFieldTypeName(int index) const
{
    switch (index)
    {
    case ID_EnableGhostNeighborsOptimization: return "bool";
    default:  return "invalid index";
    }
}

// ****************************************************************************
// Method: ConnectedComponentsAttributes::FieldsEqual
//
// Purpose: 
//   This method compares two fields and return true if they are equal.
//
// Note:       Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

bool
ConnectedComponentsAttributes::FieldsEqual(int index_, const AttributeGroup *rhs) const
{
    const ConnectedComponentsAttributes &obj = *((const ConnectedComponentsAttributes*)rhs);
    bool retval = false;
    switch (index_)
    {
    case ID_EnableGhostNeighborsOptimization:
        {  // new scope
        retval = (EnableGhostNeighborsOptimization == obj.EnableGhostNeighborsOptimization);
        }
        break;
    default: retval = false;
    }

    return retval;
}

///////////////////////////////////////////////////////////////////////////////
// User-defined methods.
///////////////////////////////////////////////////////////////////////////////
