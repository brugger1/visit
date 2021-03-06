// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#ifndef BOUNDARYOPATTRIBUTES_H
#define BOUNDARYOPATTRIBUTES_H
#include <AttributeSubject.h>


// ****************************************************************************
// Class: BoundaryOpAttributes
//
// Purpose:
//    Attributes for Boundary Operator
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class BoundaryOpAttributes : public AttributeSubject
{
public:
    // These constructors are for objects of this class
    BoundaryOpAttributes();
    BoundaryOpAttributes(const BoundaryOpAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    BoundaryOpAttributes(private_tmfs_t tmfs);
    BoundaryOpAttributes(const BoundaryOpAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~BoundaryOpAttributes();

    virtual BoundaryOpAttributes& operator = (const BoundaryOpAttributes &obj);
    virtual bool operator == (const BoundaryOpAttributes &obj) const;
    virtual bool operator != (const BoundaryOpAttributes &obj) const;
private:
    void Init();
    void Copy(const BoundaryOpAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();

    // Property setting methods
    void SetSmoothingLevel(int smoothingLevel_);

    // Property getting methods
    int GetSmoothingLevel() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;


    // IDs that can be used to identify fields in case statements
    enum {
        ID_smoothingLevel = 0,
        ID__LAST
    };

private:
    int smoothingLevel;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define BOUNDARYOPATTRIBUTES_TMFS "i"

#endif
