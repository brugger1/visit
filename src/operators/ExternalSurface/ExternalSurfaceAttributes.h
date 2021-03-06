// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#ifndef EXTERNALSURFACEATTRIBUTES_H
#define EXTERNALSURFACEATTRIBUTES_H
#include <AttributeSubject.h>


// ****************************************************************************
// Class: ExternalSurfaceAttributes
//
// Purpose:
//    This class contains attributes for the external surface operator.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class ExternalSurfaceAttributes : public AttributeSubject
{
public:
    // These constructors are for objects of this class
    ExternalSurfaceAttributes();
    ExternalSurfaceAttributes(const ExternalSurfaceAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    ExternalSurfaceAttributes(private_tmfs_t tmfs);
    ExternalSurfaceAttributes(const ExternalSurfaceAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~ExternalSurfaceAttributes();

    virtual ExternalSurfaceAttributes& operator = (const ExternalSurfaceAttributes &obj);
    virtual bool operator == (const ExternalSurfaceAttributes &obj) const;
    virtual bool operator != (const ExternalSurfaceAttributes &obj) const;
private:
    void Init();
    void Copy(const ExternalSurfaceAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();

    // Property setting methods
    void SetRemoveGhosts(bool removeGhosts_);
    void SetEdgesIn2D(bool edgesIn2D_);

    // Property getting methods
    bool GetRemoveGhosts() const;
    bool GetEdgesIn2D() const;

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
        ID_removeGhosts = 0,
        ID_edgesIn2D,
        ID__LAST
    };

private:
    bool removeGhosts;
    bool edgesIn2D;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define EXTERNALSURFACEATTRIBUTES_TMFS "bb"

#endif
