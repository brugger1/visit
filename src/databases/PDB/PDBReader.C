// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#include <DebugStream.h>

// This header file is last because it includes "scstd.h" (indirectly
// through "pdb.h"), which defines min and max, which conflict with
// "limits.h" on tru64 and aix systems.  On tru64 systems the conflict
// occurs with gcc-3.0.4 and on aix systems the conflict occurs with
// gcc-3.1.
#include <PDBReader.h>

// ****************************************************************************
//
// Method: PDBReader::PDBReader
//
// Purpose: 
//   Constructor for the PDBReader class.
//
// Arguments:
//   p : A pointer to the PDBfile object.
//
// Programmer: Brad Whitlock
// Creation:   Thu Oct 10 08:58:35 PDT 2002
//
// Modifications:
//
// ****************************************************************************

PDBReader::PDBReader(PDBFileObject *p)
{
    pdb = p;
    ownsPDBFile = false;
}

PDBReader::PDBReader(const char *filename)
{
    pdb = new PDBFileObject(filename);
    ownsPDBFile = true;
}

// ****************************************************************************
// Method: PDBReader::~PDBReader
//
// Purpose:
//   Destructor for the PDBReader class.
//
// Programmer: Brad Whitlock
// Creation:   Thu Oct 10 08:59:02 PDT 2002
//
// Modifications:
//   
// ****************************************************************************

PDBReader::~PDBReader()
{
    Close();
}

// ****************************************************************************
// Method: PDBReader::Close
//
// Purpose: 
//   Closes the reader's PDB file.
//
// Programmer: Brad Whitlock
// Creation:   Tue Sep 16 10:22:37 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

void
PDBReader::Close()
{
    if(ownsPDBFile && pdb != 0)
    {
        delete pdb;
        pdb = 0;
        ownsPDBFile = false;
    }
}

// ****************************************************************************
// Method: PDBReader::SetOwnsPDBFile
//
// Purpose: 
//   Sets a flag indicating that the reader owns its PDB file and is
//   responsible for closing it.
//
// Arguments:
//   v : The new flag value.
//
// Programmer: Brad Whitlock
// Creation:   Tue Sep 16 10:23:03 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

void
PDBReader::SetOwnsPDBFile(bool v)
{
    ownsPDBFile = v;
}

// ****************************************************************************
// Method: PDBReader::Identify
//
// Purpose: 
//   Returns true if the PDB file can be read by this reader.
//
// Returns:    True if the PDB can be read; false otherwise.
//
// Programmer: Brad Whitlock
// Creation:   Tue Sep 16 10:23:44 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

bool
PDBReader::Identify()
{
    // Make sure that the file is open.
    if(!pdb->IsOpen())
        pdb->Open();

    // Call the derived types' method for identifying the file.
    return IdentifyFormat();
}

// ****************************************************************************
// Method: PDBReader::PDB
//
// Purpose: 
//   Return the object's PDB file.
//
// Programmer: Brad Whitlock
// Creation:   Thu Nov  6 14:43:02 PST 2008
//
// Modifications:
//   
//   Mark C. Miller, Tue Apr 28 11:05:54 PDT 2009
//   Changed name of PDB() to PDBfobj() to avoid symbol collision with PDB
//   proper.
// ****************************************************************************

PDBFileObject *
PDBReader::PDBfobj()
{
    return pdb;
}

// ****************************************************************************
// Method: PDBReader::VariableData::VariableData
//
// Purpose: 
//   Constructor for the PDBReader::VariableData class.
//
// Programmer: Brad Whitlock
// Creation:   Thu Jun 26 14:55:56 PST 2003
//
// Modifications:
//   
// ****************************************************************************

PDBReader::VariableData::VariableData(const std::string &name) : varName(name)
{
    data = 0;
    dataType = NO_TYPE;
    dims = 0;
    nDims = 0;
    nTotalElements = 0;
}

// ****************************************************************************
// Method: PDBReader::VariableData::~VariableData
//
// Purpose: 
//   Destructor for the PDBReader::VariableData class.
//
// Programmer: Brad Whitlock
// Creation:   Thu Jun 26 14:56:14 PST 2003
//
// Modifications:
//   Brad Whitlock, Thu Sep 2 00:08:45 PDT 2004
//   Replaced the code with FreeData.
//
// ****************************************************************************

PDBReader::VariableData::~VariableData()
{
    FreeData();
}

// ****************************************************************************
// Method: PDBReader::VariableData::ReadValues
//
// Purpose: 
//   Reads data into the object.
//
// Arguments:
//   reader : The PDB reader object to use,
//
// Returns:    True if data was successfully read in; false otherwise.
//
// Programmer: Brad Whitlock
// Creation:   Thu Jun 26 14:56:34 PST 2003
//
// Modifications:
//   
// ****************************************************************************

bool
PDBReader::VariableData::ReadValues(PDBFileObject *pdb)
{
    if(data == 0)
    {
        data = pdb->ReadValues((char *)varName.c_str(), &dataType,
                               &nTotalElements, &dims, &nDims);
    }

    return data != 0;
}

// ****************************************************************************
// Method: PDBReader::VariableData::FreeData
//
// Purpose: 
//   Frees the object's data.
//
// Programmer: Brad Whitlock
// Creation:   Thu Sep 2 00:09:57 PDT 2004
//
// Modifications:
//   
// ****************************************************************************

void
PDBReader::VariableData::FreeData()
{
    if(data != 0)
    {
        pdb_free_void_mem(data, dataType);
        data = 0;
        dataType = NO_TYPE;
    }

    delete [] dims;
    dims = 0;
    nDims = 0;
    nTotalElements = 0;
}
