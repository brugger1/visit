// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#include "VisItDataInterfaceRuntime.h"
#include "VisItDataInterfaceRuntimeP.h"

#include "simv2_CSGMesh.h"
#include "simv2_CurveData.h"
#include "simv2_CurvilinearMesh.h"
#include "simv2_DomainBoundaries.h"
#include "simv2_DomainList.h"
#include "simv2_DomainNesting.h"
#include "simv2_MaterialData.h"
#include "simv2_PointMesh.h"
#include "simv2_RectilinearMesh.h"
#include "simv2_SimulationMetaData.h"
#include "simv2_SpeciesData.h"
#include "simv2_UnstructuredMesh.h"
#include "simv2_VariableData.h"

#include <stdlib.h>
#include <string>

#include <ImproperUseException.h>

#define ALLOC(T) (T*)calloc(1, sizeof(T))
#define FREE(PTR) if(PTR != NULL) free(PTR)

//
// Define a structure that we can use to contain all of the callback function
// pointers that the user supplied.
//
typedef struct
{
    /* Reader functions */
    int  (*cb_ActivateTimestep)(void *);
    void  *cbdata_ActivateTimestep;

    visit_handle  (*cb_GetMetaData)(void *);
    void  *cbdata_GetMetaData;

    visit_handle (*cb_GetMesh)(int, const char *, void *);
    void  *cbdata_GetMesh;

    visit_handle  (*cb_GetMaterial)(int, const char *, void *);
    void  *cbdata_GetMaterial;

    visit_handle  (*cb_GetSpecies)(int, const char *, void *);
    void  *cbdata_GetSpecies;

    visit_handle  (*cb_GetVariable)(int, const char *, void *);
    void  *cbdata_GetVariable;

    visit_handle  (*cb_GetMixedVariable)(int, const char *, void *);
    void  *cbdata_GetMixedVariable;

    visit_handle  (*cb_GetCurve)(const char *, void *);
    void  *cbdata_GetCurve;

    visit_handle  (*cb_GetDomainList)(const char *, void *);
    void  *cbdata_GetDomainList;

    visit_handle  (*cb_GetDomainBoundaries)(const char *, void *);
    void  *cbdata_GetDomainBoundaries;

    visit_handle  (*cb_GetDomainNesting)(const char *, void *);
    void  *cbdata_GetDomainNesting;

    /* Writer functions */
    int (*cb_WriteBegin)(const char *, void *);
    void *cbdata_WriteBegin;

    int (*cb_WriteEnd)(const char *, void *);
    void *cbdata_WriteEnd;

    int (*cb_WriteMesh)(const char *, int, int, visit_handle, visit_handle, void *);
    void *cbdata_WriteMesh;

    int (*cb_WriteVariable)(const char *, const char *, int, visit_handle, visit_handle, void *);
    void *cbdata_WriteVariable;

} data_callback_t;

static data_callback_t *visit_data_callbacks = NULL;

static data_callback_t *GetDataCallbacks()
{
    if(visit_data_callbacks == NULL)
        visit_data_callbacks = ALLOC(data_callback_t);
    return visit_data_callbacks;
}

void
DataCallbacksCleanup(void)
{
    if(visit_data_callbacks != NULL)
    {
        free(visit_data_callbacks);
        visit_data_callbacks = NULL;
    }
}

// *****************************************************************************
// Data Interface functions (Callable from libsimV2)
// *****************************************************************************

void
simv2_set_ActivateTimestep(int (*cb) (void *), void *cbdata)
{
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL)
    {
        callbacks->cb_ActivateTimestep = cb;
        callbacks->cbdata_ActivateTimestep = cbdata;
    }
}

void
simv2_set_GetMetaData(visit_handle (*cb) (void *), void *cbdata)
{
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL)
    {
        callbacks->cb_GetMetaData = cb;
        callbacks->cbdata_GetMetaData = cbdata;
    }
}

void
simv2_set_GetMesh(visit_handle (*cb) (int, const char *, void *), void *cbdata)
{
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL)
    {
        callbacks->cb_GetMesh = cb;
        callbacks->cbdata_GetMesh = cbdata;
    }
}

void
simv2_set_GetMaterial(visit_handle (*cb) (int, const char *, void *), void *cbdata)
{
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL)
    {
        callbacks->cb_GetMaterial = cb;
        callbacks->cbdata_GetMaterial = cbdata;
    }
}

void
simv2_set_GetSpecies(visit_handle (*cb) (int, const char *, void *), void *cbdata)
{
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL)
    {
        callbacks->cb_GetSpecies = cb;
        callbacks->cbdata_GetSpecies = cbdata;
    }
}

void
simv2_set_GetVariable(visit_handle (*cb) (int, const char *, void *), void *cbdata)
{
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL)
    {
        callbacks->cb_GetVariable = cb;
        callbacks->cbdata_GetVariable = cbdata;
    }
}

void
simv2_set_GetMixedVariable(visit_handle (*cb) (int, const char *, void *), void *cbdata)
{
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL)
    {
        callbacks->cb_GetMixedVariable = cb;
        callbacks->cbdata_GetMixedVariable = cbdata;
    }
}

void
simv2_set_GetCurve(visit_handle (*cb) (const char *, void *), void *cbdata)
{
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL)
    {
        callbacks->cb_GetCurve = cb;
        callbacks->cbdata_GetCurve = cbdata;
    }
}

void
simv2_set_GetDomainList(visit_handle (*cb) (const char *, void *), void *cbdata)
{
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL)
    {
        callbacks->cb_GetDomainList = cb;
        callbacks->cbdata_GetDomainList = cbdata;
    }
}

void
simv2_set_GetDomainBoundaries(visit_handle (*cb) (const char *, void *), void *cbdata)
{
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL)
    {
        callbacks->cb_GetDomainBoundaries = cb;
        callbacks->cbdata_GetDomainBoundaries = cbdata;
    }
}

void
simv2_set_GetDomainNesting(visit_handle (*cb) (const char *, void *), void *cbdata)
{
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL)
    {
        callbacks->cb_GetDomainNesting = cb;
        callbacks->cbdata_GetDomainNesting = cbdata;
    }
}

/* Write functions */
void
simv2_set_WriteBegin(int (*cb)(const char *, void *), void *cbdata)
{
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL)
    {
        callbacks->cb_WriteBegin = cb;
        callbacks->cbdata_WriteBegin = cbdata;
    }
}

void
simv2_set_WriteEnd(int (*cb)(const char *, void *), void *cbdata)
{
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL)
    {
        callbacks->cb_WriteEnd = cb;
        callbacks->cbdata_WriteEnd = cbdata;
    }
}

void
simv2_set_WriteMesh(int (*cb)(const char *, int, int, visit_handle, visit_handle, void *), void *cbdata)
{
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL)
    {
        callbacks->cb_WriteMesh = cb;
        callbacks->cbdata_WriteMesh = cbdata;
    }
}

void
simv2_set_WriteVariable(int (*cb)(const char *, const char *, int, visit_handle, visit_handle, void *), void *cbdata)
{
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL)
    {
        callbacks->cb_WriteVariable = cb;
        callbacks->cbdata_WriteVariable = cbdata;
    }
}

// *****************************************************************************
// Callable from SimV2 reader
// *****************************************************************************

// ****************************************************************************
// Method: simv2_invoke_ActivateTimeStep
//
// Purpose: 
//   This function invokes the simulation's ActivateTimeStep callback function.
//
// Programmer: Brad Whitlock
// Creation:   Mon Mar 16 16:12:44 PDT 2009
//
// Modifications:
//   
// ****************************************************************************

int
simv2_invoke_ActivateTimestep(void)
{
    int retval = VISIT_OKAY;
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL && callbacks->cb_ActivateTimestep != NULL)
    {
        retval = (*callbacks->cb_ActivateTimestep)(callbacks->cbdata_ActivateTimestep);
    }
    return retval;
}

// ****************************************************************************
// Method: simv2_invoke_GetMetaData
//
// Purpose: 
//   This function invokes the simulation's GetMetaData callback function and
//   returns a VisIt_SimulationMetaData structure that was populated by the
//   simulation.
//
// Arguments:
//
// Returns:    
//
// Note:       This encapsulates the code to create an object ourselves,
//             pass it to the sim along with any user callback data that was
//             provided, and check the results for errors that could threaten
//             VisIt's operation.
//
// Programmer: Brad Whitlock
// Creation:   Thu Feb 26 16:03:44 PST 2009
//
// Modifications:
//   
// ****************************************************************************

visit_handle
simv2_invoke_GetMetaData(void)
{
    visit_handle h = VISIT_INVALID_HANDLE;
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL && callbacks->cb_GetMetaData != NULL)
    {
        h = (*callbacks->cb_GetMetaData)(callbacks->cbdata_GetMetaData);

        if(h != VISIT_INVALID_HANDLE)
        {
            if(simv2_ObjectType(h) != VISIT_SIMULATION_METADATA)
            {
                simv2_FreeObject(h);
                EXCEPTION1(ImproperUseException, 
                    "The simulation returned a handle for an object other "
                    "than simulation metadata."
                );
            }

            if(simv2_SimulationMetaData_check(h) == VISIT_ERROR)
            {
                simv2_FreeObject(h);
                EXCEPTION1(ImproperUseException, 
                    "The metadata returned by the simulation did not pass "
                    "a consistency check.");
            }
        }
    }
    return h;
}

visit_handle
simv2_invoke_GetMesh(int dom, const char *name)
{
    visit_handle h = VISIT_INVALID_HANDLE;
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL && callbacks->cb_GetMesh != NULL)
    {
        h = (*callbacks->cb_GetMesh)(dom, name, callbacks->cbdata_GetMesh);

        if(h != VISIT_INVALID_HANDLE)
        {
            int msgno = 0, err = VISIT_ERROR;
            int objType = simv2_ObjectType(h);
            switch (objType)
            {
            case VISIT_CSG_MESH:
                err = simv2_CSGMesh_check(h);
                break;
            case VISIT_CURVILINEAR_MESH:
                err = simv2_CurvilinearMesh_check(h);
                break;
            case VISIT_RECTILINEAR_MESH:
                err = simv2_RectilinearMesh_check(h);
                break;
            case VISIT_POINT_MESH:
                err = simv2_PointMesh_check(h);
                break;
            case VISIT_UNSTRUCTURED_MESH:
                err = simv2_UnstructuredMesh_check(h);
                break;
            default:
                msgno = 1;
            }

            if(err == VISIT_ERROR)
            {
                simv2_FreeObject(h);
                EXCEPTION1(ImproperUseException, 
                    (msgno == 0) ? 
                    "The mesh returned by the simulation did not pass "
                    "a consistency check."
                    : 
                    "The simulation returned a handle for an object other "
                    "than a mesh."
                );
            }
        }
    }
    return h;
}

visit_handle
simv2_invoke_GetMaterial(int dom, const char *name)
{
    visit_handle h = VISIT_INVALID_HANDLE;
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL && callbacks->cb_GetMaterial != NULL)
    {
        h = (*callbacks->cb_GetMaterial)(dom, name, callbacks->cbdata_GetMaterial);

        if(h != VISIT_INVALID_HANDLE)
        {
            if(simv2_ObjectType(h) != VISIT_MATERIAL_DATA)
            {
                simv2_FreeObject(h);
                EXCEPTION1(ImproperUseException, 
                    "The simulation returned a handle for an object other "
                    "than a material data object."
                );
            }

            if(simv2_MaterialData_check(h) == VISIT_ERROR)
            {
                simv2_FreeObject(h);
                EXCEPTION1(ImproperUseException, 
                    "The material returned by the simulation did not pass "
                    "a consistency check.");
            }
        }
    }
    return h;
}

visit_handle
simv2_invoke_GetSpecies(int dom, const char *name)
{
    visit_handle h = VISIT_INVALID_HANDLE;
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL && callbacks->cb_GetSpecies != NULL)
    {
        h = (*callbacks->cb_GetSpecies)(dom, name, callbacks->cbdata_GetSpecies);

        if(h != VISIT_INVALID_HANDLE)
        {
            if(simv2_ObjectType(h) != VISIT_SPECIES_DATA)
            {
                simv2_FreeObject(h);
                EXCEPTION1(ImproperUseException, 
                    "The simulation returned a handle for an object other "
                    "than a species data object."
                );
            }

            if(simv2_SpeciesData_check(h) == VISIT_ERROR)
            {
                simv2_FreeObject(h);
                EXCEPTION1(ImproperUseException, 
                    "The species returned by the simulation did not pass "
                    "a consistency check.");
            }
        }
    }
    return h;
}

visit_handle
simv2_invoke_GetVariable(int dom, const char *name)
{
    visit_handle h = VISIT_INVALID_HANDLE;
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL && callbacks->cb_GetVariable != NULL)
    {
        h = (*callbacks->cb_GetVariable)(dom, name, callbacks->cbdata_GetVariable);

        if(h != VISIT_INVALID_HANDLE)
        {
            if(simv2_ObjectType(h) != VISIT_VARIABLE_DATA)
            {
                simv2_FreeObject(h);
                EXCEPTION1(ImproperUseException, 
                    "The simulation returned a handle for an object other "
                    "than a variable."
                );
            }
        }
    }
    return h;
}

visit_handle
simv2_invoke_GetMixedVariable(int dom, const char *name)
{
    visit_handle h = VISIT_INVALID_HANDLE;
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL && callbacks->cb_GetMixedVariable != NULL)
    {
        h = (*callbacks->cb_GetMixedVariable)(dom, name, callbacks->cbdata_GetMixedVariable);

        if(h != VISIT_INVALID_HANDLE)
        {
            if(simv2_ObjectType(h) != VISIT_VARIABLE_DATA)
            {
                simv2_FreeObject(h);
                EXCEPTION1(ImproperUseException, 
                    "The simulation returned a handle for an object other "
                    "than a variable."
                );
            }
        }
    }
    return h;
}

visit_handle
simv2_invoke_GetCurve(const char *name)
{
    visit_handle h = VISIT_INVALID_HANDLE;
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL && callbacks->cb_GetCurve != NULL)
    {
        h = (*callbacks->cb_GetCurve)(name, callbacks->cbdata_GetCurve);

        if(h != VISIT_INVALID_HANDLE)
        {
            if(simv2_ObjectType(h) != VISIT_CURVE_DATA)
            {
                simv2_FreeObject(h);
                EXCEPTION1(ImproperUseException, 
                    "The simulation returned a handle for an object other "
                    "than a curve."
                );
            }

            if(simv2_CurveData_check(h) == VISIT_ERROR)
            {
                simv2_FreeObject(h);
                EXCEPTION1(ImproperUseException, 
                    "The curve returned by the simulation did not pass "
                    "a consistency check."
                );
            }
        }
    }
    return h;
}

visit_handle
simv2_invoke_GetDomainList(const char *name)
{
    visit_handle h = VISIT_INVALID_HANDLE;
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL && callbacks->cb_GetDomainList != NULL)
    {
        h = (*callbacks->cb_GetDomainList)(name, callbacks->cbdata_GetDomainList);

        if(h != VISIT_INVALID_HANDLE)
        {
            if(simv2_ObjectType(h) != VISIT_DOMAINLIST)
            {
                simv2_FreeObject(h);
                EXCEPTION1(ImproperUseException, 
                    "The simulation returned a handle for an object other "
                    "than a domain list."
                );
            }

            if(simv2_DomainList_check(h) == VISIT_ERROR)
            {
                simv2_FreeObject(h);
                EXCEPTION1(ImproperUseException, 
                    "The domain list returned by the simulation did not pass "
                    "a consistency check."
                );
            }
        }
    }
    return h;
}

visit_handle 
simv2_invoke_GetDomainBoundaries(const char *name)
{
    visit_handle h = VISIT_INVALID_HANDLE;
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL && callbacks->cb_GetDomainBoundaries != NULL)
    {
        h = (*callbacks->cb_GetDomainBoundaries)(name, callbacks->cbdata_GetDomainBoundaries);
    }
    return h;
}

visit_handle 
simv2_invoke_GetDomainNesting(const char *name)
{
    visit_handle h = VISIT_INVALID_HANDLE;
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL && callbacks->cb_GetDomainNesting != NULL)
    {
        h = (*callbacks->cb_GetDomainNesting)(name, callbacks->cbdata_GetDomainNesting);
    }
    return h;
}

/* Writer functions. */

int
simv2_invoke_WriteBegin(const char *name)
{
    int ret = VISIT_ERROR;
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL && callbacks->cb_WriteBegin != NULL)
        ret = (*callbacks->cb_WriteBegin)(name, callbacks->cbdata_WriteBegin);
    return ret;
}

int
simv2_invoke_WriteEnd(const char *name)
{
    int ret = VISIT_ERROR;
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL && callbacks->cb_WriteEnd != NULL)
        ret = (*callbacks->cb_WriteEnd)(name, callbacks->cbdata_WriteEnd);
    return ret;
}

int
simv2_invoke_WriteMesh(const char *name, int chunk, int meshType, visit_handle md, visit_handle mmd)
{
    int ret = VISIT_ERROR;
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL && callbacks->cb_WriteMesh != NULL)
        ret = (*callbacks->cb_WriteMesh)(name, chunk, meshType, md, mmd, callbacks->cbdata_WriteMesh);
    return ret;
}

int
simv2_invoke_WriteVariable(const char *name, const char *arrName, int chunk,
    visit_handle data, visit_handle smd)
{
    int ret = VISIT_ERROR;
    data_callback_t *callbacks = GetDataCallbacks();
    if(callbacks != NULL && callbacks->cb_WriteVariable != NULL)
        ret = (*callbacks->cb_WriteVariable)(name, arrName, chunk, data, smd, callbacks->cbdata_WriteVariable);
    return ret;
}

int
simv2_ObjectType(visit_handle h)
{
    VisIt_ObjectBase *obj = VisItGetPointer(h);
    int t = -1;
    if(obj != NULL)
        t = obj->objectType();
    return t;
}

int
simv2_FreeObject(visit_handle h)
{
    int retval = VISIT_ERROR;
    VisIt_ObjectBase *obj = VisItGetPointer(h);
    if(obj != NULL)
    {
        // Rely on the virtual destructor
        delete obj;
        VisItFreePointer(h);
        retval = VISIT_OKAY;
    }
    return retval;
}
