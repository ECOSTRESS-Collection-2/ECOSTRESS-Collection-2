// metadata.c
// See metadata.h for documentation and credits.
// 10-JUN-2019 R. Freepartner Provide API to suppress _FillValue field.
// @copyright (c) 2019 JPL, All rights reserved
#include "metadata.h"
#include "bool.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUF_SIZE 4096

/// @brief Static data tables to convert HDF datatype to 
/// size in bytes
static int hdf4_type[12] = 
{
    DFNT_CHAR8, DFNT_UCHAR8, 
    DFNT_FLOAT32, DFNT_FLOAT64, 
    DFNT_INT8, DFNT_UINT8, DFNT_INT16, DFNT_UINT16,
    DFNT_INT32, DFNT_UINT32, DFNT_INT64, DFNT_UINT64
};

/// @brief Sizes for each datatype
static int hdf4_size[12] = 
{
    SIZE_NUCHAR8, SIZE_NCHAR8,
    SIZE_NFLOAT32, SIZE_NFLOAT64,
    SIZE_NINT8, SIZE_NUINT8, SIZE_NINT16, SIZE_NUINT16,
    SIZE_NINT32, SIZE_NUINT32, SIZE_NINT64, SIZE_NUINT64
};

/// @brief find size for type.
/// @param type hdf4 type code
/// @return -1=type not found, else size for type
static int get_hdf4_size(int type)
{
    int i;
    for (i = 0; i < 12; i++)
    {
        if (type == hdf4_type[i])
            return hdf4_size[i];
    }
    return -1;
}

/// @brief find hdf5 type for hdf4 type.
/// @param type hdf4 type code
/// @return -1=type not found, else size for type
static hid_t get_hdf5_type(int h4type)
{
    if (h4type == DFNT_CHAR8) return H5T_C_S1;
    if (h4type == DFNT_UCHAR8) return H5T_C_S1;
    if (h4type == DFNT_FLOAT32) return H5T_NATIVE_FLOAT;
    if (h4type == DFNT_FLOAT64) return H5T_NATIVE_DOUBLE;
    if (h4type == DFNT_INT8) return H5T_NATIVE_SCHAR;
    if (h4type == DFNT_UINT8) return H5T_NATIVE_UCHAR;
    if (h4type == DFNT_INT16) return H5T_NATIVE_SHORT;
    if (h4type == DFNT_UINT16) return H5T_NATIVE_USHORT;
    if (h4type == DFNT_INT32) return H5T_NATIVE_INT;
    if (h4type == DFNT_UINT32) return H5T_NATIVE_UINT;
    if (h4type == DFNT_INT64) return H5T_NATIVE_LONG;
    if (h4type == DFNT_UINT64) return H5T_NATIVE_ULONG;
    return -1;
}

/* not used
static int get_hdf4_type(hid_t hdf5_native)
{
    if (hdf5_native == H5T_NATIVE_CHAR) return DFNT_CHAR8;
    if (hdf5_native == H5T_NATIVE_UCHAR) return DFNT_UINT8;
    if (hdf5_native == H5T_NATIVE_FLOAT) return DFNT_FLOAT32;
    if (hdf5_native == H5T_NATIVE_DOUBLE) return DFNT_FLOAT64;
    if (hdf5_native == H5T_NATIVE_SCHAR) return DFNT_INT8;
    if (hdf5_native == H5T_NATIVE_SHORT) return DFNT_INT16;
    if (hdf5_native == H5T_NATIVE_USHORT) return DFNT_UINT16;
    if (hdf5_native == H5T_NATIVE_INT) return DFNT_INT32;
    if (hdf5_native == H5T_NATIVE_UINT) return DFNT_UINT32;
    if (hdf5_native == H5T_NATIVE_LONG) return DFNT_INT64;
    if (hdf5_native == H5T_NATIVE_ULONG) return DFNT_UINT64;
    return -1;
}
*/
void constructMetadataAttribute(MetadataAttribute *att)
{
    assert(att != NULL);
    att->index = 0;
    att->name[0] = '\0';
    att->type = 0;
    att->count = 0;
    att->value = NULL;
    att->next = NULL;
}

void destructMetadataAttribute(MetadataAttribute *att)
{
    assert(att != NULL);
    if (att->value != NULL) free(att->value);
    att->value = NULL;
    if (att->next != NULL) 
    {
        destructMetadataAttribute(att->next);
        free(att->next);
    }
    att->next = NULL;
}

void setValue(MetadataAttribute *att, int in_type, int in_count, 
                char *in_value)
{
    if (att->value) free(att->value);
    att->type = in_type;
    att->count = in_count;
    att->value = in_value;
}

void constructMetadataSet(MetadataSet *set, int32 in_objid)
{
    assert(set != NULL);
    set->objid = in_objid;
    set->natts = 0;
    set->head = NULL;
    set->tail = NULL;
}

void destructMetadataSet(MetadataSet *set)
{
    assert(set != NULL);
    set->natts = 0;
    if (set->head)
    {
        destructMetadataAttribute(set->head);
        set->head = NULL;
    }
}

char *dupuint8(uint8 ivalue)
{
    uint8 *p = (uint8*)malloc(sizeof(uint8));
    *p = ivalue;
    return (char*)p;
}

char *dupint32(int32 ivalue)
{
    int32 *p = (int32*)malloc(sizeof(int32));
    *p = ivalue;
    return (char*)p;
}

char *dupfloat32(float32 dvalue)
{
    float32 *p = (float32*)malloc(sizeof(float32));
    *p = dvalue;
    return (char*)p;
}

char *dupfloat64(double dvalue)
{
    double *p = (double*)malloc(sizeof(double));
    *p = dvalue;
    return (char*)p;
}

int loadMetadataHdf4(MetadataSet *set, int32 obj_id)
{
    // Start with known state.
    constructMetadataSet(set, obj_id);

    // For each attribute, create an attribute object and load it with the 
    // data for the attribute
    MetadataAttribute *att;
    int32 attr_index = 0;
    char attr_name[METADATA_NAME_SIZE];
    int32 data_type;
    int32 attr_count;
    size_t value_size;
    int nil_size = 0;
    int32 stat = 0;
    char *findillch;

    while (stat >= 0)
    {
        stat = SDattrinfo(obj_id, attr_index, attr_name, &data_type, &attr_count);
        if (stat == FAIL) break; // No more attributes

        // Create a new attribute object.
        att = (MetadataAttribute*)malloc(sizeof(MetadataAttribute));
        assert(att != NULL);
        constructMetadataAttribute(att);
        att->index = attr_index;
        strcpy(att->name, attr_name);
        // NetCDF and HDF5 do not like having / in names.
        findillch = strstr(att->name, "/");
        while (findillch != NULL)
        {
            *findillch = '-';
            findillch = strstr(att->name, "/");
        }
        att->type = data_type;
        att->count = attr_count;

        // Link new attribute to the end of the list for the set.
        if (set->tail == NULL)
        {
            set->head = att; // Add first attribute ot the head of the list.
        }
        else
        {
            set->tail->next = att; // Add additional attribute to the tail.
        }
        set->tail = att;
        set->natts++;

        // Allocate the buffer to contain the value.
        value_size = get_hdf4_size(att->type) * att->count;
        if (att->type == DFNT_CHAR8 || att->type == DFNT_UCHAR8)
        {
            // For char strings, add a byte for the NIL.
            nil_size = 1;
        }
        else
        {
            nil_size = 0;
        }
        att->value = (char*)malloc(value_size + nil_size);
        assert(att->value != NULL);

        // Read the attribute value into the buffer.
        stat = SDreadattr(obj_id, attr_index, att->value);
        if (stat == FAIL)
        {
            fprintf(stderr, "SDreadattr failed for obj_id=%d attr_index=%d "
                    "at %s:%d", obj_id, attr_index, __FILE__, __LINE__);
            destructMetadataSet(set);
            return -1;
        }
        if (att->type == DFNT_CHAR8 || att->type == DFNT_UCHAR8)
        {
            // Add NIL terminator to strings,
            att->value[att->count] = '\0';
        }

        // Get the next attribute index.
        attr_index++;
    }

    // Successful completion
    return 0;
}

int loadMetadataHdf5(MetadataSet *set, hid_t obj_id)
{
    static char read_buf[BUF_SIZE];
    H5O_info_t obj_info;
    // Start with known state.
    constructMetadataSet(set, obj_id);

    // For each attribute, create an attribute object and load it with the 
    // data for the attribute
    MetadataAttribute *att;
    hid_t attr;
    hsize_t attr_index = 0;
    char attr_name[METADATA_NAME_SIZE];
    hsize_t ndims;
    hsize_t dims[4];
    hsize_t maxdims[4];
    hsize_t data_count;
    hid_t data_space;
    hid_t data_type;
    hid_t dclass;
    size_t value_size;
    herr_t stat = 0;
    int readstat = 0;
    // Get number of attributes from info.
    stat = H5Oget_info(obj_id, &obj_info);
    if (stat < 0)
    {
        fprintf(stderr, "H5Oget_info fails to get info for h5 file metadata.\n");
        return -1;
    }

    // Loop for each attribute
    for (attr_index = 0; attr_index < obj_info.num_attrs; attr_index++)
    {
        attr = H5Aopen_by_idx(obj_id, ".", H5_INDEX_CRT_ORDER, H5_ITER_INC, attr_index, 
                H5P_DEFAULT, H5P_DEFAULT);
        if (attr == -1)
        {
            fprintf(stderr, "H5Aopen_by_idx failed for attribute #%llu\n", attr_index);
            continue;
        }
        H5Aget_name(attr, METADATA_NAME_SIZE-1, attr_name);

        // Obtain the data type information for the atttribute.
        data_type = H5Aget_type(attr);
        value_size = H5Tget_size(data_type);
        dclass = H5Tget_class(data_type);
        H5Tclose(data_type);

        // Use data space to get array size
        data_space = H5Aget_space(attr);
        ndims = H5Sget_simple_extent_ndims(data_space);
        assert(ndims < 5); // Only have room for 4
        if (ndims > 0)
        {
            H5Sget_simple_extent_dims(data_space, dims, maxdims);
            data_count = dims[0];
        }
        else
            data_count = 1;
        H5Sclose(data_space);

        // Create a new attribute object: att = new MetadataAttribute();
        att = (MetadataAttribute*)malloc(sizeof(MetadataAttribute));
        assert(att != NULL);
        constructMetadataAttribute(att);
        att->index = attr_index;
        strcpy(att->name, attr_name);

        if (dclass == H5T_INTEGER)
        {
            att->type = DFNT_INT32;
            value_size = sizeof(int32);
        }
        else if (dclass == H5T_FLOAT)
        {
            att->type = DFNT_FLOAT64;
            value_size = sizeof(float64);
        }
        else
        {
            att->type = DFNT_CHAR8;
        }
        att->count = data_count;

        // Link new attribute to the end of the list for the set.
        if (set->tail == NULL)
        {
            set->head = att; // Add first attribute ot the head of the list.
        }
        else
        {
            set->tail->next = att; // Add additional attribute to the tail.
        }
        set->tail = att;
        set->natts++;

        data_type = H5Tcopy(get_hdf5_type(att->type));
        // Allocate the buffer to contain the value.
        if (att->type == DFNT_CHAR8 || att->type == DFNT_UCHAR8)
        {
            // For char strings, add a byte for the NIL.
            H5Tset_size(data_type, value_size);
            att->count = value_size;
        }

        // Must fit into read_buf
        assert(value_size * data_count < BUF_SIZE);

        // Read the attribute value into the buffer.
        stat = H5Aread(attr, data_type, read_buf);
        if (stat == FAIL)
        {
            fprintf(stderr, "Warning: H5Aread failed for attr_index=%llu "
                    "name=%s type=%d size=%lu at %s:%d\n", 
                    attr_index, att->name, att->type, value_size,
                    __FILE__, __LINE__);
            readstat = -1;
        }

        if (att->type == DFNT_CHAR8 || att->type == DFNT_UCHAR8)
        {
            // Add NIL terminator to strings,
            read_buf[value_size] = '\0';
            att->value = strdup(read_buf);
        }
        else if (value_size == 4)
        {
            att->value = (char*)malloc(data_count * value_size);
            memcpy(att->value, read_buf, data_count * value_size);
        }
        else if (value_size == 8)
        {
            att->value = (char*)malloc(data_count * value_size);
            memcpy(att->value, read_buf, data_count * value_size);
        }
        else
        {
            fprintf(stderr, "Attribute %s has unknown type/size=%d/%lu\n",
                    att->name, att->type, value_size);
            readstat = -1;
        }

        H5Tclose(data_type);
        H5Aclose(attr);
    }
    return readstat;
}

int loadMetadataGroup(MetadataSet *set, hid_t hfileid, const char *groupname)
{
    // Open the existing group in the file.
    hid_t gid = H5Gopen2(hfileid, groupname, H5P_DEFAULT);
    if (gid < 0)
        return -1;

    MetadataAttribute *att;
    hsize_t  num_obj = 0;
    char objname[256];
    size_t objsize = 255;
    ssize_t objlen;
    hid_t dsetid;
    hid_t dtypeid;
    hid_t dtype;
    hid_t strtype;
    size_t ssize;
    char *valchar[1];

    herr_t stat = H5Gget_num_objs(gid, &num_obj);
    if (stat < 0)
    {
        H5Gclose(gid);
        return -1;
    }

    hsize_t idx;
    for (idx = 0; idx < num_obj; idx++)
    {
        objlen = H5Gget_objname_by_idx(gid, idx, objname, objsize);
        assert(objlen > 0);
        dsetid = H5Dopen2(gid, objname, H5P_DEFAULT);
        if (dsetid == -1)
        {
            H5Gclose(gid);
            return -1;
        }
        dtypeid = H5Dget_type(dsetid);
        dtype = H5Tget_class(dtypeid);
        ssize = H5Dget_storage_size(dsetid);

        // Prepare the new attribute.
        att = (MetadataAttribute*)malloc(sizeof(MetadataAttribute));
        assert(att != NULL);
        constructMetadataAttribute(att);
        att->index = idx;
        strcpy(att->name, objname);
        att->count = 1;

        if (dtype == H5T_INTEGER)
        {
            att->type = DFNT_INT32;
            att->value = (char*)malloc(ssize);
            assert(att->value != NULL);
            stat = H5Dread(dsetid, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    att->value);
            if (stat < 0)
            {
                H5Gclose(gid);
                return -1;
            }
        }
        else if (dtype == H5T_FLOAT)
        {
            att->type = DFNT_FLOAT64;
            att->value = (char*)malloc(ssize);
            assert(att->value != NULL);
            stat = H5Dread(dsetid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    att->value);
            if (stat < 0)
            {
                H5Gclose(gid);
                return -1;
            }
        }
        else if (dtype == H5T_STRING)
        {
            strtype = H5Tcopy(dtypeid);
            stat = H5Dread(dsetid, strtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, valchar);
            if (stat < 0)
            {
                H5Gclose(gid);
                return -1;
            }
            att->type = DFNT_CHAR8;
            att->value = strdup(valchar[0]);
            att->count = strlen(att->value);
            assert(att->value != NULL);
        }
        else //ignore dataset of unknown type
        {
            free(att);
            att = NULL;
        }
        if (att != NULL)
        {
            // Link new attribute to the end of the list for the set.
            if (set->tail == NULL)
            {
                set->head = att; // Add first attribute ot the head of the list.
            }
            else
            {
                set->tail->next = att; // Add additional attribute to the tail.
            }
            set->tail = att;
            set->natts++;

        }
        H5Tclose(dtypeid);
        H5Dclose(dsetid);
    }
    H5Gclose(gid);
    return 0;
}

MetadataAttribute *findAttribute(MetadataSet *set, const char *attrname)
{
    // This will be a simple O(N) list traversal.
    MetadataAttribute *att = set->head;
    while (att)
    {
        if (strcmp(attrname, att->name) == 0)
            return att;
        att = att->next;
    }
    return NULL; // not found
}

void addAttribute(MetadataSet *set, MetadataAttribute *att)
{
    att->next = NULL; // Adding this at the end.
    // Check for first attribute to add as head.
    if (set->head == NULL)
    {
        set->head = att;
        att->index = 0;
    }
    else
    {
        att->index = set->natts;
        set->tail->next = att;
    }
    set->natts++;
    set->tail = att;
}

void deleteAttribute(MetadataAttribute *att)
{
    att->index = -1;
}

int updateValue(MetadataSet *set, const char *attr_name, int in_type, 
                int in_count, char *in_value)
{
    MetadataAttribute *att = findAttribute(set, attr_name);
    if (att == NULL)
    {
        att = (MetadataAttribute*)malloc(sizeof(MetadataAttribute));
        constructMetadataAttribute(att);
        strcpy(att->name, attr_name);
        att->type = in_type;
        att->count = in_count;
        att->value = in_value;
        addAttribute(set, att);
        return 1;
    }
    setValue(att, in_type, in_count, in_value);
    return 0;

    // NOTE: So far there are no error conditions.
}

int writeMetadataHdf5(MetadataSet *set, hid_t loc_id)
{
    herr_t stat;
    hid_t att_id;
    hid_t space_id;
    hsize_t dims;
    hid_t type_id;

    // Loop through the list of attributes and write each value.
    MetadataAttribute *att = set->head;
    while (att != NULL)
    {
        if (att->index < 0)
        {
            // Deleted, skip it.
            att = att->next;
            continue;
        }

        // Open the attribute if it exists.
        if (H5Aexists_by_name(loc_id, ".", att->name, H5P_DEFAULT))
        {
            att_id = H5Aopen_by_name(loc_id, ".", att->name, 
                     H5P_DEFAULT, H5P_DEFAULT);
            if (att_id < 0)
            {
                fprintf(stderr, "H5Aopen_by_name failed for %s on loc_id=%d "
                        "at %s:%d\n",
                        att->name, loc_id, __FILE__, __LINE__);
                return -1;
            }
        }
        else
        {
            // This attribute doesn't exist on loc_id, so create it.
            if (att->type != DFNT_CHAR8 && att->type != DFNT_UCHAR8)
            {
                // For non-strings, create a type based on the HDF5 type.
                type_id = H5Tcopy(get_hdf5_type(att->type));
                // For non-strings, create a simple space including the dimension.
                dims = att->count;
                space_id = H5Screate_simple(1, &dims, NULL);
            }
            else
            {
                // For strings, create a type using strlen for string size.
                type_id = H5Tcopy(H5T_C_S1);
                H5Tset_size(type_id, att->count);
                H5Tset_strpad(type_id, H5T_STR_NULLTERM);
                // Strings use Scalar space.
                space_id = H5Screate(H5S_SCALAR);
            }

            att_id = H5Acreate_by_name(loc_id, ".", att->name, type_id, space_id, 
                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Sclose(space_id);
            H5Tclose(type_id);
            if (att_id < 0)
            {
                fprintf(stderr, "H5Acreate_by_name failed for %s on loc_id=%d "
                        "type_id=%d space_id=%d at %s:%d\n",
                        att->name, loc_id, type_id, space_id, 
                        __FILE__, __LINE__);
                return -1;
            }
        }

        // Write the value to the attribute.
        type_id = H5Aget_type(att_id);
        stat = H5Awrite(att_id, type_id, (const void *)att->value);
        H5Tclose(type_id);
        H5Aclose(att_id);
        if (stat < 0)
        {
            fprintf(stderr, "H5Awrite failed for %s on loc_id=%d "
                    "att_id=%d at %s:%d\n",
                    att->name, loc_id, att_id, __FILE__, __LINE__);
            return -1;
        }

        // Move to the next attribute
        att = att->next;
    }
    
    // Success
    return 0;
}

int writeMetadataGroup(MetadataSet *set, hid_t hfileid, const char *groupname)
{
    herr_t stat;
    hid_t att_id;
    hsize_t dims;
    hid_t space_id;
    hid_t type_id;
    char *writefrom;

    // Attempt to open the group in the HDF5 file.
    hid_t gid = H5Gopen2(hfileid, groupname, H5P_DEFAULT);
    if (gid < 0)
        return -1;

    // Loop through the list of attributes and write each value.
    MetadataAttribute *att = set->head;
    while (att != NULL)
    {
        if (att->index < 0)
        {
            // Deleted, skip it.
            att = att->next;
            continue;
        }

        space_id = H5Screate(H5S_SCALAR);
        if (att->type != DFNT_CHAR8 && att->type != DFNT_UCHAR8)
        {
            type_id = H5Tcopy(get_hdf5_type(att->type));
            // For non-strings, create a simple space including the dimension.
            dims = att->count;
            space_id = H5Screate_simple(1, &dims, NULL);
            writefrom = att->value; // write from the location of the data
        }
        else
        {
            // For strings, create an unlimited variable length string dataset.
            type_id = H5Tcopy(H5T_C_S1);
            H5Tset_size(type_id, H5T_VARIABLE);
            H5Tset_strpad(type_id, H5T_STR_NULLTERM);
            // Strings use Scalar space.
            space_id = H5Screate(H5S_SCALAR);
            writefrom = (char*)&att->value; // write from the location of the pointer
        }

        // Create the dataset in the group
        att_id = H5Dcreate2(gid, att->name, type_id, space_id, 
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (att_id < 0)
        {
            H5Gclose(gid);
            H5Sclose(space_id);
            H5Tclose(type_id);
            return -1;
        }

        // Write the value to the dataset.
        stat = H5Dwrite(att_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    writefrom);
        H5Sclose(space_id);
        H5Tclose(type_id);
        H5Dclose(att_id);
        if (stat < 0)
        {
            H5Gclose(gid);
            return -1;
        }

        // Move to the next attribute
        att = att->next;
    }
    
    // Success
    H5Gclose(gid);
    return 0;
}

// for ECOSTRESS
// L2_CLOUD
char StructMetadata_CLOUD[] =
"GROUP=SwathStructure\n"
"\tGROUP=SWATH_1\n"
"\t\tSwathName=\"SDS\"\n"
"\t\tGROUP=Dimension\n"
"\t\t\tOBJECT=Dimension_1\n"
"\t\t\t\tDimensionName=\"Along_Track\"\n"
"\t\t\t\tSize=5325\n"
"\t\t\tEND_OBJECT=Dimension_1\n"
"\t\t\tOBJECT=Dimension_2\n"
"\t\t\t\tDimensionName=\"Along_Scan\"\n"
"\t\t\t\tSize=5325\n"
"\t\t\tEND_OBJECT=Dimension_2\n"
"\t\tEND_GROUP=Dimension\n"
"\t\t\tOBJECT=DataField_1\n"
"\t\t\t\tDataFieldName=\"CloudMask\"\n"
"\t\t\t\tDataType=NC_UBYTE\n"
"\t\t\t\tDimList=(\"Along_Track\",\"Along_Scan\")\n"
"\t\t\tEND_OBJECT=DataField_1\n"
"\tEND_GROUP=SWATH_1\n"
"END_GROUP=SwathStructure\n"
"END\n";

char *getCLOUDStructMetadata()
{
    return StructMetadata_CLOUD;
}

// L2_LSTE
char StructMetadata_LSTE[] =
"GROUP=SwathStructure\n"
"\tGROUP=SWATH_1\n"
"\t\tSwathName=\"LSTE\"\n"
"\t\tGROUP=Dimension\n"
"\t\t\tOBJECT=Dimension_1\n"
"\t\t\t\tDimensionName=\"Along_Track\"\n"
"\t\t\t\tSize=####\n"
"\t\t\tEND_OBJECT=Dimension_1\n"
"\t\t\tOBJECT=Dimension_2\n"
"\t\t\t\tDimensionName=\"Along_Scan\"\n"
"\t\t\t\tSize=####\n"
"\t\t\tEND_OBJECT=Dimension_2\n"
"\t\tEND_GROUP=Dimension\n"
"\t\tGROUP=DimensionMap\n"
"\t\tEND_GROUP=DimensionMap\n"
"\t\tGROUP=IndexDimensionMap\n"
"\t\tEND_GROUP=IndexDimensionMap\n"
"\t\tGROUP=GeoField\n"
"\t\t\tOBJECT=GeoField_1\n"
"\t\t\t\tGeoFieldName=\"Latitude\"\n"
"\t\t\t\tDataType=NC_FLOAT\n"
"\t\t\t\tDimList=(\"Along_Track\",\"Along_Scan\")\n"
"\t\t\tEND_OBJECT=GeoField_1\n"
"\t\t\tOBJECT=GeoField_2\n"
"\t\t\t\tGeoFieldName=\"Longitude\"\n"
"\t\t\t\tDataType=NC_FLOAT\n"
"\t\t\t\tDimList=(\"Along_Track\",\"Along_Scan\")\n"
"\t\t\tEND_OBJECT=GeoField_2\n"
"\t\tEND_GROUP=GeoField\n"
"\t\tGROUP=DataField\n"
"\t\t\tOBJECT=DataField_1\n"
"\t\t\t\tDataFieldName=\"Emis_1\"\n"
"\t\t\t\tDataType=NC_UBYTE\n"
"\t\t\t\tDimList=(\"Along_Track\",\"Along_Scan\")\n"
"\t\t\tEND_OBJECT=DataField_1\n"
"\t\t\tOBJECT=DataField_2\n"
"\t\t\t\tDataFieldName=\"Emis_2\"\n"
"\t\t\t\tDataType=NC_USHORT\n"
"\t\t\t\tDimList=(\"Along_Track\",\"Along_Scan\")\n"
"\t\t\tEND_OBJECT=DataField_2\n"
"\t\t\tOBJECT=DataField_3\n"
"\t\t\t\tDataFieldName=\"Emis_3\"\n"
"\t\t\t\tDataType=NC_UBYTE\n"
"\t\t\t\tDimList=(\"Along_Track\",\"Along_Scan\")\n"
"\t\t\tEND_OBJECT=DataField_3\n"
"\t\t\tOBJECT=DataField_4\n"
"\t\t\t\tDataFieldName=\"Emis_4\"\n"
"\t\t\t\tDataType=NC_UBYTE\n"
"\t\t\t\tDimList=(\"Along_Track\",\"Along_Scan\")\n"
"\t\t\tEND_OBJECT=DataField_4\n"
"\t\t\tOBJECT=DataField_5\n"
"\t\t\t\tDataFieldName=\"Emis_5\"\n"
"\t\t\t\tDataType=NC_UBYTE\n"
"\t\t\t\tDimList=(\"Along_Track\",\"Along_Scan\")\n"
"\t\t\tEND_OBJECT=DataField_5\n"
"\t\t\tOBJECT=DataField_6\n"
"\t\t\t\tDataFieldName=\"EmisWB\"\n"
"\t\t\t\tDataType=NC_UBYTE\n"
"\t\t\t\tDimList=(\"Along_Track\",\"Along_Scan\")\n"
"\t\t\tEND_OBJECT=DataField_6\n"
"\t\t\tOBJECT=DataField_7\n"
"\t\t\t\tDataFieldName=\"Emis_1_err\"\n"
"\t\t\t\tDataType=NC_UBYTE\n"
"\t\t\t\tDimList=(\"Along_Track\",\"Along_Scan\")\n"
"\t\t\tEND_OBJECT=DataField_7\n"
"\t\t\tOBJECT=DataField_8\n"
"\t\t\t\tDataFieldName=\"Emis_2_err\"\n"
"\t\t\t\tDataType=NC_USHORT\n"
"\t\t\t\tDimList=(\"Along_Track\",\"Along_Scan\")\n"
"\t\t\tEND_OBJECT=DataField_8\n"
"\t\t\tOBJECT=DataField_9\n"
"\t\t\t\tDataFieldName=\"Emis_3_err\"\n"
"\t\t\t\tDataType=NC_UBYTE\n"
"\t\t\t\tDimList=(\"Along_Track\",\"Along_Scan\")\n"
"\t\t\tEND_OBJECT=DataField_9\n"
"\t\t\tOBJECT=DataField_10\n"
"\t\t\t\tDataFieldName=\"Emis_4_err\"\n"
"\t\t\t\tDataType=NC_USHORT\n"
"\t\t\t\tDimList=(\"Along_Track\",\"Along_Scan\")\n"
"\t\t\tEND_OBJECT=DataField_10\n"
"\t\t\tOBJECT=DataField_11\n"
"\t\t\t\tDataFieldName=\"Emis_5_err\"\n"
"\t\t\t\tDataType=NC_USHORT\n"
"\t\t\t\tDimList=(\"Along_Track\",\"Along_Scan\")\n"
"\t\t\tEND_OBJECT=DataField_11\n"
"\t\t\tOBJECT=DataField_12\n"
"\t\t\t\tDataFieldName=\"LST\"\n"
"\t\t\t\tDataType=NC_USHORT\n"
"\t\t\t\tDimList=(\"Along_Track\",\"Along_Scan\")\n"
"\t\t\tEND_OBJECT=DataField_12\n"
"\t\t\tOBJECT=DataField_13\n"
"\t\t\t\tDataFieldName=\"LST_err\"\n"
"\t\t\t\tDataType=NC_UBYTE\n"
"\t\t\t\tDimList=(\"Along_Track\",\"Along_Scan\")\n"
"\t\t\tEND_OBJECT=DataField_13\n"
"\t\t\tOBJECT=DataField_14\n"
"\t\t\t\tDataFieldName=\"PWV\"\n"
"\t\t\t\tDataType=NC_SHORT\n"
"\t\t\t\tDimList=(\"Along_Track\",\"Along_Scan\")\n"
"\t\t\tEND_OBJECT=DataField_14\n"
"\t\t\tOBJECT=DataField_15\n"
"\t\t\t\tDataFieldName=\"QC\"\n"
"\t\t\t\tDataType=NC_USHORT\n"
"\t\t\t\tDimList=(\"Along_Track\",\"Along_Scan\")\n"
"\t\t\tEND_OBJECT=DataField_15\n"
"\t\t\tOBJECT=DataField_16\n"
"\t\t\t\tDataFieldName=\"View_angle\"\n"
"\t\t\t\tDataType=NC_UBYTE\n"
"\t\t\t\tDimList=(\"Along_Track\",\"Along_Scan\")\n"
"\t\t\tEND_OBJECT=DataField_16\n"
"\t\t\tOBJECT=DataField_17\n"
"\t\t\t\tDataFieldName=\"oceanpix\"\n"
"\t\t\t\tDataType=NC_UBYTE\n"
"\t\t\t\tDimList=(\"Along_Track\",\"Along_Scan\")\n"
"\t\t\tEND_OBJECT=DataField_17\n"
"\t\tEND_GROUP=DataField\n"
"\t\tGROUP=SwathAttributes\n"
"\t\tEND_GROUP=SwathAttributes\n"
"\tEND_GROUP=SWATH_1\n"
"END_GROUP=SwathStructure\n"
"GROUP=GridStructure\n"
"END_GROUP=GridStructure\n"
"GROUP=PointStructure\n"
"END_GROUP=PointStructure\n"
"END\n";


char *getLSTEStructMetadata(size_t rows, size_t cols)
{
    int i;
    char buff[32];
    char *findstr = strstr(StructMetadata_LSTE, "####");
    if (findstr)
    {
        sprintf(buff, "%4lu", rows);
        for (i = 0; i < 4; i++)
            *findstr++ = buff[i];
    }
    findstr = strstr(StructMetadata_LSTE, "####");
    if (findstr)
    {
        sprintf(buff, "%4lu", cols);
        for (i = 0; i < 4; i++)
            *findstr++ = buff[i];
    }
    return StructMetadata_LSTE;
}

int writeDatasetMetadataHdf5(hid_t h5fid, const char * dataset, 
        const char *long_name, const char *units,
        int32 range_type, const void *range_min, const void *range_max, 
        int32 fill_type, void *fill_value, 
        double scale, double offset)
{
    hid_t dataset_id = H5Dopen2(h5fid, dataset, H5P_DEFAULT);
    if (dataset_id < 0) return dataset_id;

    // Make a set of data for the dataset.
    MetadataSet dsmeta;
    constructMetadataSet(&dsmeta, dataset_id);

    // Attribute for long_name
    MetadataAttribute *dsatt =
        (MetadataAttribute*)malloc(sizeof(MetadataAttribute));;
    constructMetadataAttribute(dsatt);
    strcpy(dsatt->name, "long_name");
    dsatt->type = DFNT_CHAR8;
    dsatt->value = strdup(long_name);
    dsatt->count = strlen(dsatt->value);
    addAttribute(&dsmeta, dsatt);

    dsatt = (MetadataAttribute*)malloc(sizeof(MetadataAttribute));;
    constructMetadataAttribute(dsatt);
    strcpy(dsatt->name, "units");
    dsatt->type = DFNT_CHAR8;
    dsatt->value = strdup(units);
    dsatt->count = strlen(dsatt->value);
    addAttribute(&dsmeta, dsatt);

    dsatt = (MetadataAttribute*)malloc(sizeof(MetadataAttribute));;
    constructMetadataAttribute(dsatt);
    strcpy(dsatt->name, "format");
    dsatt->type = DFNT_CHAR8;
    bool scaled = false;
    if (scale != 0.0 && scale != 1.0)
    {
        dsatt->value = strdup("scaled");
        scaled = true;
    }
    else
    {
        dsatt->value = strdup("unscaled");
    }
    dsatt->count = strlen(dsatt->value);
    addAttribute(&dsmeta, dsatt);

    dsatt = (MetadataAttribute*)malloc(sizeof(MetadataAttribute));;
    constructMetadataAttribute(dsatt);
    strcpy(dsatt->name, "coordsys");
    dsatt->type = DFNT_CHAR8;
    dsatt->value = strdup("cartesian");
    dsatt->count = strlen(dsatt->value);
    addAttribute(&dsmeta, dsatt);

    dsatt = (MetadataAttribute*)malloc(sizeof(MetadataAttribute));;
    constructMetadataAttribute(dsatt);
    strcpy(dsatt->name, "valid_range");
    dsatt->type = range_type;
    size_t rsize = get_hdf4_size(range_type);
    dsatt->value = (char *)malloc(2 * rsize);
    if (range_type == DFNT_FLOAT32)
    {
       float32 *pf32 = (float32*)dsatt->value;
       pf32[0] = *(float32*)range_min;
       pf32[1] = *(float32*)range_max;
    }
    else if (range_type == DFNT_FLOAT64)
    {
       double *pf64 = (double*)dsatt->value;
       pf64[0] = *(double*)range_min;
       pf64[1] = *(double*)range_max;
    }
    else if (range_type == DFNT_UINT16)
    {
       uint16 *pu16 = (uint16*)dsatt->value;
       pu16[0] = *(uint16*)range_min;
       pu16[1] = *(uint16*)range_max;
    }
    else if (range_type == DFNT_INT16)
    {
       int16 *pi16 = (int16*)dsatt->value;
       pi16[0] = *(int16*)range_min;
       pi16[1] = *(int16*)range_max;
    }
    else if (range_type == DFNT_UINT8)
    {
       uint8 *pu8 = (uint8*)dsatt->value;
       pu8[0] = *(uint8*)range_min;
       pu8[1] = *(uint8*)range_max;
    }
    else
    {
        fprintf(stderr, "Unrecognized type for valid_range attribute of %s\n", dataset);
        return -1;
    }
    dsatt->count = 2;
    addAttribute(&dsmeta, dsatt);

    // If a fill type is provided:
    if (fill_type != DFNT_NONE)
    {
        dsatt = (MetadataAttribute*)malloc(sizeof(MetadataAttribute));;
        constructMetadataAttribute(dsatt);
        strcpy(dsatt->name, "_FillValue");
        dsatt->type = fill_type;
        size_t dsize;
        if (fill_type == DFNT_CHAR8)
        {
            dsatt->count = strlen((const char*)fill_value);
            dsize = dsatt->count + 1;
        }
        else
        {
            dsatt->count = 1;
            dsize = get_hdf4_size(fill_type);
        }
        dsatt->value = (char *)malloc(dsize);
        memcpy(dsatt->value, fill_value, dsize);
        addAttribute(&dsmeta, dsatt);
    }

    if (scaled)
    {
        dsatt = (MetadataAttribute*)malloc(sizeof(MetadataAttribute));;
        constructMetadataAttribute(dsatt);
        strcpy(dsatt->name, "scale_factor");
        dsatt->type = DFNT_FLOAT64;
        dsatt->value = (char *)malloc(sizeof(double));
        double *vscale = (double*)dsatt->value;
        vscale[0] = scale;
        dsatt->count = 1;
        addAttribute(&dsmeta, dsatt);
    }

    if (scaled || offset != 0.0)
    {
        dsatt = (MetadataAttribute*)malloc(sizeof(MetadataAttribute));;
        constructMetadataAttribute(dsatt);
        strcpy(dsatt->name, "add_offset");
        dsatt->type = DFNT_FLOAT64;
        dsatt->value = (char *)malloc(sizeof(double));
        double *voff = (double*)dsatt->value;
        voff[0] = offset;
        dsatt->count = 1;
        addAttribute(&dsmeta, dsatt);
    }

    int stat = writeMetadataHdf5(&dsmeta, dataset_id);
    H5Dclose(dataset_id);
    destructMetadataSet(&dsmeta);
    return stat;
}
