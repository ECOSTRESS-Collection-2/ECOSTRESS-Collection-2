#pragma once
/**
 * @file Contains definitions and functions used in order to process the 
 * metadata for NPP/VIIRS and other processing systems.
 *
 * @author Robert Freepartner, JPL/Raytheon/JaDa Systems, January 2016
 * 10-JUN-2019 R. Freepartner Provide API to suppress _FillValue field.
 *
 * @copyright (c) 2019 JPL, All rights reserved
 */

#include <hdf.h>
#include <hdf5.h>
#include <mfhdf.h>
#define METADATA_NAME_SIZE FIELDNAMELENMAX
 
 /**
  * @brief Contains information and value for one attribute.
  */
typedef struct MetadataAttribute
{
    int     index;  ///< attribute index
    char    name[METADATA_NAME_SIZE]; ///< Attribute name
    int32   type;   ///< HDF4 data type
    int32   count;  ///< Array size (strlen for char*)
    char   *value;  ///< Buffer for value, technically used as VOIDP
    struct MetadataAttribute *next; ///< Link for list.
} MetadataAttribute;

void constructMetadataAttribute(MetadataAttribute *att);
void destructMetadataAttribute(MetadataAttribute *att);
void setValue(MetadataAttribute *att, int in_type, int in_count, 
                char *in_value);

/**
 * @brief Contains a set of attributes for an object;
 */
typedef struct
{
    int32   objid;      ///< HDF object id
    int     natts;      ///< Number of attributes
    MetadataAttribute *head; ///< List of attributes
    MetadataAttribute *tail; ///< Last attribute
} MetadataSet;

void constructMetadataSet(MetadataSet *set, int32 in_objid);
void destructMetadataSet(MetadataSet *set);

/**
 * @brief Allocates memory and copies a uint8 value to it.
 *
 * @param   ivalue  Value to duplicate in new memory
 * @return  pointer to duplicated int32 value cast as (char*)
 */
char *dupuint8(uint8 ivalue);

/**
 * @brief Allocates memory and copies an int32 value to it.
 *
 * @param   ivalue  Value to duplicate in new memory
 * @return  pointer to duplicated int32 value cast as (char*)
 */
char *dupint32(int32 ivalue);

/**
 * @brief Allocates memory and copies a float32 value to it.
 *
 * @param   dvalue  Value to duplicate in new memory
 * @return  pointer to duplicated float64 value cast as (char*)
 */
char *dupfloat32(float32 dvalue);

/**
 * @brief Allocates memory and copies a float64 value to it.
 *
 * @param   dvalue  Value to duplicate in new memory
 * @return  pointer to duplicated float64 value cast as (char*)
 */
char *dupfloat64(double dvalue);

/**
 * @brief Loads the set of metadata for an object from hdf4.
 *
 * @param   obj_id      Object id in file for metadata 
 * @param[in,out] set   Object to contain set of metadata
 * @return  0=success  -1=error
 */
int loadMetadataHdf4(MetadataSet *set, int32 obj_id);

/**
 * @brief Loads the set of metadata for an object from an
 * hdf5 or netcdf file.
 *
 * @param[in,out] set   Object to contain set of metadata
 * @param   obj_id      Object id in file for metadata 
 * @return  0=success  -1=error
 */
int loadMetadataHdf5(MetadataSet *set, hid_t obj_id);

/**
 * @brief load scalar datasets from an HDF5 group as metadata.
 *
 * @param[in,out] set   Object to contain set of metadata
 * @param   hfileid     File ID for hdf5 file.
 * @param   groupname   Name of group containing metadata
 * @return  0=success  -1=error
 */
int loadMetadataGroup(MetadataSet *set, hid_t hfileid, const char *groupname);

 /**
  * @brief Finds an attribute by name.
  *
  * @param  set         Metadata set
  * @param  attrname    Name to find
  * @return NULL=not found  else, MetadataAttribute address.
  */
MetadataAttribute *findAttribute(MetadataSet *set, const char *attrname);

/**
 * @brief Adds an attribute to a metadata set.  Updates att->index.
 *
 * @param   set     Metadata set to get a new attribute
 * @param   att     Metadata attribute to be added
 */
void addAttribute(MetadataSet *set, MetadataAttribute *att);

/**
 * @brief Marks an attribute as deleted by setting index to -1.
 *
 * @param   att     Metadata attribute to be added
 */
void deleteAttribute(MetadataAttribute *att);

/**
 * @brief Find an attribute in a set. If found, replace the value.
 * Else, add a new attribute.
 *
 * @param   set         Set containing the attribute to be used
 * @param   attr_name   The name of the attribute to look for
 * @param   in_type     Type of the new data
 * @param   in_count    Size of the new data
 * @param   in_value    Buffer containing the new value
 * @return -1=error 0=replaced 1=added
 */
int updateValue(MetadataSet *set, const char *attr_name, int in_type, 
                int in_count, char *in_value);

/**
 * @brief Writes metadata to an object.
 * 
 * @param   obj_id      HDF5 object handle
 * @param   set         MetadataSet to write
 * @return 0=success -1=error
 */
int writeMetadataHdf5(MetadataSet *set, hid_t loc_id);

/**
 * @brief writes meatadata as scalar datasets in an HDF5 group.
 *
 * The group must exist, so create it before calling this function.
 *
 * @param   set         Object to contain set of metadata
 * @param   hfileid     File ID for hdf5 file.
 * @param   groupname   Name of group to contain metadata
 * @return  0=success  -1=error
 */
int writeMetadataGroup(MetadataSet *set, hid_t hfileid, const char *groupname);

/**
 * @return pointer to the StructMetadata.0 string for CLOUD data.
 */
char *getCLOUDStructMetadata();

/**
 * @return pointer to the StructMetadata.0 string for LSTE data.
 */
char *getLSTEStructMetadata(size_t rows, size_t cols);

/**
 * @brief Writes a standard set of attributes to a dataset.
 *
 * @param   h5fid       HDF5 file ID of file open for output.
 * @param   dataset     Fully qualified name of the dataset
 * @param   long_name   Value of long_name attribute
 * @param   units       Name for units, or "n/a"
 * @param   range_type  HDF4 type for range min/max
 * @param   range_min   Minimum valid value after scale/offset
 * @param   range_max   Maximum valid value
 * @param   fill_type   HDF4 type for fill_value, DFNT_NONE = no _FillValue
 * @param   fill_value  Pointer to fill value of type fill_type.
 * @param   scale       Amount to multiply the stored value
 * @param   offset      Amount to add to the stored value after scale.
 * @return -1=error, else success
 */
int writeDatasetMetadataHdf5(hid_t h5fid, const char *dataset, 
        const char *long_name, const char *units, 
        int32 range_type, const void *range_min, const void *range_max, 
        int32 fill_type, void *fill_value, double scale, double offset);
