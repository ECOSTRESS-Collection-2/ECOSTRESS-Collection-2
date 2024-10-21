// config.h
#pragma once
/** 
 * @file This is the header for the config package. This package uses
 * expat to parse the XML RunConfig file and provides an API to search 
 * for values by group or overall.
 *
 * @author Robert Freepartner, JPL/Raytheon/JaDa Systems, April 2016
 *
 * @copyright (c) 2016 JPL, All rights reserved
 *
 * Example: Get the L1B_RAD input file name from group InputFileGroup.
 *
 * ConfigSet *RunConfig = parseConfig("RunConfig.xml");
 * assert(RunConfig != NULL);
 * assert (RunConfig->status >= 0);
 * ConfigElement *l1b_rad = 
 *     findGroupElement(RunConfig, "InputFileGroup", "L1B_RAD");
 * assert(l1b_rad != NULL);
 * const char *input_filename = getValue(l1b_rad, 0);
 * // Get other values.
 * destructConfigSet(RunConfig);
 */
#include "bool.h"
#define CONFIG_MAX_VALUES 32
#define CONFIG_MAX_NAMELEN 64
#define CONFIG_VALUELEN 256
#ifndef PATH_MAX
#define PATH_MAX 1024
#endif

/**
 * @brief Container for a single config element
 */
typedef struct ConfigElement
{
    char    *name;  ///< name attribute from scalar or vector tag.
    int     count;  ///< Number of values. For scalar this is always 1.
    char    *value[CONFIG_MAX_VALUES];  ///< Array of pointers to input values.
    struct ConfigElement *next;         ///< Link pointer.
} ConfigElement;

void constructConfigElement(ConfigElement *ce, const char *_name);
void destructConfigElement(ConfigElement *ce);

/**
 * @brief Container for a config group
 */
typedef struct ConfigGroup
{
    char    *name;  ///< the name attribute from the group tag
    int     count;  ///< number of elements in the group
    ConfigElement   *head;  ///< first element in the group
    ConfigElement   *tail;  ///< last element in the group
    struct ConfigGroup *next;   ///< link pointer
} ConfigGroup;

void constructConfigGroup(ConfigGroup *cg, const char *_name);
void destructConfigGroup(ConfigGroup *cg);

/**
 * @brief Head structure for a configuration set
 */
typedef struct
{
    int             count;   ///< number of groups in config
    int             verbose; ///< 0=no output, 3=max
    int             status;  ///< 0=success -1=I/O error -2=parse error
    ConfigGroup     *head_group; ///< first element in the group
    ConfigGroup     *tail_group; ///< last element in the group
} ConfigSet;

void constructConfigSet(ConfigSet *cs, int in_verbose);
void destructConfigSet(ConfigSet *cs);

/**
 * @brief Add a new group to a set.
 *
 * @param   cs          ConfigSet pointer
 * @param   group_name  Name of the group
 * @return  pointer to new ConfigGroup
 */
ConfigGroup *addGroup(ConfigSet *cs, const char *group_name);

/**
 * @brief Add a new element to a group.
 *
 * @param   cg              ConfigGroup pointer
 * @param   element_name    Name of the element
 * @return  pointer to new ConfigElement
 */
ConfigElement *addElement(ConfigGroup *cg, const char *element_name);

/**
 * @brief Get a value for an element in a group.
 *
 * @param   cg              Pointer to ConfigGroup
 * @param   element_name    Name of the element.
 * @return  NULL=not found, else, pointer to ConfigElement
 */
ConfigElement *findElementInGroup(ConfigGroup *cg, const char *element_name);

/**
 * @brief find an element in a specific group.
 *
 * @param   ce          ///< ConfigElement pointer
 * @param   n           ///< Index to element, 0 for scalar
 * @return  NULL=not found, else pointer to the element
 */
ConfigElement *findGroupElement(ConfigSet *cs, const char *group_name,
    const char *element_name);

/**
 * @brief Search groups to get a value for an element in the set.
 *
 * @param   cs              ConfigSet pointer
 * @param   group_name      Name of group. NULL for default.
 * @param   element_name    Name of the element.
 * @return NULL=not found, else, pointer to ConfigElement
 */
ConfigElement *findElement(ConfigSet *cs, const char *element_name);

/**
 * @brief Get a named group
 *
 * @param   cs          ///< ConfigSet pointer
 * @param   group_name  ///< Name of group
 * @return  NULL=invalid index, else pointer to ConfigGroup
 */
ConfigGroup *getGroup(ConfigSet *cs, const char *group_name);

/**
 * @brief Get value number n from an element, where n < count
 *
 * @param   ce          ///< ConfigElement pointer
 * @param   n           ///< Index to element, 0 for scalar
 * @return  NULL=invalid index, else pointer to value string
 */
char *getValue(ConfigElement *ce, int n);

 /**
  * @brief Parses the RunConfig file and stores the config elements 
  * by group in the set. Any elements not in a group will be stored
  * in a default group.
  *
  * @param  filename    RunConfig filename
  * @return pointer to new set. Call deconstructConfigSet before 
  * freeing the pointer when done using it.
  */
ConfigSet *parseConfig(const char *filename);

/**
 * @brief Set the verbose level for parsing
 *
 * @param   _verbose    0=no messages, 1=errors only, 2=parser progress, 3=all
 */
void setConfigVerbose(int _verbose);
