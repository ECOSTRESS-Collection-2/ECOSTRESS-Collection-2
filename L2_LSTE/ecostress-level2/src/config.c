/**
 * See config.h for documentation.
 * @author Robert Freepartner, JPL/Raytheon/JaDa Systems, April 2016
 * @copyright (c) 2016 JPL, All rights reserved
 *
 */
 
#include "config.h"
#include <expat.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Static data for the parsing procedure.
static ConfigSet *config_set = NULL;
static ConfigGroup *current_group = NULL;
static ConfigElement *current_element = NULL;
static bool snag_value = false;
static int verbose_level = 1; // default = only error messages

bool is_white(char ch)
{
    if (ch > ' ' && ch <= '~') return false;
    return true;
}

const char *strtrim(const char *str, int *len)
{
    // trim front
    while (*len > 0 && is_white(*str))
    {
        str++;
        (*len)--;
    }
    // trim back
    while (*len > 0 && is_white(*(str + *len - 1)))
    {
        (*len)--;
    }
    return str;
}
void constructConfigElement(ConfigElement *ce, const char *_name)
{
   assert(ce != NULL);
   ce->name = strdup(_name);;
   ce->count = 0;
   ce->next = NULL;    
}

void destructConfigElement(ConfigElement *ce)
{
    if (ce == NULL)
        return;
    int i;
    for (i = 0; i < ce->count; i++)
    {
        free(ce->value[i]);
        ce->value[i] = NULL;
    }
    if (ce->name != NULL) free(ce->name);
    ce->name = NULL;
    if (ce->next != NULL)
        destructConfigElement(ce->next);
    ce->next = NULL;
    ce->count = 0;
}

void constructConfigGroup(ConfigGroup *cg, const char *_name)
{
    assert(cg != NULL);
    cg->name = strdup(_name);
    cg->count = 0;
    cg->head = NULL;
    cg->tail = NULL;
    cg->next =  NULL;
}

void destructConfigGroup(ConfigGroup *cg)
{
    if (cg == NULL);
        return;
    if (cg->name != NULL) free(cg->name);
    cg->name = NULL;
    if (cg->count > 0)
        destructConfigElement(cg->head);
    cg->head = NULL;
    cg->tail = NULL;
    cg->count = 0;
    if (cg->next != NULL)
        destructConfigGroup(cg->next);
    cg->next = NULL;
}

void constructConfigSet(ConfigSet *cs, int in_verbose)
{
    assert(cs != NULL);
    cs->count = 0;
    cs->verbose = in_verbose;
    cs->status = 0;
    cs->head_group = NULL;
    cs->tail_group = NULL;
}

void destructConfigSet(ConfigSet *cs)
{
    if (cs == NULL)
        return;
    if (cs->head_group != NULL)
        destructConfigGroup(cs->head_group);
    constructConfigSet(cs, 0); // This NULLs out the pointers for safety.
}

ConfigGroup *addGroup(ConfigSet *cs, const char *group_name)
{
    if (cs == NULL)
        return NULL;
    ConfigGroup *cg = (ConfigGroup*)malloc(sizeof(ConfigGroup));
    assert(cg != NULL);
    constructConfigGroup(cg, group_name);
    if (cs->head_group == NULL)
    {
        cs->head_group = cg;
    }
    else
    {
        cs->tail_group->next = cg;
    }
    cs->tail_group = cg;
    cs->count++;
    
    if (cs->verbose > 1)
    {
        printf("Added <group> %s\n", group_name);
    }
    return cg;
}

ConfigElement *addElement(ConfigGroup *cg, const char *element_name)
{
    if (cg == NULL)
        return NULL;
    ConfigElement *ce = (ConfigElement*)malloc(sizeof(ConfigElement));
    assert(ce != NULL);
    constructConfigElement(ce, element_name);
    if (cg->head == NULL)
    {
        cg->head = ce;
    }
    else
    {
        cg->tail->next = ce;
    }
    cg->tail = ce;
    cg->count++;

    return ce;
}

ConfigElement *findElementInGroup(ConfigGroup *cg, const char *element_name)
{
    if (cg == NULL)
        return NULL;
    ConfigElement *ce = cg->head;
    while (ce != NULL)
    {
        if (strcmp(ce->name, element_name) == 0)
            return ce;
        ce = ce->next;
    }
    return NULL;
}

ConfigElement *findGroupElement(ConfigSet *cs, const char *group_name,
    const char *element_name)
{
    if (cs == NULL)
        return NULL;
    ConfigGroup *cg = cs->head_group;
    while (cg != NULL)
    {
        if (strcmp(cg->name, group_name) == 0)
        {
            return findElementInGroup(cg, element_name);
        }
        cg = cg->next;
    }
    return NULL;
}

ConfigElement *findElement(ConfigSet *cs, const char *element_name)
{
    if (cs == NULL)
        return NULL;
    ConfigGroup *cg = cs->head_group;
    ConfigElement *ce;
    while (cg != NULL)
    {
        ce = findElementInGroup(cg, element_name);
        if (ce != NULL) return ce;
        cg = cg->next;
    }
    return NULL;
}

ConfigGroup *getGroup(ConfigSet *cs, const char *group_name)
{
    if (cs == NULL)
        return NULL;
    ConfigGroup *cg = cs->head_group;
    while (cg)
    {
        if (strcmp(cg->name, group_name) == 0)
        {
            return cg;
        }
        cg = cg->next;
    }
    return NULL;
}

char *getValue(ConfigElement *ce, int n)
{
    if (ce == NULL)
        return NULL;
    if (n >= ce->count)
        return NULL;
    return ce->value[n];
}


// EXPAT XML Parser handlers:


// This handler is called when a start tag is found, e.g., <group>
void start_handler(void *data, const char *el, const char **attr)
{
    snag_value = false;
    if (strcmp(el, "group") == 0)
    {
        // <group> -- add a new group
        current_group = addGroup(config_set, attr[1]);
    }
    else if (strcmp(el, "scalar") == 0)
    {
        // Single value element
        snag_value = true;
        if (current_group == NULL)
        {
            // Found an element before a group, so create a default group.
            current_group = getGroup(config_set, "default");
            if (current_group == NULL)
                current_group = addGroup(config_set, "default");
        }
        current_element = addElement(current_group, attr[1]);
        if (config_set->verbose > 1)
        {
            printf("Added <scalar> %s to group %s\n",
                    attr[1], current_group->name);
        }
    }
    else if (strcmp(el, "vector") == 0)
    {
        // Multiple value element
        if (current_group == NULL)
        {
            // Found an element outside of a group, so create/use a default group.
            current_group = getGroup(config_set, "default");
            if (current_group == NULL)
                current_group = addGroup(config_set, "default");
        }
        current_element = addElement(current_group, attr[1]);
        if (config_set->verbose > 1)
        {
            printf("Added <vector> %s to group %s\n",
                    attr[1], current_group->name);
        }
    }
    else if (strcmp(el, "element") == 0)
    {
        // One of the values of the multi-value element
        if (current_element == NULL)
        {
            config_set->status = -2;
            if (config_set->verbose > 0)
                fprintf(stderr, "XML Parser error: found tag <element> not "
                    "preceded by a <vector> tag.\n");
        }
        else
            snag_value = 1;
    }
    else if (strcmp(el, "input") == 0)
    {
        // Ignore this tag
    }
    else
    {
        if (config_set->verbose > 2)
        {
            printf("Unrecognized tag in input file: <%s>\n", el);
        }
    }
}

// This handler is called when the value following a start tag is found.
void value_handler(void *user_data, const char *val, int len)
{
    if (!snag_value) return;
    const char *trimval = strtrim(val, &len);
    if (len < 1) return; // Skip blank strings
    snag_value = 0;
    assert(current_element != NULL); // This can "never" happen.
    int i = current_element->count;
    current_element->count++;
    if (current_element->count > CONFIG_MAX_VALUES)
    {
        config_set->status = -2;
        if (config_set->verbose > 0)
            fprintf(stderr, 
                    "XML Parser error: too many values for %s in group %s\n",
                    current_element->name, current_group->name);
        current_element->count--;
        return;
    }
    current_element->value[i] = (char*)malloc(len+1);
    assert(current_element->value[i] != NULL);
    memcpy(current_element->value[i], trimval, len);
    current_element->value[i][len] = '\0';
    if (config_set->verbose > 2)
    {
        printf("  Value[%d] set to %s\n",
                current_element->count, current_element->value[i]);
    }
}

// This handler is called for an end tag, e.g., </group>
void end_handler(void *data, const char *el)
{
    if (strcmp(el, "scalar") == 0 ||
        strcmp(el, "vector") == 0)
    {
        // End of current_element
        current_element = NULL;
    }
    if (strcmp(el, "group") == 0)
    {
        // End of current_group
        current_group = NULL;
    }
}

ConfigSet *parseConfig(const char *filename)
{
    // This is the big one. It reads in the RunConfig file and parses
    // the contents using expat.
    int stat;

    // Initialize the ConfigSet
    ConfigSet *cs = (ConfigSet*)malloc(sizeof(ConfigSet));
    assert(cs != NULL);
    config_set = cs; // Set the static pointer for the parser handlers.
    constructConfigSet(config_set, verbose_level);
    current_group = NULL;
    current_element = NULL;

    // Initialize the parser.
    XML_Parser ParseRunConfig = XML_ParserCreate(NULL);
    XML_SetElementHandler(ParseRunConfig, start_handler, end_handler);
    XML_SetCharacterDataHandler(ParseRunConfig, value_handler);
    
    // Open the input RunConfig file.
    FILE *pf = fopen(filename, "r");
    if (pf == NULL)
    {
        if (cs->verbose > 0)
            perror(filename);
        cs->status = -1;
        config_set = NULL;
        return cs;
    }

    // Find out the size of the file.
    stat = fseek(pf, 0L, SEEK_END);
    assert(stat >= 0);
    int fsize = ftell(pf);
    stat = fseek(pf, 0L, SEEK_SET);
    assert(stat >= 0);

    // Allocate a buffer to hold the file's contents.    
    char *buffer = (char*)malloc(fsize);
    assert(buffer != NULL);

    // Read the entire file into the buffer.
    size_t nbytes = fread(buffer, fsize, 1, pf);
    if (nbytes < 1)
    {
        if (cs->verbose > 0)
            perror(filename);
        cs->status = -1;
        config_set = NULL;
        return cs;
    }
    fclose(pf);

    // Parse the XML data in the buffer
    stat = XML_Parse(ParseRunConfig, buffer, fsize, 1);
    if (stat == 0)
    {
        if (cs->verbose != 0)
        {
            fprintf(stderr, "XML Parsing error: input file contains bad XML data.\n");
        }
        cs->status = -2;
    }
    XML_ParserFree(ParseRunConfig);
    free(buffer);
    config_set = NULL; // This parser could be called again for a different file.
    return cs;
}


void setConfigVerbose(int _verbose)
{
    verbose_level = _verbose;
}
