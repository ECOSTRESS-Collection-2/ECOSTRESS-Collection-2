#pragma once
/**
 * @file This header file contains definitions and prototypes for the functions 
 * that handle error, warning, and info messages.
 *
 * @author Robert Freepartner, JaDa Systems/Raytheon/JPL
 *
 * @copyright (c) Copyright 2015, Jet Propulsion Laboratories, Pasadena, CA
 */

/**
 * @brief Log message severity
 */
typedef enum
{
    Info,
    Warning,
    Error
} Severity;

/**
 * @brief sets the log file name
 *
 * @param logname   name of output log file
 */
void set_log_filename(const char *logname);

/**
 * @brief enable echo all log mesages to stdout.
 */
void echo_log_to_stdout();

/**
 * @brief Output a message. Exit if severity is error.
 *
 * This function accepts a variable number of args so that
 * the message string may be a fprintf-like format string 
 * followed by arguments to match the format items. Note that
 * newline characters are not needed.
 *
 * @param   severe  Error, Warning, or Log. Error causes exit.
 * @param   srcfile source file name
 * @param   lineno  line number in source file
 * @param   ernum   error number (used with exit for Error severity.)
 * @param   message fprintf style message format string
 * @param   ...     variable arguments for message format string
 */
void logmsg(Severity severe, const char *srcfile, int linenum, 
            int ernum, const char *message, ...);

// In plain old C, the ... part of the macro cannot be empty.
// NULL can be used as a dumy, e.g., 
//    ERROR("Bad thing happened", NONE);
#define NONE 0

// ERREXIT Macro
#define ERREXIT(ERRNO, message, args...) \
    logmsg(Error, (__FILE__), (__LINE__), (ERRNO), (message), args)

// ERROR Macro
#define ERROR(message, args...) \
    logmsg(Error, (__FILE__), (__LINE__), (__LINE__), (message), args)

// WARNING Macro
#define WARNING(message, args...) \
    logmsg(Warning, (__FILE__), (__LINE__), 0, (message), args)

// INFO Macro
#define INFO(message, args...) \
    logmsg(Info, (__FILE__), (__LINE__), 0, (message), args)
