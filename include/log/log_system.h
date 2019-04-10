/*
 * ErrorHandling.h
 *
 *  Created on: Apr 29, 2013
 *      Author: ivan
 */

#ifndef ERRORHANDLING_H_
#define ERRORHANDLING_H_

#include <stdarg.h>
#include <string>
#include <fstream>
#include <time.h>
#include <stdarg.h>

#define LOG_NEWLINE         LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, LogSystem::FormatString("\n"), std::string("[]"));
#define LOG_NOHEADER(...)         LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, LogSystem::FormatString(__VA_ARGS__), std::string("[]"));
#define LOG_ALL(...)        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, LogSystem::FormatString(__VA_ARGS__), std::string(__FUNCTION__)); // , std::string(__FILE__), __LINE__);
#define LOG_MEDHIGH(...)        LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED | VERBOSE_LEVEL_HIGH, true, LogSystem::FormatString(__VA_ARGS__), std::string(__FUNCTION__)); // , std::string(__FILE__), __LINE__);
#define LOG_MEDHIGH_NOHEADER(...)         LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED | VERBOSE_LEVEL_HIGH, true, LogSystem::FormatString(__VA_ARGS__), std::string("[]"));
#define LOG_HIGH(...)        LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH, true, LogSystem::FormatString(__VA_ARGS__), std::string(__FUNCTION__)); // , std::string(__FILE__), __LINE__);
#define LOG_HIGH_NOHEADER(...)         LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH, true, LogSystem::FormatString(__VA_ARGS__), std::string("[]"));
#define LOG_DEBUG(...)      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, true, LogSystem::FormatString(__VA_ARGS__), std::string(__FUNCTION__), std::string(__FILE__), __LINE__);
#define LOG_DEBUG_MEDHIGH(...)      LogSystem::GetInstance().Log(VERBOSE_LEVEL_MEDHIGH_DEBUG, true, LogSystem::FormatString(__VA_ARGS__), std::string(__FUNCTION__), std::string(__FILE__), __LINE__);
#define LOG_DEBUG_HIGH(...)      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, true, LogSystem::FormatString(__VA_ARGS__), std::string(__FUNCTION__), std::string(__FILE__), __LINE__);
#define LOG_DEBUG_NOHEADER(...)      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, true, LogSystem::FormatString(__VA_ARGS__), std::string("[]"), std::string(__FILE__), __LINE__);
#define LOG_DEBUG_SPEC(...) LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, LogSystem::FormatString(__VA_ARGS__), std::string(__FUNCTION__), std::string(__FILE__), __LINE__);
#define LOG_DEBUG_SPEC_NO_HEADER(...) LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, LogSystem::FormatString(__VA_ARGS__), std::string("[]"), std::string(__FILE__), __LINE__);
#define LOG_DEBUG_SPEC_NEWLINE LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, "\n", std::string("[]"));

#define WARNING_REPORT(error_code, ...)   LogSystem::GetInstance().Error(SEVERITY_INT_WARNING, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(error_code, __VA_ARGS__))
#define ERROR_REPORT(error_code, ...)   LogSystem::GetInstance().Error(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(error_code, __VA_ARGS__))
#define FATAL_REPORT(error_code, ...)   LogSystem::GetInstance().Error(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(error_code, __VA_ARGS__))

// LogSystem class is used to standardize the error reporting process.
// Sample usage:
//    LOG(FATAL) << LogSystem::GenerateErrorMessage(ERR_MEMORY, "Parameter that caused the crash: fp = %ld", (int64_t) fp);
//
class LogSystem {
 public:
  LogSystem();
  ~LogSystem();

  static LogSystem& GetInstance();

  // GenerateErrorMessage takes the error number and an additional message
  // formatted as C-style string (same as printf's formatted string). It then
  // compiles an error message, converted to a C++ string, that can be given to
  // a logging unit for output.
  // Sample usage:
  //  std::string error_message = ErrorHandling::GenerateErrorMessage (ERR_OPENING_FILE, "Could not open file: '%s'", file_path.c_str());
  std::string GenerateErrorMessage(uint32_t error_type,
                                          const char *additional_message, ...);
  int WriteLog(std::string log_entry, bool always_output_to_std=false);
  int Error(int severity, std::string function, std::string message);
  int Log(uint32_t verbose_level, bool trigger_condition, std::string message, std::string calling_func="", std::string calling_file="", int32_t calling_line=-1);
  void SetProgramVerboseLevelFromInt(int64_t verbose_level);

  std::string GetLocalTime();
  std::string GetUTCTime(std::string fmt="%a, %d %b %y %T %z");
  // The function uses vasprintf for formatting the function arguments.
  // This option was chosen because it preallocates the memory needed to store
  // the output formatted string, unlike most other *printf alternatives
  static std::string FormatString(const char* additional_message, ...);

  std::string LOG_FILE;
  uint32_t LOG_VERBOSE_TYPE;
  uint32_t PROGRAM_VERBOSE_LEVEL;
};

#define LOG_VERBOSE_STD     1 << 0
#define LOG_VERBOSE_FILE    1 << 1
#define LOG_VERBOSE_FULL    LOG_VERBOSE_STD | LOG_VERBOSE_FILE

#define VERBOSE_LEVEL_LOW       1 << 0
#define VERBOSE_LEVEL_MED       1 << 1
#define VERBOSE_LEVEL_HIGH      1 << 2
#define VERBOSE_LEVEL_DEBUG     1 << 3
#define VERBOSE_LEVEL_ALL       (VERBOSE_LEVEL_LOW | VERBOSE_LEVEL_MED | VERBOSE_LEVEL_HIGH)
#define VERBOSE_LEVEL_ALL_DEBUG (VERBOSE_LEVEL_HIGH | VERBOSE_LEVEL_MED | VERBOSE_LEVEL_LOW | VERBOSE_LEVEL_DEBUG)
#define VERBOSE_LEVEL_LOW_DEBUG (VERBOSE_LEVEL_LOW | VERBOSE_LEVEL_DEBUG)
#define VERBOSE_LEVEL_MED_DEBUG (VERBOSE_LEVEL_MED | VERBOSE_LEVEL_DEBUG)
#define VERBOSE_LEVEL_HIGH_DEBUG (VERBOSE_LEVEL_HIGH | VERBOSE_LEVEL_DEBUG)
#define VERBOSE_LEVEL_MEDHIGH_DEBUG (VERBOSE_LEVEL_MED | VERBOSE_LEVEL_HIGH | VERBOSE_LEVEL_DEBUG)
#define VERBOSE_LEVEL_FORCE     1 << 4
#define VERBOSE_FREQ_LOW        1 << 5
#define VERBOSE_FREQ_MED        1 << 6
#define VERBOSE_FREQ_HIGH       1 << 7
#define VERBOSE_FREQ_ALL        (VERBOSE_FREQ_LOW | VERBOSE_FREQ_MED | VERBOSE_FREQ_HIGH)

#define SEVERITY_INT_INFO     1 << 10
#define SEVERITY_INT_WARNING  1 << 11
#define SEVERITY_INT_ERROR    1 << 12
#define SEVERITY_INT_FATAL    1 << 13



#define SEVERITY_INFO     ((std::string) "INFO")
#define SEVERITY_WARNING  ((std::string) "WARNING")
#define SEVERITY_ERROR    ((std::string) "ERROR")
#define SEVERITY_FATAL    ((std::string) "FATAL")

#define ERR_MEMORY            1
#define ERR_OPENING_FILE        2
#define ERR_CLOSING_FILE        3
#define ERR_OTHER               4
#define ERR_UNEXPECTED_VALUE    5

#define ERR_FILE_WRITE_DATA       6
#define ERR_FILE_READ_DATA        7
#define ERR_WRONG_FILE_TYPE       8
#define ERR_NOT_IMPLEMENTED       9
#define ERR_SEQUENCE_MISMATCH     10
#define ERR_WRONG_PARAMS        11
#define ERR_FILE_DEFORMED_FORMAT 12
#define ERR_FILE_NOT_FOUND      13
#define ERR_FOLDER_NOT_FOUND    14
#define ERR_GENERIC    15

#endif /* ERRORHANDLING_H_ */
