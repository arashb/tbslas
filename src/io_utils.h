#include <iostream>
#include <string>
#include <stdio.h>
#include <time.h>

namespace tbslas {

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string get_current_datetime() {
  time_t     now = time(0);
  struct tm  tstruct;
  char       buf[80];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
  return buf;
}

const std::string get_result_dir() {
  char res_dir_buffer[200];
  snprintf(res_dir_buffer, sizeof(res_dir_buffer),
           "%s/%s/", getenv("WORK"), getenv("TBSLAS_RES_DIR_NAME"));

  return res_dir_buffer;
}

}  // namespace tbslas
