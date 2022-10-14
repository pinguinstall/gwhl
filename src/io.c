#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

/**
* @brief writes like fprintf but uses two separate outputs
* 
* @param stream1 p_stream1: first stream to write to
* @param stream2 p_stream2: second stream to write to
* @param format p_format: format string
* @return int - 1 or 2 if failing to stream 1 or 2 fails, 0 on success
*/
int GWHL_fprintf2(FILE* stream1, FILE* stream2, const char *format, ...){
  va_list args1;
  va_list args2;
  // we need two arglists here, since vfprintf clobbers that list
  va_start(args1, format);
  va_start(args2, format);
  
  if(vfprintf(stream1, format, args1) < 0 ){
      perror("could not write to file");
      return 1;
  }
    
  if(vfprintf(stream2, format, args2) < 0 ){
      perror("could not write to file");
      return 2;
  }
  
  va_end(args1);
  va_end(args2);
  
  return 0;
}
