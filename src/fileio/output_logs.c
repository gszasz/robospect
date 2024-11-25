/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"
#include <stdarg.h>

#define LOG_COMMENT_LENGTH 256

void log_comment(opts *options, int level, char *format, ...) {
  va_list arg;
  char *comment;
  int status = 0;

  if (options->verbose & level) {
    comment = malloc(sizeof(char) * LOG_COMMENT_LENGTH);
    if (comment == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    va_start(arg,format);
    status = vsnprintf(comment,LOG_COMMENT_LENGTH,format,arg);
    if (status > LOG_COMMENT_LENGTH) {
      fprintf(stderr,"log had leftovers: %d\n",status - LOG_COMMENT_LENGTH);
    }
    
    if (options->log) {
      fprintf(options->log,"%s\n",comment);
      fflush(options->log);
    }

    if (options->verbose & ROBO_VERBOSE_SCREEN) {
      fprintf(stderr,"%s\n",comment);
    }
    va_end(arg);
    free(comment);
  }
}
  
