/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "robo.h"

char *generate_filename(opts *options, char *rule) {
  static char root[1024];
  char order[16];
  if ((options->max_order != 0)&&(strcasecmp(rule,"LOG") != 0)) {
    snprintf(order,sizeof(order),"_order%03d",options->order);
  }
  else {
    order[0] = 0;
  }

  if ((options->iteration < options->max_iterations - 1)&&(strcasecmp(rule,"LOG") != 0)) {
    if (options->path_base) {
      snprintf(root,sizeof(root),"%s%s_iter%ld",options->path_base,order,options->iteration + 1);
    }
    else {
      snprintf(root,sizeof(root),"%s%s_iter%ld",basename(options->infilename),order,options->iteration + 1);
    }
  }
  else {
    if (options->path_base) {
      snprintf(root,sizeof(root),"%s%s",options->path_base,order);
    }
    else {
      snprintf(root,sizeof(root),"%s%s",basename(options->infilename),order);
    }
  }

  if (strcasecmp(rule,"SPECTRA") == 0) {
    return(strncat(root,".robospect",32));
  }
  else if (strcasecmp(rule,"LINES") == 0) {
    return(strncat(root,".robolines",32));
  }
  else if (strcasecmp(rule,"PLOT") == 0) {
    return(strncat(root,".robo.ps",32));
  }
  else if (strcasecmp(rule,"LOG") == 0) {
    return(strncat(root,".robo.log",32));
  }
  else {
    return(NULL);
  }
}    


       
