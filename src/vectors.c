/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "robo.h"

/* Lazy handler for array math. */


double *vector_add(double *A, double *B, int N) {
  double *O = malloc(sizeof(double) * N);
  if (O == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  int i;

  for (i = 0; i < N ; i++) {
    O[i] = A[i] + B[i];
  }

  return(O);
}

double *vector_subtract(double *A, double *B, int N) {
  double *O = malloc(sizeof(double) * N);
  if (O == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  int i;

  for (i = 0; i < N ; i++) {
    O[i] = A[i] - B[i];
  }

  return(O);
}

double *vector_subtract_constant(double *A, double B, int N) {
  double *O = malloc(sizeof(double) * N);
  if (O == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  int i;

  for (i = 0; i < N ; i++) {
    O[i] = A[i] - B;
  }

  return(O);
}

double *vector_multiply(double *A, double *B, int N) {
  double *O = malloc(sizeof(double) * N);
  if (O == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory for O\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  int i;

  for (i = 0; i < N ; i++) {
    O[i] = A[i] * B[i];
  }

  return(O);
}

double *vector_divide(double *A, double *B, int N) {
  double *O = malloc(sizeof(double) * N);
  if (O == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  int i;

  for (i = 0; i < N ; i++) {
    O[i] = A[i] / B[i];
  }

  return(O);
}

double *vector_constant(double V, int N) {
  double *O = malloc(sizeof(double) * N);
  if (O == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  int i;
  for (i = 0; i < N ; i++) {
    O[i] = V;
  }

  return(O);
}

double *vector_copy(double *A, int N) {
  double *O = malloc(sizeof(double) * N);
  if (O == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  int i;
  for (i = 0; i < N; i++) {
    O[i] = A[i];
  }

  return(O);
}
