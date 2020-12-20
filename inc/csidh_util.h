#ifndef CSIDH_UTIL_H
#define CSIDH_UTIL_H

#include "fp.h"
#include "edwards_curve.h"

extern const fp R_mod_p;

#define VERSION 0.1

void pprint_pk(void *x);

void pprint_sk(void *y);

void pprint_ss(uint64_t *x);

void save_file(char *file, void *buf, size_t len);

void error_exit(char *str);

int read_file(const char *file, uint8_t *buf, size_t len);

int read_stdin(uint8_t *buf, int len);
#endif 
