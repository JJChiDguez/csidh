#include <getopt.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>

#include "csidh_util.h"

// slightly modified csidh from main/csidh.c
static uint8_t csidh(proj out, const uint8_t sk[], const proj in)
{
  if (!validate(in)) {
    return 0;
  };

  action_evaluation(out, sk, in);
  return 1;
};

// slightly modified fp_print from inc/fp.h
void pprint_ss(uint64_t *x)
{
    int NUM = NUMBER_OF_WORDS;
    int TYPE = 0;
    int i;
    for(i=NUM-1; i > -1; i--){
        printf("%.16" PRIX64 "", x[i]);
    }
    printf("\n");
}

void pprint_pk(void *x) {
  for (size_t i = 0; i < sizeof(proj); ++i) {
    printf("%02hhx", ((uint8_t *)x)[i]);
  }
  printf("\n");
};

void pprint_sk(void *x) {
  for (size_t i = 0; i < sizeof(uint8_t[N]); ++i) {
    printf("%02hhx", ((uint8_t *)x)[i]);
  }
  printf("\n");
};

void save_file(char *file, void *buf, size_t len) {
  FILE *fhandle;
  if (file != NULL) {
    fhandle = fopen(file, "w");
    if (fhandle == NULL) {
      fprintf(stderr, "Unable to open %s", file);
      exit(1);
    }
    for (size_t i = 0; i < len; ++i)
      fprintf(fhandle, "%02hhx", ((uint8_t *)buf)[i]);
    fprintf(fhandle, "\n");
    fclose(fhandle);
  }
}

void error_exit(char *str) {
  fprintf(stderr, "%s\n", str);
  exit(1);
}

int read_file(const char *file, uint8_t *buf, size_t len) {
  size_t c = 0;
  FILE *fhandle;
  fhandle = fopen(file, "r");
  if (fhandle == NULL) {
    fprintf(stderr, "Unable to open %s", file);
    exit(3);
  }
  for (size_t i = 0; i < len; ++i) {
    c += fscanf(fhandle, "%02hhx", &buf[i]);
  }
  fclose(fhandle);
  return c;
}

int read_stdin(uint8_t *buf, int len) {
  size_t c = 0;
  for (size_t i = 0; i < len; ++i) {
    c += fscanf(stdin, "%02hhx", &buf[i]);
  }
  return c;
}

int main(int argc, char **argv) {
  uint8_t private_key[N];
  proj public_key;
  proj shared_secret_key;

  //explicit_bzero(&private_key, sizeof(private_key));
  //explicit_bzero(&public_key, sizeof(private_key));
  //explicit_bzero(&shared_secret_key, sizeof(shared_secret_key));

  int csidh_validate = 0;
  int option = 0;
  size_t verbose = 0;
  size_t generation_mode = 0;
  size_t derivation_mode = 0;
  size_t error = 0;
  char *priv_key_file = NULL;
  char *pub_key_file = NULL;

  while ((option = getopt(argc, argv, "hvVgdp:s:")) != -1) {
    switch (option) {
    case 'V':
      fprintf(stderr, "csidh-p%i-util version: %f\n", BITS, VERSION);
      return 0;
    case 'h':
      fprintf(stderr, "csidh-p%i-util version: %f\n", BITS, VERSION);
      fprintf(stderr, "  -V: print version\n");
      fprintf(stderr, "  -verbose: increase verbosity\n");
      fprintf(stderr, "  -g: key generation mode\n");
      fprintf(stderr, "  -d: key derivation mode\n");
      fprintf(stderr, "  -p: public key file name\n");
      fprintf(stderr, "  -s: private key file name\n");
      return 0;
    case 'v':
      verbose += 1;
      break;
    case 'g':
      generation_mode = 1;
      if (derivation_mode) {
        error += 1;
      };
      break;
    case 'd':
      derivation_mode = 1;
      if (generation_mode) {
        error = 1;
      };
      break;
    case 'p':
      pub_key_file = optarg;
      if (verbose) {
        fprintf(stderr, "pub_key_file=%s\n", pub_key_file);
      };
      break;
    case 's':
      priv_key_file = optarg;
      if (verbose) {
        fprintf(stderr, "priv_key_file=%s\n", priv_key_file);
      };
      break;
    default:
      exit(1);
    }
  }

  if (error) {
    error_exit("Mutually exclusive options chosen");
  }

  if (generation_mode) {
    if (verbose) {
      fprintf(stderr, "Key generation mode\n");
    }
    random_key(private_key);
    csidh_validate = csidh(public_key, private_key, E);
    if (!csidh_validate && !validate(public_key)) {
      error_exit("csidh_validate: failed");
    }
    if (priv_key_file == NULL){
      pprint_sk(private_key);
    } else {
      save_file(priv_key_file, &private_key, sizeof(uint8_t [N]));
    }
    if (pub_key_file == NULL){
      pprint_pk(public_key);
    } else {
      save_file(pub_key_file, public_key, sizeof(proj));
    }
    return 0;
  }

  if (derivation_mode) {
    if (verbose) {
      fprintf(stderr, "DH mode\n");
    }

    if (sizeof(private_key) !=
        ((priv_key_file != NULL)
             ? read_file(priv_key_file, (uint8_t *)&private_key,
                         sizeof(private_key))
             : read_stdin((uint8_t *)&private_key, sizeof(private_key)))) {
      error_exit("Unable to read correct number of bytes for private key");
    }

    if (sizeof(public_key) !=
        ((pub_key_file != NULL)
             ? read_file(pub_key_file, (uint8_t *)public_key,
                         sizeof(proj))
             : read_stdin((uint8_t *)public_key, sizeof(proj)))) {
      error_exit("Unable to read correct number of bytes for public key");
    }

    if (!validate(public_key)) {
      error_exit("csidh_validate: failed");
    }
    csidh_validate = csidh(shared_secret_key, private_key, public_key);
    if (!csidh_validate) {
      error_exit("csidh_validate: failed");
    }
    fp shared_secret;
    fp_inv(shared_secret_key[1]);
    fp_mul(shared_secret, shared_secret_key[0], shared_secret_key[1]);

    if (verbose) {
      pprint_sk(private_key);
      pprint_pk(public_key);
    }
    pprint_ss(shared_secret);
    return 0;
  }

  return 1;
}
