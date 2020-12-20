#include <getopt.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>

#include "csidh_util.h"

void normalize_public_key(proj public_key, fp *out) {
    /* Compress the public key from x,y (256 bits) to x/y,NULL (128 bits) */
    // public_key has entries x,y or public_key[0] and public_key[1]
    fp n_public_key;
    fp_inv(public_key[1]); // fp_inv becomes x,1/y or public_key[0] = x and public_key[1] = 1/y
    fp_mul(n_public_key, public_key[0], public_key[1]); // normalized_public_key is x*1/y or x/y aka public_key[0]/public_key[1]
    /* Convert to Montgomery form */
    fp_mul(n_public_key, n_public_key, E[1]); // x/y * 4
    fp_sub(n_public_key, n_public_key, E[0]); // x/y - 2
    memcpy(out, n_public_key, sizeof(fp));
}

// slightly modified csidh from main/csidh.c
static uint8_t csidh(proj out, const uint8_t sk[], const proj in)
{
  if (!validate(in)) {
    return 0;
  };

  action_evaluation(out, sk, in);
  return 1;
};

void pprint_ss(uint64_t *x)
{
    /* for p512 we print 8 64bit little endian values as hex. */
    int ceiling = 7;
    int i;
    for(i=ceiling; i >= 0; --i){
        printf("%.16" PRIX64 "", x[i]);
    }
    printf("\n");
}

void pprint_pk(void *x) {
  for (size_t i = 0; i < (sizeof(proj)/2); ++i) {
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
      fprintf(stderr, "Unable to open %s\n", file);
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
    fprintf(stderr, "Unable to open %s\n", file);
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
  proj expanded_public_key;
  proj public_key;
  proj shared_secret_key;

  explicit_bzero(&private_key, sizeof(private_key));
  explicit_bzero(&public_key, sizeof(proj));
  explicit_bzero(&expanded_public_key, sizeof(proj));
  explicit_bzero(&shared_secret_key, sizeof(proj));

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
    if (!csidh_validate) {
      error_exit("csidh_validate: failed");
    }
    if (priv_key_file == NULL){
      pprint_sk(private_key);
    } else {
      if (verbose) { pprint_sk(private_key); }
      save_file(priv_key_file, &private_key, sizeof(uint8_t [N]));
    }

    /* Normalize (aka compress) the public key in size from 256 bits to 128 bits. */
    fp normalized_public_key;
    explicit_bzero(&normalized_public_key, sizeof(fp));
    normalize_public_key(public_key, &normalized_public_key);
    if (pub_key_file == NULL){
      pprint_pk(normalized_public_key);
    } else {
      if (verbose) { pprint_pk(normalized_public_key);}
      save_file(pub_key_file, normalized_public_key, (sizeof(proj)/2));
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

    if ((sizeof(public_key)/2) !=
        ((pub_key_file != NULL)
             ? read_file(pub_key_file, (uint8_t *)public_key,
                         (sizeof(proj)/2))
             : read_stdin((uint8_t *)public_key, (sizeof(proj)/2)))) {
      error_exit("Unable to read correct number of bytes for public key");
    }

    /* Expand the normalized key and convert from Montgomery to Edwards form. */
    /*
     * If one wanted to use Edwards rather than Montgomery, the three lines
     * below could be replaced with the following two lines:
     *
     *  memcpy(expanded_public_key[1], R_mod_p, (sizeof(expanded_public_key[1])));
     *  memcpy(expanded_public_key[0], public_key, (sizeof(expanded_public_key[0])));
     *
     * */
    memcpy(expanded_public_key[1], E[1], (sizeof(expanded_public_key[1]))); // E[1] ==2
    memcpy(expanded_public_key[0], public_key, (sizeof(expanded_public_key[0]))); // Original value from user in Montgomery form
    fp_add(expanded_public_key[0], expanded_public_key[0], E[0]); // E[0] == 4

    /* Operate on the expanded public key. */
    csidh_validate = csidh(shared_secret_key, private_key, expanded_public_key);
    if (!csidh_validate) {
      error_exit("csidh_validate: failed");
    }

    /* Normalize our shared secret. */
    fp shared_secret;
    fp_inv(shared_secret_key[1]);
    fp_mul(shared_secret, shared_secret_key[0], shared_secret_key[1]);
    /* Convert from Edwards to Montgomery. */
    fp_mul(shared_secret, shared_secret, E[1]); // x/y * 4
    fp_sub(shared_secret, shared_secret, E[0]); // x/y - 2

    if (verbose) {
      pprint_sk(private_key);
      pprint_pk(expanded_public_key);
    }
    /* Output shared secret in Little Endian to match the other known CSIDH
     * implementations. */
    pprint_ss(shared_secret);
    return 0;
  }

  return 1;
}
