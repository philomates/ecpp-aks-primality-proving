#include "miller-rabin.h"
#include <gmp.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>

#define FALSE 0
#define TRUE 1

/**
 * Returns a number that is extremely likely to be prime.
 * If the -n option is specified the number will not be prime.
 */
int main(int argc, char** argv) {
  if (argc < 3) {
    gmp_printf("usage: gprime min maxi [-n]\n");
    return 1;
  }

  int generate_prime = TRUE;
  if (argc > 3 && strcmp(argv[3], "-n") == 0) {
    generate_prime = FALSE;
  } 

  gmp_randstate_t state;
  mpz_t n, min, max, count;

  gmp_randinit_default(state); 
  struct timeval t;
  gettimeofday(&t, NULL);
  gmp_randseed_ui(state, t.tv_usec); 

  mpz_init(n);
  mpz_init(min);
  mpz_init(max);
  mpz_init(count);

  mpz_set_str(min, argv[1], 0);
  mpz_set_str(max, argv[2], 0);
  mpz_sub(max, max, min);
  mpz_add_ui(max, max, 1);

  int is_prime;
  for (mpz_set_ui(count, 0); mpz_cmp(count, max) < 0; mpz_add_ui(count, count, 1)) {
    mpz_urandomm(n, state, max);
    mpz_add(n, n, min);
    is_prime = miller_rabin_is_prime(n, DEFAULT_K);
    if ((is_prime && generate_prime) || (!is_prime && !generate_prime)) {
      break;
    }
    else {
      mpz_set_ui(n, 0);
    }
  }

  gmp_printf("%Zd\n", n);

  gmp_randclear(state);  

  mpz_clear(n);
  mpz_clear(min);
  mpz_clear(max);
  mpz_clear(count);

  return 0;
}

