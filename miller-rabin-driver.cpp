#include "miller-rabin.h"
#include <gmp.h>
#include <string.h>
#include <stdio.h>

#define FALSE 0
#define TRUE 1

/**
 * If the -gmp option is specied then the GMP implementation is used.
 */
int main (int argc, char** argv) {
  int use_gmp_solution = FALSE;
  if (argc > 1 && strcmp(argv[1], "-gmp") == 0) {
    use_gmp_solution = TRUE;
  }

  mpz_t n;
  mpz_init(n);

  while (!feof(stdin)) {
    gmp_scanf("%Zd", &n);
    if(feof(stdin))
      break;
    int prime = use_gmp_solution ? mpz_probab_prime_p(n, DEFAULT_K) : miller_rabin_is_prime(n, DEFAULT_K);
    gmp_printf("%d\n", prime);
  }
  mpz_clear(n);
  return 0;
}

