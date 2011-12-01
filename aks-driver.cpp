#include "aks.h"

#include <gmp.h>
#include <mpfr.h>
#include <stdlib.h>
#include <string.h>

/**
 * Debug statements will be used if the -d option is used.
 */
int main(int argc, char** argv) {
  if (argc > 1 && strcmp(argv[1], "-d") == 0) {
    aks_debug = 1;
  }

  mpz_t n;
  mpz_init(n);

  while (!feof(stdin)) {
    gmp_scanf("%Zd", &n);
    int prime = aks_is_prime(n);
    gmp_printf("%d\n", prime);
  }

  mpz_clear(n);
  mpfr_free_cache();
  return 0;
}
