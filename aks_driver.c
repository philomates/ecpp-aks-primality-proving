#include "aks.h"

#include <gmp.h>
#include <mpfr.h>
#include <stdlib.h>

int main(int argc, char** argv) {
  mpz_t n;
  mpz_init(n);
  gmp_scanf("%Zd", &n);
  int prime = is_prime(n);
  gmp_printf("%d\n", prime);
  mpz_clear(n);
  mpfr_free_cache();
  return 0;
}
