#include "miller-rabin.h"
#include <gmp.h>

int main (int argc, char** argv) {
  mpz_t n;
  mpz_init(n);
  gmp_scanf("%Zd", &n);
  int prime = miller_rabin_is_prime(n);
  gmp_printf("%d\n", prime);
  mpz_clear(n);
  return 0;
}

