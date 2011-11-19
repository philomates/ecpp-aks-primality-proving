#include <gmp.h>

int is_prime(const mpz_t n);

int main (int argc, char** argv) {
  mpz_t n;
  mpz_init(n);
  gmp_scanf("%Zd", &n);
  int prime = is_prime(n);
  gmp_printf("%d\n", prime);
  mpz_clear(n);
  return 0;
}

