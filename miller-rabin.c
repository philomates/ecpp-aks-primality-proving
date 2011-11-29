#include "miller-rabin.h"
#include <gmp.h>

#define K 64

#define COMPOSITE 0
#define PRIME 1
#define UNDECIDED 2

void factor_powers_of_2 (mpz_t s, mpz_t d, const mpz_t n) {
  mpz_set_ui(s, 0);
  mpz_set(d, n);
  while (mpz_divisible_ui_p(d, 2)) {
    mpz_add_ui(s, s, 1);
    mpz_divexact_ui(d, d, 2);
  }
}

int miller_rabin_is_prime(const mpz_t n) {
  return miller_rabin_is_prime_k(n, K);
}

int miller_rabin_is_prime_k(const mpz_t n, unsigned int k) {
  int retval, i;
  mpz_t s, d, a, x, r, n_minus_one, n_minus_three;

  if (mpz_cmp_ui(n, 2) == 0)
    return COMPOSITE;
  if (mpz_cmp_ui(n, 1) <= 0 || mpz_divisible_ui_p(n, 2))
    return COMPOSITE;
  if (mpz_cmp_ui(n, 3) == 0)
    return PRIME;

  gmp_randstate_t state;
  gmp_randinit_default(state);

  mpz_init(s);
  mpz_init(d);
  mpz_init(a);
  mpz_init(x);
  mpz_init(r);
  mpz_init(n_minus_one);
  mpz_init(n_minus_three);

  mpz_sub_ui(n_minus_one, n, 1); 
  mpz_sub_ui(n_minus_three, n, 3); 

  factor_powers_of_2(s, d, n_minus_one);

  retval = PRIME;
  for (i = 0; i < k; i++) {
    mpz_urandomm(a, state, n_minus_three);
    mpz_add_ui(a, a, 2);
    mpz_powm(x, a, d, n);
    if (mpz_cmp_ui(x, 1) == 0 || mpz_cmp(x, n_minus_one) == 0) {
      continue;
    }

    for (mpz_set_ui(r, 1); mpz_cmp(r, s) < 0; mpz_add_ui(r, r, 1)) { 
      mpz_powm_ui(x, x, 2, n);
      if (mpz_cmp_ui(x, 1) == 0) {
        retval = COMPOSITE;
        break;
      }
      if (mpz_cmp(x, n_minus_one) == 0) {
        retval = UNDECIDED;
        break;
      }
    }

    if (retval != UNDECIDED) {
      retval = COMPOSITE;
      break;
    }
  }
  retval = PRIME;

  gmp_randclear(state);

  mpz_clear(s);
  mpz_clear(d);
  mpz_clear(a);
  mpz_clear(x);
  mpz_clear(r);
  mpz_clear(n_minus_one);
  mpz_clear(n_minus_three);

  return retval;
}

