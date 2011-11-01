#include <gmp.h>

#define k 64

#define FALSE 0
#define TRUE  1

#define COMPOSITE 0
#define PRIME 1


void factor_powers_of_2 (mpz_t s, mpz_t d, const mpz_t n) {
  mpz_set_ui(s, 0);
  mpz_set(d, n);
  while (mpz_divisible_ui_p(d, 2)) {
    mpz_add_ui(s, s, 1);
    mpz_divexact_ui(d, d, 2);
  }
}


int is_prime(const mpz_t n) {
  int i, should_continue;
  mpz_t s, d, a, x, r, n_minus_one, n_minus_three;

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

  for (i = 0; i < k; i++) {
    mpz_urandomm(a, state, n_minus_three);
    mpz_add_ui(a, a, 2);
    mpz_powm(x, a, d, n);
    if (mpz_cmp_ui(x, 1) == 0 || mpz_cmp(x, n_minus_one) == 0)
      continue;

    should_continue = FALSE; 
    for (mpz_set_ui(r, 1); mpz_cmp(r, s) < 0; mpz_add_ui(r, r, 1)) { 
      mpz_powm_ui(x, x, 2, n);
      if (mpz_cmp_ui(x, 1) == 0) {
        return FALSE;
      }
      if (mpz_cmp(x, n_minus_one) == 0) {
        should_continue = TRUE;
        break;
      }
    }
    if (should_continue)
      continue;
    return FALSE;
  }
  return TRUE;
}

int main (int argc, char** argv) {
  if (argc < 2) {
    gmp_printf("usage: miller-rabin num\n");
    return 1;
  }
 
  mpz_t n;
  mpz_init_set_str (n, argv[1], 0);
  int prime = is_prime(n);
  gmp_printf("%d\n", prime);
  return 0;
}

