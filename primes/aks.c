#include <gmp.h>
#include <mpfr.h>
#include <stdlib.h>

#define FALSE 0
#define TRUE 1

#define COMPOSITE 0
#define PRIME 1


void compute_logn(mpz_t rop, mpz_t n) {
  mpfr_t tmp;
  mpfr_init(tmp);
  mpfr_set_z(tmp, n, MPFR_RNDN);
  mpfr_log(tmp, tmp, MPFR_RNDN);
  mpfr_get_z(rop, tmp, MPFR_RNDN);
}

void compute_logn2(mpz_t rop, mpz_t n) {
  mpfr_t tmp;
  mpfr_init(tmp);
  
  mpfr_set_z(tmp, n, MPFR_RNDN);
  mpfr_log(tmp, tmp, MPFR_RNDA);
  mpfr_pow_ui(tmp, tmp, 2, MPFR_RNDA);
  mpfr_ceil(tmp, tmp);

  mpfr_get_z(rop, tmp, MPFR_RNDA);
}

void find_smallest_r(mpz_t r, mpz_t n) {
  mpz_t logn2, i, tmp;
  mpz_init(logn2);
  mpz_init(i);
  mpz_init(tmp);

  compute_logn2(logn2, n);
  mpz_set(r, logn2);

  int found_r = FALSE;
  while (!found_r) {
    found_r = TRUE;
    for (mpz_set_ui(i, 1); mpz_cmp(i, logn2) <= 0; mpz_add_ui(i, i, 1)) {
       mpz_powm(tmp, n, i, r);
       if (mpz_cmp_ui(tmp, 1) == 0) {
         found_r = FALSE;
         break;
       }
    }
    if (!found_r) {
      mpz_add_ui(r, r, 1);
    }
  }
}

int check_a_exists(mpz_t n, mpz_t r) {
  mpz_t a, gcd;
  mpz_init(a);
  mpz_init(gcd);

  for (mpz_set_ui(a, 1); mpz_cmp(a, r) <= 0; mpz_add_ui(a, a, 1)) {
    mpz_gcd(gcd, a, n);
    if (mpz_cmp_ui(gcd, 1) > 0 && mpz_cmp(gcd, n) < 0)
      return TRUE;
  }

  return FALSE;
}

void totient(mpz_t rop, mpz_t op) {
  mpz_t i, gcd;
  mpz_init(i);
  mpz_init(gcd);
  mpz_set_ui(rop, 0);

  for (mpz_set(i, op); mpz_cmp_ui(i, 0) != 0; mpz_sub_ui(i, i, 1)) {
    mpz_gcd(gcd, i, op);
    if (mpz_cmp_ui(gcd, 1) == 0)
      mpz_add_ui(rop, rop, 1);
  } 
}

void compute_upper_limit(mpz_t rop, mpz_t r, mpz_t n) {
  mpz_t tot, logn;
  mpz_init(tot);
  mpz_init(logn);

  totient(tot, r);
  gmp_printf("tot=%Zd\n", tot);
  mpz_sqrt(tot, tot);

  compute_logn(logn, n);
  mpz_mul(rop, tot, logn);
}

int check_poly(mpz_t n, mpz_t a, mpz_t r) {
  int i, terms, equality_holds;

  mpz_t tmp, loop;
  mpz_init(tmp);
  mpz_init(loop);

  terms = mpz_get_ui(r) + 1;
  mpz_t* poly = malloc(sizeof(mpz_t) * terms);

  for (i = 0; i < terms; i++)
    mpz_init(poly[i]);

  mpz_set(poly[0], a);
  mpz_set_ui(poly[1], 1);

  for (mpz_set_ui(loop, 2); mpz_cmp(loop, n) <= 0; mpz_add_ui(loop, loop, 1)) {
    mpz_set(tmp, poly[terms-1]);

    for (i = terms - 1; i > 0; i--) {
      mpz_mul(poly[i], poly[i], a);
      mpz_add(poly[i], poly[i], poly[i-1]);
      mpz_mod(poly[i], poly[i], n);
    }

    mpz_mul(poly[0], poly[0], a);
    mpz_add(poly[0], poly[0], poly[i-1]);
    mpz_mod(poly[0], poly[0], n);
  }

  equality_holds = TRUE;
  if (mpz_cmp(poly[0], a) != 0 || mpz_cmp_ui(poly[terms-1], 1) != 0) {
    equality_holds = FALSE;
  }
  else {
    for (i = 1; i < terms - 1; i++) {
      if (mpz_cmp_ui(poly[i], 0) != 0)
        equality_holds = FALSE;
        break;
    }
  }

  for (i = 0; i < terms; i++)
    mpz_clear(poly[i]);
  mpz_clear(tmp);
  mpz_clear(loop);

  return equality_holds;
}

int check_polys(mpz_t r, mpz_t n) {
  mpz_t a, lim;
  mpz_init(a);
  mpz_init(lim);

  gmp_printf("computing upper limit\n");
  compute_upper_limit(lim, r, n);
  gmp_printf("lim=%Zd\n", lim);
  for (mpz_set_ui(a, 1); mpz_cmp(a, lim) <= 0; mpz_add_ui(a, a, 1)) {
    if (!check_poly(n, a, r))
      return COMPOSITE;
  }
}

int is_prime(mpz_t n) {
  if (mpz_cmp_ui(n, 1) <= 0 || mpz_divisible_ui_p(n, 2))
    return COMPOSITE; 

  if (mpz_perfect_power_p(n))
    return COMPOSITE;

  mpz_t r;
  mpz_init(r);
  find_smallest_r(r, n);

  gmp_printf("r=%Zd\n", r);

  if (check_a_exists(n, r))
   return COMPOSITE;

  gmp_printf("a does not exist\n");

  if (mpz_cmp(n, r) <= 0)
    return PRIME;

  gmp_printf("checking polynomial equation\n");

  if (check_polys(r, n))
    return COMPOSITE;

  return PRIME; 
}

int main(int argc, char** argv) {
  if (argc < 2) {
    gmp_printf("usage: aks num\n");
    return 1;
  }

  mpz_t n;
  mpz_init(n);
  mpz_set_str(n, argv[1], 0);
  int prime = is_prime(n);
  gmp_printf("%d\n", prime);

  return 0;
}

