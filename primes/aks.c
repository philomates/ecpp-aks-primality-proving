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
  mpfr_clear(tmp);  
}

void compute_logn2(mpz_t rop, mpz_t n) {
  mpfr_t tmp;
  mpfr_init(tmp);
  
  mpfr_set_z(tmp, n, MPFR_RNDN);
  mpfr_log(tmp, tmp, MPFR_RNDA);
  mpfr_pow_ui(tmp, tmp, 2, MPFR_RNDA);
  mpfr_ceil(tmp, tmp);

  mpfr_get_z(rop, tmp, MPFR_RNDA);
  mpfr_clear(tmp);
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

  mpz_clear(logn2);
  mpz_clear(i);
  mpz_clear(tmp);
}

int check_a_exists(mpz_t n, mpz_t r) {
  mpz_t a, gcd;
  mpz_init(a);
  mpz_init(gcd);

  int exists = FALSE;

  for (mpz_set_ui(a, 1); mpz_cmp(a, r) <= 0; mpz_add_ui(a, a, 1)) {
    mpz_gcd(gcd, a, n);
    if (mpz_cmp_ui(gcd, 1) > 0 && mpz_cmp(gcd, n) < 0) {
      exists = TRUE;
      break;
    }
  }

  mpz_clear(a);
  mpz_clear(gcd);

  return exists;
}

void totient(mpz_t rop, mpz_t op) {
  mpz_t i, gcd;
  mpz_init(i);
  mpz_init(gcd);
  mpz_set_ui(rop, 0);

  for (mpz_set(i, op); mpz_cmp_ui(i, 0) != 0; mpz_sub_ui(i, i, 1)) {
    mpz_gcd(gcd, i, op);
    if (mpz_cmp_ui(gcd, 1) == 0) {
      mpz_add_ui(rop, rop, 1);
    }
  }

  mpz_clear(i);
  mpz_clear(gcd);
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

  mpz_clear(tot);
  mpz_clear(logn);
}

void polymul(mpz_t* rop, mpz_t* op1, unsigned int len1, mpz_t* op2, unsigned int len2, mpz_t n) {
  int i, j, t;

  for (i = 0; i < len1; i++) {
    for (j = 0; j < len2; j++) {
      t = (i + j) % len1;
      mpz_addmul(rop[t], op1[i], op2[j]);
      mpz_mod(rop[t], rop[t], n);
    }
  }
}

mpz_t* init_poly(unsigned int terms) {
  int i;
  mpz_t* poly = malloc(sizeof(mpz_t) * terms);
  for (i = 0; i < terms; i++) {
    mpz_init(poly[i]);
  }
  return poly;
}

void clear_poly(mpz_t* poly, unsigned int terms) {
  int i;
  for (i = 0; i < terms; i++) {
    mpz_clear(poly[i]);
  }
  free(poly);
}

int check_poly(mpz_t n, mpz_t a, mpz_t r) {
  unsigned int i, terms, equality_holds;

  mpz_t tmp, neg_a, loop;
  mpz_init(tmp);
  mpz_init(neg_a);
  mpz_init(loop);

  terms = mpz_get_ui(r) + 1;
  mpz_t* poly = init_poly(terms);
  mpz_t* ptmp = init_poly(terms);
  mpz_t* stmp;

  mpz_mul_ui(neg_a, a, -1);

  mpz_set(poly[0], neg_a);
  mpz_set_ui(poly[1], 1);

  for (mpz_set_ui(loop, 2); mpz_cmp(loop, n) <= 0; mpz_mul(loop, loop, loop)) {
    polymul(ptmp, poly, terms, poly, terms, n);
    stmp = poly;
    poly = ptmp;
    ptmp = stmp;
  }

  mpz_t* xMinusA = init_poly(2);
  mpz_set(ptmp[0], neg_a);
  mpz_set_ui(ptmp[1], 1);
  for (; mpz_cmp(loop, n) <= 0; mpz_add_ui(loop, loop, 1)) {
    polymul(ptmp, poly, terms, xMinusA, 2, n);
    stmp = poly;
    poly = ptmp;
    ptmp = stmp;
  }
  clear_poly(xMinusA, 2);

  equality_holds = TRUE;
  if (mpz_cmp(poly[0], neg_a) != 0 || mpz_cmp_ui(poly[terms - 1], 1) != 0) {
    equality_holds = FALSE;
  }
  else {
    for (i = 1; i < terms - 1; i++) {
      if (mpz_cmp_ui(poly[i], 0) != 0) {
        equality_holds = FALSE;
        break;
      }
    }
  }

  clear_poly(poly, terms);
  clear_poly(ptmp, terms);

  mpz_clear(tmp);
  mpz_clear(neg_a);
  mpz_clear(loop);

  return equality_holds;
}

int check_polys(mpz_t r, mpz_t n) {
  mpz_t a, lim;
  mpz_init(a);
  mpz_init(lim);

  int status = PRIME;
  gmp_printf("computing upper limit\n");
  compute_upper_limit(lim, r, n);
  gmp_printf("lim=%Zd\n", lim);
  for (mpz_set_ui(a, 1); mpz_cmp(a, lim) <= 0; mpz_add_ui(a, a, 1)) {
    if (!check_poly(n, a, r)) {
      status = COMPOSITE;
      break;
    }
  }

  mpz_clear(a);
  mpz_clear(lim);

  return status;
}

int is_prime(mpz_t n) {
  if (mpz_cmp_ui(n, 2) == 0) {
    return PRIME;
  }

  if (mpz_cmp_ui(n, 1) <= 0 || mpz_divisible_ui_p(n, 2)) {
    return COMPOSITE;
  }

  if (mpz_perfect_power_p(n)) {
    return COMPOSITE;
  }

  mpz_t r;
  mpz_init(r);
  find_smallest_r(r, n);

  gmp_printf("r=%Zd\n", r);

  if (check_a_exists(n, r)) {
    mpz_clear(r);
    return COMPOSITE;
  }

  gmp_printf("a does not exist\n");

  if (mpz_cmp(n, r) <= 0) {
    mpz_clear(r);
    return PRIME;
  }

  gmp_printf("checking polynomial equation\n");

  if (check_polys(r, n)) {
    mpz_clear(r);
    return COMPOSITE;
  }

  mpz_clear(r);
  return PRIME; 
}

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

