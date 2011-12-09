#include "aks.h"

#include <gmp.h>
#include <mpfr.h>
#include <stdlib.h>

#define FALSE 0
#define TRUE 1

#define COMPOSITE 0
#define PRIME 1

int aks_debug = 0;

/**
 * Wrapper function to find the log of a number of type mpz_t.
 */
void compute_logn(mpz_t rop, mpz_t n) {
  mpfr_t tmp;
  mpfr_init(tmp);
  mpfr_set_z(tmp, n, MPFR_RNDN);
  mpfr_log(tmp, tmp, MPFR_RNDN);
  mpfr_get_z(rop, tmp, MPFR_RNDN);
  mpfr_clear(tmp);
}

/**
 * Wrapper function to find the square of the log of a number of type mpz_t.
 */
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

/**
 * Finds the smallest r such that order of the a modula r which is the 
 * smallest number k such that n ^ k = 1 (mod r) is greater than log(n) ^ 2.
 */
void find_smallest_r(mpz_t r, mpz_t n) {
  mpz_t logn2, k, tmp;
  mpz_init(logn2);
  mpz_init(k);
  mpz_init(tmp);

  // Compute log(n) ^ 2 in order to do the comparisons
  compute_logn2(logn2, n);
  // R must be at least log(n) ^ 2
  mpz_set(r, logn2);

  int found_r = FALSE;
  while (!found_r) {
    found_r = TRUE;
    // Check several values of k from 1 up to log(n) ^ 2 to find one that satisfies the equality
    for (mpz_set_ui(k, 1); mpz_cmp(k, logn2) <= 0; mpz_add_ui(k, k, 1)) {
       // Compute n ^ k % r
       mpz_powm(tmp, n, k, r);
       // If it is not equal to 1 than the equality n ^ k = 1 (mod r) does not hold
       // and we must find test a different value of k
       if (mpz_cmp_ui(tmp, 1) == 0) {
         found_r = FALSE;
         break;
       }
    }
    // All possible values of k were checked so we must start looking for a new r
    if (!found_r) {
      mpz_add_ui(r, r, 1);
    }
  }

  mpz_clear(logn2);
  mpz_clear(k);
  mpz_clear(tmp);
}

/**
 * Return true if there exists an a such that 1 < gcd(a, n) < n for some a <= r. 
 */
int check_a_exists(mpz_t n, mpz_t r) {
  mpz_t a, gcd;
  mpz_init(a);
  mpz_init(gcd);

  int exists = FALSE;

  // Simply iterate for values of a from 1 to r and see if equations hold for the gcd of a and n
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

/**
 * Return the totient of op, which is the count of numbers less than op which are coprime to op.
 */
void totient(mpz_t rop, mpz_t op) {
  mpz_t i, gcd;
  mpz_init(i);
  mpz_init(gcd);
  mpz_set_ui(rop, 0);

  // Simply iterate through all values from op to 1 and see if the gcd of that number and op is 1.
  // If it is then it is coprime and added to the totient count.
  for (mpz_set(i, op); mpz_cmp_ui(i, 0) != 0; mpz_sub_ui(i, i, 1)) {
    mpz_gcd(gcd, i, op);
    if (mpz_cmp_ui(gcd, 1) == 0) {
      mpz_add_ui(rop, rop, 1);
    }
  }

  mpz_clear(i);
  mpz_clear(gcd);
}

/**
 * Returns sqrt(totient(r)) * log(n) which is used by step 5.
 */
void compute_upper_limit(mpz_t rop, mpz_t r, mpz_t n) {
  mpz_t tot, logn;
  mpz_init(tot);
  mpz_init(logn);

  totient(tot, r);
  if (aks_debug) gmp_printf("tot=%Zd\n", tot);
  mpz_sqrt(tot, tot);

  compute_logn(logn, n);
  mpz_mul(rop, tot, logn);

  mpz_clear(tot);
  mpz_clear(logn);
}

/**
 * Multiplies two polynomials with appropriate mod where their coefficients are indexed into the array. 
 */
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

/**
 * Allocates an array where each element represents a coeffecient of the polynomial.
 */
mpz_t* init_poly(unsigned int terms) {
  int i;
  mpz_t* poly = (mpz_t*) malloc(sizeof(mpz_t) * terms);
  for (i = 0; i < terms; i++) {
    mpz_init(poly[i]);
  }
  return poly;
}

/**
 * Frees the array and clears each element in the array.
 */
void clear_poly(mpz_t* poly, unsigned int terms) {
  int i;
  for (i = 0; i < terms; i++) {
    mpz_clear(poly[i]);
  }
  free(poly);
}

/**
 * Test if (X + a) ^ n != X ^ n + a (mod X ^ r - 1,n)
 */
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

/**
 * Run step 5 of the AKS algorithm.
 */
int check_polys(mpz_t r, mpz_t n) {
  mpz_t a, lim;
  mpz_init(a);
  mpz_init(lim);

  int status = PRIME;
  if (aks_debug) gmp_printf("computing upper limit\n");
  compute_upper_limit(lim, r, n);
  if (aks_debug) gmp_printf("lim=%Zd\n", lim);
  // For values of a from 1 to sqrt(totient(r)) * log(n)
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

int aks_is_prime(mpz_t n) {

  // Peform simple checks before running the AKS algorithm
  if (mpz_cmp_ui(n, 2) == 0) {
    return PRIME;
  }

  if (mpz_cmp_ui(n, 1) <= 0 || mpz_divisible_ui_p(n, 2)) {
    return COMPOSITE;
  }

  // Step 1: Check if n is a perfect power, meaning n = a ^ b where a is a natural number and b > 1
  if (mpz_perfect_power_p(n)) {
    return COMPOSITE;
  }

  // Step 2: Find the smallest r such that or(n) > log(n) ^ 2
  mpz_t r;
  mpz_init(r);
  find_smallest_r(r, n);

  if (aks_debug) gmp_printf("r=%Zd\n", r);

  // Step 3: Check if there exists an a <= r such that 1 < (a,n) < n
  if (check_a_exists(n, r)) {
    mpz_clear(r);
    return COMPOSITE;
  }

  if (aks_debug) gmp_printf("a does not exist\n");

  // Step 4: Check if n <= r
  if (mpz_cmp(n, r) <= 0) {
    mpz_clear(r);
    return PRIME;
  }

  if (aks_debug) gmp_printf("checking polynomial equation\n");

  // Step 5
  if (check_polys(r, n)) {
    mpz_clear(r);
    return COMPOSITE;
  }

  mpz_clear(r);

  // Step 6
  return PRIME;
}

