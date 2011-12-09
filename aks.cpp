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
  mpfr_log2(tmp, tmp, MPFR_RNDN);
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
  mpfr_log2(tmp, tmp, MPFR_RNDA);
  mpfr_pow_ui(tmp, tmp, 2, MPFR_RNDA);
  mpfr_ceil(tmp, tmp);

  mpfr_get_z(rop, tmp, MPFR_RNDA);
  mpfr_clear(tmp);
}

/**
 * Finds the smallest r such that order of the a modula r which is the 
 * smallest number k such that n ^ k = 1 (mod r) is greater than log(n) ^ 2.
 * 
 * Will also return composite if the number is guaranteed to be composite.
 */
int find_smallest_r(mpz_t r, mpz_t n) {
  int retval = PRIME;
  int found_r;

  mpz_t logn2, k, tmp;
  mpz_init(logn2);
  mpz_init(k);
  mpz_init(tmp);

  // Compute log(n) ^ 2 in order to know when r is found
  compute_logn2(logn2, n);

  // The value r <= ceil(log(n) ^ 5), but can go up to n for n <= 5,690,034
  for (mpz_set_ui(r, 2); mpz_cmp(r, n) < 0; mpz_add_ui(r, r, 1)) {

    if (mpz_divisible_p(n, r)) {
      retval = COMPOSITE;
      break;
    }

    mpz_gcd(tmp, n, r);
    if (mpz_cmp_ui(tmp, 1) != 0) {
      continue;
    }

    found_r = TRUE;
 
    // Check values of k from 1 up to log(n) ^ 2 to see if any value satisfies the equality
    for (mpz_set_ui(k, 1); mpz_cmp(k, logn2) <= 0; mpz_add_ui(k, k, 1)) {
 
       // Compute n ^ k % r
       mpz_powm(tmp, n, k, r);

       if (mpz_cmp_ui(tmp, 1) == 0) {
         found_r = FALSE;
         break;
       }
    }

    // k is greater than log(n) ^ 2 so we have found r
    if (found_r == TRUE) {
      break;
    }
  }

  mpz_clear(logn2);
  mpz_clear(k);
  mpz_clear(tmp);

  return retval;
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
 * Return the totient of op, which is the count of numbers less than or equal to op which are coprime to op.
 */
void totient(mpz_t rop, mpz_t op) {
  mpz_t i, gcd;

  mpz_init(i);
  mpz_init(gcd);

  mpz_set_ui(rop, 0);

  for (mpz_sub_ui(i, op, 1); mpz_cmp_ui(i, 0) != 0; mpz_sub_ui(i, i, 1)) {
    mpz_gcd(gcd, i, op);
    if (mpz_cmp_ui(gcd, 1) == 0) {
      mpz_add_ui(rop, rop, 1);
    }
  }

  mpz_clear(i);
  mpz_clear(gcd);
}

/**
 * Returns sqrt(r) * log(n) which is used by step 5.
 * Uses the Lenstra and Pomerance improvements
 */
void compute_upper_limit(mpz_t rop, mpz_t r, mpz_t n) {
  mpz_t logn;
  mpz_init(logn);

  totient(rop, r);
  mpz_sqrt(rop, rop);
  compute_logn(logn, n);
  mpz_mul(rop, rop, logn);

  mpz_clear(logn);
}

/**
 * Prints a coefficient array representing a polynomial.
 */ 
void polyprint(mpz_t* poly, unsigned int len) {
  unsigned int i;
  for (i = 0; i < len; i++) {
    gmp_printf("%Zd ", poly[i]);
  }
  gmp_printf("\n");
}

/**
 * Copies the op polynomial coefficient array into rop.
 */
void polycopy(mpz_t* rop, mpz_t* op, unsigned int len) {
  unsigned int i;
  for (i = 0; i < len; i++) {
    mpz_set(rop[i], op[i]);
  }
}

/**
 * Multiplies two polynomials (mod x^r - 1, n) where the ops are coefficient arrays.
 */
void polymul(mpz_t* rop, mpz_t* op1, mpz_t* op2, unsigned int len, mpz_t n, mpz_t* tmp) {
  int i, j, t;

  for (i = 0; i < len; i++) {
    mpz_set_ui(tmp[i], 0);
  }

  for (i = 0; i < len; i++) {
    for (j = 0; j < len; j++) {
      t = (i + j) % len;
      mpz_addmul(tmp[t], op1[i], op2[j]);
      mpz_mod(tmp[t], tmp[t], n);
    }
  }

  polycopy(rop, tmp, len);  
}

/**
 * Raise the given polynomial to the given power (mod x^r - 1, n). This is done
 * by repeated squaring.
 */
void polypow(mpz_t* rop, mpz_t* op, unsigned int len, mpz_t n, mpz_t* tmp, mpz_t* tmp2) {
  mpz_t s, i, remainder;
  mpz_init(s);
  mpz_init(i);
  mpz_init(remainder);

  mpz_set_ui(s, 0);
  polycopy(rop, op, len);

  while (mpz_cmp(s, n) < 0) {
    mpz_sub(remainder, n, s);
    polycopy(tmp, op, len);

    if (mpz_cmp_ui(remainder, 1) == 0) {
      polymul(rop, rop, op, len, n, tmp);
      break;
    }

    for (mpz_set_ui(i, 2); mpz_cmp(i, remainder) <= 0; mpz_mul_ui(i, i, 2)) {
      polymul(tmp, tmp, tmp, len, n, tmp2);
    }

    if (mpz_cmp_ui(s, 0) == 0) {
      polycopy(rop, tmp, len);
    }
    else {
      polymul(rop, rop, tmp, len, n, tmp2);
    }

    mpz_divexact_ui(i, i, 2);
    mpz_add(s, s, i);
  }

  mpz_clear(s);
  mpz_clear(i);
  mpz_clear(remainder);
}

/**
 * Allocates an array where each element represents a coeffecient of the polynomial.
 */
mpz_t* init_poly(unsigned int terms) {
  int i;
  mpz_t* poly = (mpz_t*) malloc(sizeof(mpz_t) * terms);
  for (i = 0; i < terms; i++) {
    mpz_init(poly[i]);
    mpz_set_ui(poly[i], 0);
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
int check_poly(mpz_t n, mpz_t a, mpz_t r, mpz_t* poly, mpz_t* xpa, mpz_t* tmp1, mpz_t* tmp2) {
  unsigned int i, terms, retval;

  terms = mpz_get_ui(r);

  mpz_set(xpa[0], a);
  mpz_set_ui(xpa[1], 1);

  polypow(poly, xpa, terms, n, tmp1, tmp2);

  retval = PRIME;

  mpz_t nmodr;
  mpz_init(nmodr);
  mpz_mod(nmodr, n, r);

  unsigned int index = mpz_get_ui(nmodr);

  if (mpz_cmp(poly[0], a) != 0 || mpz_cmp_ui(poly[index], 1) != 0) {
    retval = COMPOSITE;
  }
  if (retval == PRIME) {
    for (i = 1; i < index; i++) {
      if (mpz_cmp_ui(poly[i], 0) != 0) {
        retval = COMPOSITE;
        break;
      }
    }
    if (retval == PRIME) {
      for (i = index + 1; i < terms; i++) {
        if (mpz_cmp_ui(poly[i], 0) != 0) {
          retval = COMPOSITE;
          break;
        }
      }
    }
  }

  mpz_clear(nmodr);

  if (aks_debug) gmp_printf("check_poly returning %d\n", retval); 
  return retval;
}

/**
 * Run step 5 of the AKS algorithm.
 */
int check_polys(mpz_t r, mpz_t n) {
  int retval = PRIME;

  mpz_t a, lim;
  mpz_init(a);
  mpz_init(lim);

  unsigned int terms = mpz_get_ui(r);

  mpz_t* pol1 = init_poly(terms);
  mpz_t* pol2 = init_poly(terms);
  mpz_t* tmp1 = init_poly(terms);
  mpz_t* tmp2 = init_poly(terms);

  if (aks_debug) gmp_printf("computing upper limit\n");

  compute_upper_limit(lim, r, n);

  if (aks_debug) gmp_printf("lim=%Zd\n", lim);

  for (mpz_set_ui(a, 1); mpz_cmp(a, lim) <= 0; mpz_add_ui(a, a, 1)) {

    if (aks_debug) gmp_printf("checking a=%Zd\n", a);

    if (check_poly(n, a, r, pol1, pol2, tmp1, tmp2) == COMPOSITE) {
      if (aks_debug) gmp_printf("proven composite\n");
      retval = COMPOSITE;
      break;
    }
  }

  clear_poly(pol1, terms);
  clear_poly(pol2, terms);
  clear_poly(tmp1, terms);
  clear_poly(tmp2, terms);

  mpz_clear(a);
  mpz_clear(lim);

  return retval;
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
  if (find_smallest_r(r, n) == COMPOSITE) {
    mpz_clear(r);
    return COMPOSITE;
  }

  if (aks_debug) gmp_printf("r=%Zd\n", r);

  // Step 3: Check if there exists an a <= r such that 1 < (a,n) < n
  if (check_a_exists(n, r)) {
    mpz_clear(r);
    return COMPOSITE;
  }

  if (aks_debug) gmp_printf("a does not exist\n");

  // Step 4: Check if n <= r, is only relevant for n <= 5,690,034
  if (mpz_cmp(n, r) <= 0) {
    mpz_clear(r);
    return PRIME;
  }

  if (aks_debug) gmp_printf("checking polynomial equation\n");

  // Step 5
  if (check_polys(r, n) == COMPOSITE) {
    mpz_clear(r);
    return COMPOSITE;
  }

  mpz_clear(r);

  // Step 6
  return PRIME;
}

