// <PROJECT>_<PATH>_<FILE>_H_
#ifndef PRIMES_AKS_H_
#define PRIMES_AKS_H_

#include <gmp.h>

void compute_logn(mpz_t, mpz_t);
void compute_logn2(mpz_t, mpz_t);
void find_smallest_r(mpz_t, mpz_t);
int check_a_exists(mpz_t, mpz_t);
void totient(mpz_t, mpz_t op);
void compute_upper_limit(mpz_t, mpz_t, mpz_t);
void polymul(mpz_t*, mpz_t*, unsigned int, mpz_t*, unsigned int, mpz_t);
void clear_poly(mpz_t*, unsigned int);
int check_poly(mpz_t, mpz_t, mpz_t);
int check_polys(mpz_t, mpz_t);
int aks_is_prime(mpz_t);

#endif // PRIMES_AKS_H_
