// <PROJECT>_<PATH>_<FILE>_H_
#ifndef PRIMES_MILLER_RABIN_H_
#define PRIMES_MILLER_RABIN_H_

#include <gmp.h>

int miller_rabin_is_prime(const mpz_t);
int miller_rabin_is_prime_k(const mpz_t, unsigned int);

#endif // PRIMES_MILLER_RABIN_H_

