/**
 * @author Ryan Lindeman
 * @class  CS6150
 * @description Implementation of the Atkin-Morain primality test using
 * Elliptic curve primality proving (ECPP)
 *
 * History
 * @date 20111029 - Initial release
 */

#include <time.h>
#include <stdio.h>
#include <stdint.h>
#include <gmp.h>
#include <assert.h>

// Global constants
const unsigned int MAX_DISCRIMINANTS = 28;
const unsigned int MAX_POINTS = 100;
const unsigned int MAX_PRIMES = 4294967295U;

// Global Variables
gmp_randstate_t gRandomState;  ///< Holds random generator state and algorithm type
mpz_t gD[MAX_DISCRIMINANTS];

// Point structure
struct Point
{
  mpz_t x;
  mpz_t y;
  mpz_t z;
};

/**
 * InitDiscriminants will initialize the gD array the appropriate
 * values.
 */
void InitDiscriminants(void)
{
  int32_t anD[MAX_DISCRIMINANTS] = {
  /* Place holder */ 0, 
  /* Class 1 */ -3,-4, -7,-8,-11,-19,-43,-67,-163,
  /* Class 2 */ -15,-20,-24,-35,-40,-51,-52,-88,-91,-115,-123,-148,-187,-232,
                -235,-267,-403,-427};

  // Call mpz_init on each array element and assign its value
  for(unsigned int i=0;i<MAX_DISCRIMINANTS;i++)
  {
    // Initialize and set each gD array element
    mpz_init_set_si(gD[i], anD[i]);

    // Display discriminants available
    //gmp_printf("%d=%Zd\n", i, gD[i]);
  }
}

/**
 * Computes s and t that satisfy 2^s * t = n, with t odd.
 */
void FactorPow2(mpz_t* s, mpz_t* t, mpz_t& n)
{
  mpz_set_ui(*s, 0);
  mpz_set(*t, n);
  while (mpz_even_p(*t))
  {
    mpz_add_ui(*s, *s, 1);
    mpz_divexact_ui(*t, *t, 2);
  }
}

/**
 * SquareMod returns the solution x to x^2 === a (mod p) which is used by the
 * ModifiedCornacchia to find the initial square root.
 */
bool SquareMod(mpz_t* theX, mpz_t& theA, mpz_t& theP)
{
  bool anResult = false; // Was a valid X found?
  
  // Perform the simple Jacobi test.  If Jacobi returns -1 or there is no
  // solution
  if(1 == mpz_jacobi(theA, theP))
  {
    mpz_t anMod4;
    mpz_t anMod8;

    // Initialize some temporary variables
    mpz_init(anMod4);
    mpz_init(anMod8);

    // Compute anMod4 first
    mpz_mod_ui(anMod4, theP, 4);

    // Compute anMod8 next
    mpz_mod_ui(anMod8, theP, 8);
    
    // Test the simplest cases p === 3, 7 (mod 8)
    if(mpz_cmp_ui(anMod4,3) == 0)
    {
      mpz_t anExp;  // The exponent

      // Initialize and set the exponent
      mpz_init_set(anExp, theP);

      // Add 1 to the exponent
      mpz_add_ui(anExp, anExp, 1);

      // Divide the exponent by 4
      mpz_div_ui(anExp, anExp, 4);

      // Compute x = a^(p+1)/4 mod p
      mpz_powm(*theX, theA, anExp, theP);
      
      // Clear the exponent
      mpz_clear(anExp);

      // We found a valid x
      anResult = true;
    }
    // Test the next simplest case p === 5 (mod 8)
    else if(mpz_cmp_ui(anMod8, 5) == 0)
    {
      mpz_t anExp1;  // The exponent (p-1)/4
      mpz_t anModP;  // Temporary value for second test

      // Initialize and set the exponent
      mpz_init_set(anExp1, theP);
      mpz_init(anModP);

      // Subtract from the exponent 1
      mpz_sub_ui(anExp1, anExp1, 1);

      // Divide the exponent by 4
      mpz_tdiv_q_ui(anExp1, anExp1, 4);

      // Compute x = a^(p-1)/4 mod p
      mpz_powm(anModP, theA, anExp1, theP);
      mpz_mod(anModP, anModP, theP);

      // See if the anModP result is 1
      if(mpz_cmp_ui(anModP,1) == 0)
      {
        mpz_t anExp2;  // The exponent (p+3)/8

        // Initialize and set the exponent
        mpz_init_set(anExp2, theP);

        // Add to the exponent 3
        mpz_add_ui(anExp2, anExp2, 3);

        // Divide the exponent by 8
        mpz_tdiv_q_ui(anExp2, anExp2, 8);

        // Compute x = a^(p+3)/8 mod p
        mpz_powm(*theX, theA, anExp2, theP);

        // Clear our values we don't need anymore
        mpz_clear(anExp2);
      }
      else
      {
        mpz_t anExp2;  // The exponent (p-5)/8
        mpz_t anBase;  // The base a*4

        // Initialize and set the exponent
        mpz_init_set(anExp2, theP);
        mpz_init(anBase);

        // Subtract from the exponent 5
        mpz_sub_ui(anExp2, anExp2, 5);

        // Divide the exponent by 8
        mpz_tdiv_q_ui(anExp2, anExp2, 8);

        // Compute 4*a as the base
        mpz_mul_ui(anBase, theA, 4);

        // Compute x = 4*a^(p-5)/8 mod p
        mpz_powm(anBase, anBase, anExp2, theP);

        // Multiply the result by 2*a*x mod p
        mpz_mul_ui(anBase, anBase, 2);
        mpz_mul(anBase, anBase, theA);
        mpz_mod(*theX, anBase, theP);

        // Clear our values we don't need anymore
        mpz_clear(anBase);
        mpz_clear(anExp2);
      }

      // Clear our values we don't need anymore
      mpz_clear(anModP);
      mpz_clear(anExp1);

      // We found a valid x
      anResult = true;
    }
    // Still nothing? Try the hardest case p === 1 (mod 8)
    else
    {
      // Case 2 of algorithm
      mpz_t d, s, t, A, D, m, i;
      mpz_t anExp1, anExp2, anModP, two;

      mpz_init(d);
      mpz_init(s);
      mpz_init(t);
      mpz_init(A);
      mpz_init(D);
      mpz_init(m);
      mpz_init(i);
      mpz_init(anExp1);
      mpz_init(anExp2);
      mpz_init(anModP);
      mpz_init(two);

      // Find a random integer d = [2, p - 1] with Jacobi -1
      mpz_sub_ui(anExp1, theP, 3);
      do {
        mpz_urandomm(d, gRandomState, anExp1);
        mpz_add_ui(d, d, 2);
      } while (mpz_jacobi(d, theP) != -1);
      
      // Represent p - 1 = 2 ^ s * t, with t odd
      mpz_sub_ui(anExp1, theP, 1);
      FactorPow2(&s, &t, anExp1);

      // Compute A = a ^ t mod p
      mpz_powm(A, theA, t, theP);

      // Compute D = d ^ t mod p
      mpz_powm(D, d, t, theP);

      // Compute -1 mod p
      mpz_sub_ui(anModP, theP, 1);
      
      // Compute m
      mpz_set_ui(two, 2);
      mpz_set_ui(m, 0);
      for (mpz_set_ui(i, 0); mpz_cmp(i, s) < 0; mpz_add_ui(i, i, 1))
      {
        // Compute 2 ^ (s - 1 - i)
        mpz_sub_ui(anExp1, s, 1);
        mpz_sub(anExp1, anExp1, i);
        mpz_powm(anExp1, two, anExp1, theP);

        // Compute A * D ^ m
        mpz_powm(anExp2, D, m, theP);
        mpz_mul(anExp2, anExp2, A);
        mpz_mod(anExp2, anExp2, theP);

        // Compute (A * D ^ m) ^ (2 ^ (s - 1 - i))
        mpz_powm(anExp1, anExp2, anExp1, theP);

        if (mpz_cmp(anExp1, anModP) == 0)
        {
          // Compute m = m + 2 ^ i
          mpz_powm(anExp1, two, i, theP);
          mpz_add(m, m, anExp1);
        }
      }

      // Compute a ^ ((t + 1) / 2)
      mpz_add_ui(anExp1, t, 1);
      mpz_divexact_ui(anExp1, anExp1, 2);
      mpz_powm(anExp1, theA, anExp1, theP);

      // Compute D ^ (m / 2)
      mpz_div_ui(anExp2, m, 2);
      mpz_powm(anExp2, D, anExp2, theP);
      
      // Compute x = a ^ ((t + 1) / 2) * D ^ (m / 2) mod p
      mpz_mul(anExp1, anExp1, anExp2);
      mpz_mod(anExp1, anExp1, theP);
      mpz_set(*theX, anExp1);
      
      mpz_clear(two);
      mpz_clear(anModP);
      mpz_clear(anExp2);
      mpz_clear(anExp1);
      mpz_clear(i);
      mpz_clear(m);
      mpz_clear(D);
      mpz_clear(A);
      mpz_clear(t);
      mpz_clear(s);
      mpz_clear(d);

      anResult = true;
    }

    // Clear our temporary variables before returning
    mpz_clear(anMod8);
    mpz_clear(anMod4);
  }

  // If no valid result was found, then set theX to 0
  if(false == anResult)
  {
    // Is this necessary?
    mpz_set_ui(*theX, 0);
  }

  // Return the result found if any
  return anResult;
}

/**
 * ModifiedCornacchia will either report that no solution exists for
 * 4p = u^2 + abs(D)v^2 (where p is a given prime and -4p < D < 0 or returns
 * the solution (u, v)
 */
bool ModifiedCornacchia(mpz_t* theU, mpz_t* theV, mpz_t& theP, mpz_t& theD)
{
  bool anResult = false; // No solution found yet
  mpz_t x0; // Initial square root value for x
  mpz_t a;  // Euclid chain value a
  mpz_t b;  // Euclid chain value b
  mpz_t c;  // Euclid chain value c
  mpz_t t;  // Temporary value for final report

  // Initalize square root value
  mpz_init(x0);
  mpz_init(a);
  mpz_init(b);
  mpz_init(c);
  mpz_init(t);

  // Obtain the initial square root value for x0
  bool anFound = SquareMod(&x0, theD, theP); // uses 2.3.8 algorithm

  // If x0 was found, then ensure x0^2 !=== D (mod 2)
  if(anFound)
  {
    mpz_t xt; // Temporary value x0 mod 2
    mpz_t dt; // Temporary value theD mod 2

    // Initalize some temporary values
    mpz_init(xt);
    mpz_init(dt);

    // Compute x0 mod 2 and D mod 2
    mpz_mod_ui(dt, theD, 2);
    mpz_mod_ui(xt, x0, 2);

    // If x0^2 !=== D (mod 2) then adjust x0
    if(mpz_cmp(dt,xt) != 0)
    {
      mpz_sub(x0, theP, x0);
    }

    // Clear our temporary values
    mpz_clear(dt);
    mpz_clear(xt);
  } else {
    // Warn the user and proceed anyway with x0 = 0
    printf("x0 not found!?\n");
  }

  // Initialize Euclid chain
  // Set a = 2p
  mpz_set(a, theP);
  mpz_mul_ui(a, a, 2);

  // Set b = x0
  mpz_set(b, x0);

  // Set c = lower_bound(2*sqrt(p)) (see 9.2.11 algorithm)
  mpz_sqrt(c, theP);
  mpz_mul_ui(c, c, 2);

  // Euclid chain
  while(mpz_cmp(b, c) > 0)
  {
    mpz_t anModB; // Temporary result for a mod b

    // Initialize our temporary value for a mod b
    mpz_init(anModB);

    // Compute a mod b first
    mpz_mod(anModB, a, b);

    // Set a = b
    mpz_set(a, b);

    // Set b = anModB
    mpz_set(b, anModB);

    // Clear our temporary result for a mod b now
    mpz_clear(anModB);
  } // while(b > c)

  // Final report/check
  // Compute a = b^2 to use later
  mpz_pow_ui(a, b, 2);

  // Compute c = 4p to use later
  mpz_mul_ui(c, theP, 4);

  // t = 4p - y^2;
  mpz_sub(t, c, a);

  // Compute a = abs(D)
  mpz_abs(a, theD);

  // Compute c = t / a
  mpz_tdiv_q(c, t, a);

  // Compute a = t mod a
  mpz_mod(a, t, a);

  // Now test if solution was found
  if(mpz_cmp_ui(a, 0) == 0 && mpz_perfect_square_p(c) != 0)
  {
    // Set theU value for our solution
    mpz_set(*theU, b);

    // Set theV value for our solution
    mpz_sqrt(*theV, c);

    // We found a solution
    anResult = true;
  }

  // Clear our values used above
  mpz_clear(t);
  mpz_clear(c);
  mpz_clear(b);
  mpz_clear(a);
  mpz_clear(x0);

  // Return anResult which is true if solution was found
  return anResult;
}

/**
 * Preforms ellpitic addtion of two affine coordinates. The argument d 
 * will also be passed back with the divisor used in the addtion.
 */
void EllipticAdd(Point& p, Point& p1, Point& p2, mpz_t d, mpz_t a) {  
  if (mpz_cmp_ui(p1.z, 0) == 0) // z1 == 0
  {
    mpz_set(p.x, p2.x);
    mpz_set(p.y, p2.y);
    mpz_set(p.z, p2.z);
    return;
  }
  else if (mpz_cmp_ui(p2.z, 0) == 0) // z2 == 0
  {
    mpz_set(p.x, p1.x);
    mpz_set(p.y, p1.y);
    mpz_set(p.z, p1.z);
    return;
  }
  
  mpz_t m;
  mpz_init(m);
  
  if (mpz_cmp(p1.x, p2.x) == 0) // x1 == x2
  {
    mpz_add(m, p1.y, p2.y);
    if (mpz_cmp_ui(m, 0) == 0)
    {
      mpz_set_ui(p.x, 0);
      mpz_set_ui(p.y, 1);
      mpz_set_ui(p.z, 0);
      mpz_clear(m);
      return;
    }
    // m = (3 * x1 ^ 2 + a) / (2 * y1 )
    mpz_mul_ui(d, p1.y, 2); 
    mpz_pow_ui(m, p1.x, 2);
    mpz_mul_ui(m, m, 3);
    mpz_add(m, m, a);
    mpz_div(m, m, d);
  }
  else
  {
    // m = (y2 - y1) / (x2 - x1)
    mpz_sub(d, p2.x, p1.x); 
    mpz_sub(m, p2.y, p1.y);
    mpz_div(m, m, d);
  }
  
  mpz_t x3, y3;
  mpz_init(x3);
  mpz_init(y3);
  
  // x3 = m ^ 2 - x1 - x2
  mpz_pow_ui(x3, m, 2);
  mpz_sub(x3, x3, p1.x);
  mpz_sub(x3, x3, p2.x);
  
  // m * ( x1 - x3) - y1
  mpz_sub(y3, p1.x, x3);
  mpz_mul(y3, y3, m);
  mpz_sub(y3, y3, p1.y);
  
  mpz_set(p.x, x3);
  mpz_set(p.y, y3);
  mpz_set_ui(p.z, 1);
  
  mpz_clear(m);
  mpz_clear(x3);
  mpz_clear(y3);
}

/**
 * Sets n to the next prime that does exceed the limit. If there is
 * not one, n is set to the limit. This is very inefficient so 
 * should only be used with small numbers.
 * 
 * We could possibly add a precomputed table of the first <limit> primes.
 */
void NextPrime(mpz_t n, unsigned long limit)
{
  mpz_t factor, sqrtN, remainder;
  
  mpz_init(factor);
  mpz_init(sqrtN);
  mpz_init(remainder);
  
  mpz_add_ui(n, n, 1);
  while (mpz_cmp_ui(n, limit) < 0)
  {
    bool isPrime = true;
    mpz_sqrt(sqrtN, n);
    for (mpz_set_ui(factor, 2); mpz_cmp(factor, sqrtN) <= 0; mpz_add_ui(factor, factor, 1))
    {
      mpz_mod(remainder, n, factor);
      if (mpz_cmp_ui(remainder, 0) == 0)
      {
        isPrime = false;
        break;
      }
    }
    if (isPrime)
    {
      break;
    }
    mpz_add_ui(n, n, 1);
  }
  
  mpz_clear(factor);
  mpz_clear(sqrtN);
  mpz_clear(remainder);
}

/**
 * LenstraECM will attempt to find the largest non-trivial prime number that
 * will factor theN provided. This is called by FindFactor to find a probable
 * prime theQ that if proven to be prime will prove theN to be prime.
 */
void LenstraECM(mpz_t* theQ, mpz_t& theN)
{
  unsigned long B1 = 10000;
  Point P;  // Point (x,y) on curve E
  mpz_t a;  // Random value a for curve E
  mpz_t b;  // Random value b for curve E
  mpz_t g;  // The factor found if any
  mpz_t t;  // Temporary value for computing b
  mpz_t p, d;

  // Initialize our values
  mpz_init(P.x);
  mpz_init(P.y);
  mpz_init(a);
  mpz_init(b);
  mpz_init(g);
  mpz_init(t);
  mpz_init(p);
  mpz_init(d);

  do {
    // Pick a random x from 0 to N-1
    mpz_urandomm(P.x, gRandomState, theN);
    // Pick a random y from 0 to N-1
    mpz_urandomm(P.y, gRandomState, theN);
    // Pick a random a from 0 to N-1
    mpz_urandomm(a, gRandomState, theN);

    // Compute b = (y^2 - x^3 - ax) mod n
    mpz_pow_ui(b, P.y, 2); // y^2
    mpz_pow_ui(t, P.x, 3); // x^3
    mpz_sub(b, b, t); // y^2 - x^3
    mpz_mul(t, a, P.x); // ax
    mpz_sub(b, b, t); // y^2 - x^3 - ax
    mpz_mod(b, b, theN);  // (y^2 - x^3 - ax) mod n

    // Compute 4a^3 + 27b^2
    mpz_pow_ui(a, a, 3);  // a^3
    mpz_mul_ui(a, a, 4);  // 4*(a^3)
    mpz_pow_ui(b, b, 2);  // b^2
    mpz_mul_ui(b, b, 27); // 27*(b^2)
    mpz_add(t, a, b);     // (4a^3) + (27b^2)
    mpz_gcd(g, t, theN);  // Compute g = gcd(4a^3+27b^2, n)
    // Continue until we find a g != n
    if(mpz_cmp(g, theN) == 0)
      continue;
    // If G is greater than 1 but not equal to n then we're done
    else if(mpz_cmp_ui(g, 1) > 0)
    {
      gmp_printf("g=%Zd\n", g);
      break;
    }
    
    bool factorFound = false;
    for (mpz_set_ui(p, 1); mpz_cmp_ui(p, B1) < 0; NextPrime(p, B1))
    {
      mpz_pow_ui(t, p, 1);
      for (unsigned long j = 1; mpz_cmp_ui(t, B1) <= 0; mpz_pow_ui(t, p, ++j))
      {
        for (int k = 0; mpz_cmp_ui(p, k) > 0; k++)
        {
          EllipticAdd(P, P, P, d, a);
          mpz_gcd(g, theN, d);
          if (mpz_cmp_ui(g, 1) != 0 && mpz_cmp(g, theN) != 0) {
            mpz_set(theN, g);
            factorFound = true;
            break;
          }
        }
        if (factorFound) break;
      }
      if (factorFound) break;
    }
    if (factorFound) break;
    
    // Possibly increment B1
    
  } while(true);

  // Clear our values used above
  mpz_clear(t);
  mpz_clear(g);
  mpz_clear(b);
  mpz_clear(a);
  mpz_clear(P.y);
  mpz_clear(P.x);
  mpz_clear(p);
  mpz_clear(d);
}

/**
 * FindFactor will attempt to reduce theM provided to a potential prime theQ
 * that is also less than theT by factoring all easy primes out of theM.  If
 * theM is a composite then FindFactor will return false meaning theQ wasn't
 * found because m is a composite.
 */
bool FindFactor(mpz_t* theQ, mpz_t& theM, mpz_t& theT)
{
  bool anResult = true; // A suitable theQ was found
  unsigned long count = MAX_PRIMES; // Primes to try to remove from theQ
  mpz_t prime;  // Prime numbers to remove from theQ

  // If theM is prime or probably prime then stop now
  if(mpz_probab_prime_p(theM, 10))
    return false;

  // Starting with theM value provided, try to reduce it to
  // some form of m = k*q by finding k factors that can be
  // removed from m leaving q as a possible prime factor
  mpz_set(*theQ, theM);

  // TODO: Replace the following with LenstraECM implementation
  // Initialize our first prime number to 2
  mpz_init_set_ui(prime, 2);

  // Loop through each prime number from 2 up and remove each
  // of these from theQ
  do 
  {
    // Remove as many of this prime as possible
    mpz_remove(*theQ, *theQ, prime);

    // Find the next prime to try to remove from theQ
    mpz_nextprime(prime, prime);

    // Is theQ prime now? then exit our loop
    if(mpz_probab_prime_p(*theQ, 10))
      break;
  } while(mpz_cmp(*theQ, theT) >= 0 && count-- > 0);

  // Let the user know that our FindFactors algorithm might have failed!
  if(count == 0)
  {
    printf("Warning q might not be prime!\n");
  }
  
  // Clear our prime value used above
  mpz_clear(prime);

  // Make sure theQ is still larger than theT
  if(mpz_cmp(*theQ, theT) < 0)
    return false; // theQ is smaller than theT

  // One more final check, m != q
  if(mpz_cmp(*theQ, theM) == 0)
    return false;

  return anResult;
}

/**
 * FactorOrders attempts to find a possible order m that factors as m = kq
 * where k > 1 and q is a probable prime > (n^0.25 + 1)^2. If this can't be
 * done after K_max iterations than return FALSE and choose a new discriminant
 * D and curve m.
 */
bool FactorOrders(mpz_t* theM, mpz_t* theQ, mpz_t& theU, mpz_t& theV,
                  mpz_t& theN, mpz_t& theD)
{
  bool anResult = false;  // Was factor theQ found?
  mpz_t t;  // t = (n^0.25 + 1)^2 to test theQ with
  mpz_t m0;  // m0 = n + 1
  mpz_t m1;  // m1 = m0 + u
  mpz_t m2;  // m2 = m0 - u

  // Initialize and set our temporary value t to theN
  mpz_init(t);
  mpz_init_set(m0, theN);
  mpz_init(m1);
  mpz_init(m2);

  // Add 1 to m0
  mpz_add_ui(m0, m0, 1);

  // Now take the double square root of theN and store in T
  mpz_sqrt(t, theN);
  mpz_sqrt(t, t);

  // Now add one to the result
  mpz_add_ui(t, t, 1);

  // Now square the result
  mpz_pow_ui(t, t, 2);

  // Test initial special cases of m0 = n + 1 +/- u
  mpz_add(m1, m0, theU);  // Add u to m0 for m1
  mpz_sub(m2, m0, theU);  // Subtract u from m0 for m2

  if(true == FindFactor(theQ, m1, t))
  {
    // Set m = m1
    mpz_set(*theM, m1);

    // Curve order m and factor q found
    anResult = true;
  }
  else if(true == FindFactor(theQ, m2, t))
  {
    // Set m = m2
    mpz_set(*theM, m2);

    // Curve order m and factor q found
    anResult = true;
  }
  else if(mpz_cmp_si(theD, -4) == 0)
  {
    // Test next special cases m1 = m0 + 2v and m2 = m0 - 2v
    mpz_set(m1, m0);
    mpz_set(m2, m0);
    mpz_addmul_ui(m1, theV, 2);
    mpz_submul_ui(m2, theV, 2);
    if(true == FindFactor(theQ, m1, t))
    {
      // Set m = m1
      mpz_set(*theM, m1);

      // Curve order m and factor q found
      anResult = true;
    }
    else if(true == FindFactor(theQ, m2, t))
    {
      // Set m = m2
      mpz_set(*theM, m2);

      // Curve order m and factor q found
      anResult = true;
    }
    else
    {
      // No curve found for D = -4
    }
  }
  else if(mpz_cmp_si(theD, -3) == 0)
  {
    mpz_t m3; // (u + 3v) / 2
    mpz_t m4; // (u - 3v) / 2
    mpz_init_set(m3, theU);
    mpz_init_set(m4, theU);

    // Prepare factor m3 = (u + 3v) / 2
    mpz_addmul_ui(m3, theV, 3);
    mpz_tdiv_q_ui(m3, m3, 2);

    // Prepare factor m4 = (u - 3v) / 2
    mpz_submul_ui(m4, theV, 3);
    mpz_tdiv_q_ui(m4, m4, 2);

    // Test special case m1 = m0 + m3 and m2 = m0 + m4
    mpz_add(m1, m0, m3);
    mpz_add(m2, m0, m4);
    if(true == FindFactor(theQ, m1, t))
    {
      // Set m = m1
      mpz_set(*theM, m1);

      // Curve order m and factor q found
      anResult = true;
    }
    else if(true == FindFactor(theQ, m2, t))
    {
      // Set m = m2
      mpz_set(*theM, m2);

      // Curve order m and factor q found
      anResult = true;
    }
    else
    {
      // Test special case m1 = m0 - m3 and m2 = m0 - m4
      mpz_sub(m1, m0, m3);
      mpz_sub(m2, m0, m4);
      if(true == FindFactor(theQ, m1, t))
      {
        // Set m = m1
        mpz_set(*theM, m1);

        // Curve order m and factor q found
        anResult = true;
      }
      else if(true == FindFactor(theQ, m2, t))
      {
        // Set m = m2
        mpz_set(*theM, m2);

        // Curve order m and factor q found
        anResult = true;
      }
      else
      {
        // No curve found for D = -3
      }
    }

    // Clear our values used above
    mpz_clear(m4);
    mpz_clear(m3);
  }
  else
  {
    // No curve found for D < -4
  }

  // Clear our values used above
  mpz_clear(m2);
  mpz_clear(m1);
  mpz_clear(m0);
  mpz_clear(t);

  // Return anResult found if any
  return anResult;
}

/**
 * CalculateNonresidue will find a random quadratic nonresidue g mod p and if
 * theD == -3 then g must also be a cubic nonresidue. This is based on step 1
 * of algorithm (7.5.9) or step 2b of algorithm (7.5.10)
 */
void CalculateNonresidue(mpz_t* theG, mpz_t& theN, mpz_t& theD)
{
  mpz_t t;  // t = theG^(n_div_3) mod n
  mpz_t n_div_3;  // (N - 1) / 3

  // Initialize our temporary value of (N - 1) / 3
  mpz_init(t);
  mpz_init_set(n_div_3, theN);

  // Compute (N - 1) / 3
  mpz_sub_ui(n_div_3, n_div_3, 1);
  mpz_tdiv_q_ui(n_div_3, n_div_3, 3);

  do
  {
    // Pick a random x from 0 to N-1
    mpz_urandomm(*theG, gRandomState, theN);
    //mpz_set_ui(*theG, 1000021176);

    // Eliminate 0 as a possible theG value
    if(mpz_cmp_ui(*theG, 0) == 0)
      continue;

    // Make sure it passes the Jacobi test
    if(-1 != mpz_jacobi(*theG, theN))
      continue; // Jacobi returned -1, try another g

    // Make sure g is a cubic nonresidue
    if(mpz_cmp_si(theD, -3) == 0)
    {
      mpz_powm(t, *theG, n_div_3, theN);

      // If our result is equal to 1, then its not cubic nonresidue
      if(mpz_cmp_ui(t, 1) == 0)
        continue;
    }

    // If we got to here we can quit, theG is good to use
    break;
  } while (true);

  // Clear our temporary value
  mpz_clear(n_div_3);
  mpz_clear(t);
}

// Lookup Table 7.1
bool LookupCurveParameters(mpz_t* r, mpz_t* s, mpz_t& p, mpz_t& d)
{
  signed long int D = mpz_get_si(d);
  bool isValid = true;
  mpz_t tmp1, tmp2;

  // Initialize our temporary variables
  mpz_init(tmp1);
  mpz_init(tmp2);

  // Find the discriminant in our table of class 1 and 2 discriminants
  switch(D) {
    case -7:
      mpz_set_ui(*r, 125);
      mpz_set_ui(*s, 189);
      break;
    case -8:
      mpz_set_ui(*r, 125);
      mpz_set_ui(*s, 98);
      break;
    case -11:
      mpz_set_ui(*r, 512);
      mpz_set_ui(*s, 539);
      break;
    case -19:
      mpz_set_ui(*r, 512);
      mpz_set_ui(*s, 513);
      break;
    case -43:
      mpz_set_ui(*r, 512000);
      mpz_set_ui(*s, 512001);
      break;
    case -67:
      mpz_set_ui(*r, 85184000);
      mpz_set_ui(*s, 85184001);
      break;
    case -163:
      mpz_set_str(*r, "151931373056000", 10);
      mpz_set_str(*s, "151931373056001", 10);
      break;
    case -15:
      // 1225 - 2080 * sqrt(5)
      mpz_set_ui(tmp1, 5);
      isValid = SquareMod(&tmp1, tmp1, p);
      mpz_set_ui(*r, 1225);
      mpz_set_ui(*s, 5929);
      mpz_submul_ui(*r, tmp1, 2080);
      break;
    case -20:
      // 108250 + 29835 * sqrt(5)
      mpz_set_ui(tmp1, 5);
      isValid = SquareMod(&tmp1, tmp1, p);
      mpz_set_ui(*r, 108250);
      mpz_set_ui(*s, 174724);
      mpz_addmul_ui(*r, tmp1, 29835);
      break;
    case -24:
      // 1757 - 494 * sqrt(2)
      mpz_set_ui(tmp1, 2);
      isValid = SquareMod(&tmp1, tmp1, p);
      mpz_set_ui(*r, 1757);
      mpz_set_ui(*s, 1058);
      mpz_submul_ui(*r, tmp1, 494);
      break;
    case -35:
      // -1126400 - 1589760 * sqrt(5)
      mpz_set_ui(tmp1, 5);
      isValid = SquareMod(&tmp1, tmp1, p);
      mpz_set_si(*r, -1126400);
      mpz_set_ui(*s, 2428447);
      mpz_submul_ui(*r, tmp1, 1589760);
      break;
    case -40:
      // 54175 - 1020 * sqrt(5)
      mpz_set_ui(tmp1, 5);
      isValid = SquareMod(&tmp1, tmp1, p);
      mpz_set_ui(*r, 54175);
      mpz_set_ui(*s, 51894);
      mpz_submul_ui(*r, tmp1, 1020);
      break;
    case -51:
      // 75520 - 7936 * sqrt(5)
      mpz_set_ui(tmp1, 17);
      isValid = SquareMod(&tmp1, tmp1, p);
      mpz_set_ui(*r, 75520);
      mpz_set_ui(*s, 108241);
      mpz_submul_ui(*r, tmp1, 7936);
      break;
    case -52:
      // 1778750 + 5125 * sqrt(13)
      mpz_set_ui(tmp1, 13);
      isValid = SquareMod(&tmp1, tmp1, p);
      mpz_set_ui(*r, 1778750);
      mpz_set_ui(*s, 1797228);
      mpz_addmul_ui(*r, tmp1, 5125);
      break;
    case -88:
      // 181713125 - 44250  * sqrt(2)
      mpz_set_ui(tmp1, 2);
      isValid = SquareMod(&tmp1, tmp1, p);
      mpz_set_ui(*r, 181713125);
      mpz_set_ui(*s, 181650546);
      mpz_submul_ui(*r, tmp1, 44250);
      break;
    case -91:
      // 74752 - 36352 * sqrt(13)
      mpz_set_ui(tmp1, 13);
      isValid = SquareMod(&tmp1, tmp1, p);
      mpz_set_ui(*r, 74752);
      mpz_set_ui(*s, 205821);
      mpz_submul_ui(*r, tmp1, 36352);
      break;
    case -115:
      // 269593600 - 89157120  * sqrt(5)
      mpz_set_ui(tmp1, 5);
      isValid = SquareMod(&tmp1, tmp1, p);
      mpz_set_ui(*r, 269593600);
      mpz_set_ui(*s, 468954981);
      mpz_submul_ui(*r, tmp1, 89157120);
      break;
    case -123:
      // 1025058304000 - 1248832000 * sqrt(41)
      mpz_set_ui(tmp1, 41);
      isValid = SquareMod(&tmp1, tmp1, p);
      mpz_set_str(*r, "1025058304000", 10);
      mpz_set_str(*s, "1033054730449", 10);
      mpz_submul_ui(*r, tmp1, 1248832000);
      break;
    case -148:
      // 499833128054750 + 356500625 * sqrt(37)
      mpz_set_ui(tmp1, 37);
      isValid = SquareMod(&tmp1, tmp1, p);
      mpz_set_str(*r, "499833128054750", 10);
      mpz_set_str(*s, "499835296563372", 10);
      mpz_addmul_ui(*r, tmp1, 356500625);
      break;
    case -187:
      // 91878880000 - 1074017568000 * sqrt(17)
      mpz_set_ui(tmp1, 17);
      mpz_set_str(tmp2, "1074017568000", 10);
      isValid = SquareMod(&tmp1, tmp1, p);
      mpz_set_str(*r, "91878880000", 10);
      mpz_set_str(*s, "4520166756633", 10);
      mpz_submul(*r, tmp1, tmp2);
      break;
    case -232:
      // 1728371226151263375 - 11276414500 * sqrt(29)
      mpz_set_ui(tmp1, 29);
      mpz_set_str(tmp2, "11276414500", 10);
      isValid = SquareMod(&tmp1, tmp1, p);
      mpz_set_str(*r, "1728371226151263375", 10);
      mpz_set_str(*s, "1728371165425912854", 10);
      mpz_submul(*r, tmp1, tmp2);
      break;
    case -235:
      // 7574816832000 - 190341944320 * sqrt(5)
      mpz_set_ui(tmp1, 5);
      mpz_set_str(tmp2, "190341944320", 10);
      isValid = SquareMod(&tmp1, tmp1, p);
      mpz_set_str(*r, "7574816832000", 10);
      mpz_set_str(*s, "8000434358469", 10);
      mpz_submul(*r, tmp1, tmp2);
      break;
    case -267:
      // 3632253349307716000000 - 12320504793376000 * sqrt(89)
      mpz_set_ui(tmp1, 89);
      mpz_set_str(tmp2, "12320504793376000", 10);
      isValid = SquareMod(&tmp1, tmp1, p);
      mpz_set_str(*r, "3632253349307716000000", 10);
      mpz_set_str(*s, "3632369580717474122449", 10);
      mpz_submul(*r, tmp1, tmp2);
      break;
    case -403:
      // 16416107434811840000 - 4799513373120384000 * sqrt(13)
      mpz_set_ui(tmp1, 13);
      mpz_set_str(tmp2, "4799513373120384000", 10);
      isValid = SquareMod(&tmp1, tmp1, p);
      mpz_set_str(*r, "16416107434811840000", 10);
      mpz_set_str(*s, "33720998998872514077", 10);
      mpz_submul(*r, tmp1, tmp2);
      break;
    case -427:
      // 564510997315289728000 - 5784785611102784000 * sqrt(61)
      mpz_set_ui(tmp1, 61);
      mpz_set_str(tmp2, "5784785611102784000", 10);
      isValid = SquareMod(&tmp1, tmp1, p);
      mpz_set_str(*r, "564510997315289728000", 10);
      mpz_set_str(*s, "609691617259594724421", 10);
      mpz_submul(*r, tmp1, tmp2);
      break;
    default:
      isValid = false;
      break;
  }

  // Did we set an r and s value above?
  if(isValid)
  {  
    // keep things mod p
    mpz_mod(*r, *r, p);
    mpz_mod(*s, *s, p);
  }

  // Clear our temporary variables used above
  mpz_clear(tmp2);
  mpz_clear(tmp1);

  // Return true if *r and *s were set above
  return isValid;
}

/**
 * ObtainCurveParameters will attempt to obtain the curve parameters a and b
 * for an elliptic curve that would have order m if n is indeed prime. This is
 * based on steps 5-7 of algorithm (7.5.9) or steps 3-5 of algorithm (7.5.10).
 */
bool ObtainCurveParameters(mpz_t* theA, mpz_t* theB, mpz_t& theN,
                           mpz_t& theD, mpz_t& theG, unsigned int theK)
{
  // theA, theB: curve parameters a and b computed by this method
  // theG: 'suitable nonresidue of p'
  // theK: k-th iteration of ObtainCurveParameters
  // theN: number we are currently testing for primality
  // theD: current discriminant from global discriminant array gD

  bool anResult = true; // True if curve parameters a and b were found

  // Special case for D=-3 and D=-4
  if(mpz_cmp_si(theD, -3) == 0)
  {
    // If theK >= 6 then no curve parameters will be returned
    if(theK >= 6)
    {
      // No curve parameters available, try another discriminant
      anResult = false;
    }
    else
    {
      // a = 0 in all k cases below
      mpz_set_ui(*theA,0);

      // Quickly handle special case of k=0: b=-1 mod n
      if(theK == 0)
      {
        // Set b = -1 mod n = n - 1
        mpz_set(*theB,theN);
        mpz_sub_ui(*theB, *theB, 1);
        //mpz_set_str(*theB, "273047738243666261637604519", 10);
      }
      // Handle other cases for theK < 6
      else
      {
        // Set b = b * g   (b = -g^k)
        mpz_mul(*theB, *theB, theG);
        mpz_mod(*theB, *theB, theN);

        // Begin Sanity Check
        mpz_t tmpG;
        mpz_init(tmpG);

        // tmpG = -g^k mod n
        mpz_pow_ui(tmpG, theG, theK);
        mpz_mul_si(tmpG, tmpG, -1);
        mpz_mod(tmpG, tmpG, theN);

        assert(mpz_cmp(*theB, tmpG) == 0);
        mpz_clear(tmpG);
        // End Sanity Check
      }
    }
  }
  else if(mpz_cmp_si(theD, -4) == 0)
  {
    // If theK >= 4 then no curve parameters will be returned
    if(theK >= 4)
    {
      // No curve parameters available, try another discriminant
      anResult = false;
    }
    else
    {
      // b = 0 in all k cases below
      mpz_set_ui(*theB,0);

      // Quickly handle special case of k=0: a=-1 mod n
      if(theK == 0)
      {
        // a=-1 mod n, b=0
        mpz_set(*theA,theN);
        mpz_sub_ui(*theA, *theA, 1);
      }
      // Handle other cases for theK < 4
      else
      {
        // Set a = a * g
        mpz_mul(*theA, *theA, theG);
        mpz_mod(*theA, *theA, theN);
      }
    }
  }
  else
  {
    // If theK >= 2 then no curve parameters will be returned
    if(theK >= 2)
    {
      // No curve parameters available, try another discriminant
      anResult = false;
    }
    // Case 1: k==0, call LookupCurveParameters to obtain the r and s values
    else if(theK == 0)
    {
      mpz_t r; // r value returned by LookupCurveParameters
      mpz_t s; // s value returned by LookupCurveParameters
      
      // Initialize our r and s values before calling LookupCurveParameters
      mpz_init(r);
      mpz_init(s);

      // Determine if theD can be found by LookupCurveParameters
      anResult = LookupCurveParameters(&r, &s, theN, theD);
      gmp_printf("r=%Zd, s=%Zd\n", r, s);

      // Did LookupCurveParameters yield an r and s value?
      if(anResult)
      {
        // Compute a = -3rs^3g^(2k) with k==0
        mpz_powm_ui(*theA, s, 3, theN);
        mpz_mul(*theA, *theA, r);
        mpz_mul_si(*theA, *theA, -3);
        mpz_mod(*theA, *theA, theN);

        // Compute b = 2rs^5g^(3k) with k==0
        mpz_powm_ui(*theB, s, 5, theN);
        mpz_mul(*theB, *theB, r);
        mpz_mul_ui(*theB, *theB, 2);
        mpz_mod(*theB, *theB, theN);
      }

      // Clear our r and s values provided by LookupCurveParameters
      mpz_clear(s);
      mpz_clear(r);
    }
    // Case 2: k==1, just multiply previous a and b values with g provided
    else
    {
      mpz_t tmpG; // g^2 value used in both multiplications below
      
      // Initialize our temporary g^2 value
      mpz_init(tmpG);

      // Compute g^2 which is used for both a and b below
      mpz_pow_ui(tmpG, theG, 2);
      
      // Compute a = a*g^(2k) with k==1
      mpz_mul(*theA, *theA, tmpG);
      mpz_mod(*theA, *theA, theN);

      // Compute b = b*g^(3k) with k==1
      mpz_mul(*theB, *theB, tmpG);
      mpz_mul(*theB, *theB, theG);
      mpz_mod(*theB, *theB, theN);

      // Clear our temporary g^2 value
      mpz_clear(tmpG);
    }
  }

  // Return true if we found our curve parameters a and b above
  return anResult;
}

/**
 * ChoosePoint will try and find a point (x,y) on the curve using the a and b
 * values provided according to algorithm (7.2.1). This will return true if
 * during the process of finding P, ChoosePoint discovers that N is composite,
 * false otherwise.
 */
bool ChoosePoint(struct Point* theP, mpz_t& theN, mpz_t& theA, mpz_t& theB)
{
  bool anResult = false;  // True if we discover theN is composite
  bool anFound = false;  // True when we find a suitable P
  mpz_t Q;
  mpz_t t;

  // Initialize our temporary values
  mpz_init(Q);
  mpz_init(t);

  do
  {
    // Step 4a: Choose a random x such that Q = (x^3 + ax + b) mod n and
    // Jacobi(Q/N) != -1.
    do
    {
      // Pick a random x from 0 to N-1
      mpz_urandomm(theP->x, gRandomState, theN);
      //mpz_set_ui(theP->x, 946781885);

      // Compute Q = x^3
      mpz_powm_ui(Q, theP->x, 3, theN);

      // Compute Q = Q + ax
      mpz_addmul(Q, theA, theP->x);

      // Compute Q = Q + b
      mpz_add(Q, Q, theB);

      // Compute Q = Q mod n
      mpz_mod(Q, Q, theN);
    } while(-1 == mpz_jacobi(Q,theN));

    // Step 4b: Apply Algorithm 2.3.8 or 2.3.9 (with a = Q and p = n) to find
    // an integer y that would satisfy y2 = Q (mod n) if n were prime
    anFound = SquareMod(&theP->y, Q, theN);
  } while(false == anFound);

  // Step 4c: If y^2 mod n != Q then N is composite
  mpz_powm_ui(t, theP->y, 2, theN);
  if(mpz_cmp(Q, t) != 0)
  {
    // N is composite
    anResult = true;
  }

  // Clear the temporary values used above
  mpz_clear(t);
  mpz_clear(Q);

  // Return true if we found a composite, false otherwise
  return anResult;
}

/**
 * Add implements the Elliptic add method described by algorithm (7.2.2). This
 * will return true if an illegal inversion occurred and will not set theR
 * unless the inversion succeeded.
 */
bool Add(struct Point* theR, struct Point& theP1, struct Point& theP2,
         mpz_t& theN)
{
  bool anResult = false;  // True if illegal inversion occurred
  mpz_t m;  // Elliptic slope value m

  // Initialize our m value first
  mpz_init(m);

  // Does P1 and P2 possibly represent inverse points?
  if(mpz_cmp(theP1.x, theP2.x) == 0)
  {
    // Compute (y2 + y1) mod n
    mpz_add(m, theP1.y, theP2.y);
    mpz_mod(m, m, theN);

    // Are we at the infinite point O? then return O
    if(mpz_cmp_ui(m, 0) == 0)
    {
      // Set our return result as O
      mpz_set_ui(theR->x, 0);
      mpz_set_ui(theR->y, 1);

      // Clear our m value used above
      mpz_clear(m);

      // Illegal inversion at the infinite point O
      return true;
    }
  }

  // Compute (x2 - x1) portion of m = (y2-y1)(x2-x1)^-1
  mpz_sub(m, theP2.x, theP1.x);

  // Compute inverse (x2-x1)^-1 first
  anResult = (mpz_invert(m, m, theN) == 0);

  // Only proceed if our inversion succeeded
  if(false == anResult)
  {
    Point R;  // Point R to return via theR

    // Initialize our working point R
    mpz_init(R.x);
    mpz_init(R.y);

    // Compute (y2 - y1)
    mpz_sub(R.x, theP2.y, theP1.y);

    // Compute m = ((y2-y1)(x2-x1)^-1) mod n
    mpz_mul(m, R.x, m);
    mpz_mod(m, m, theN);

    // Compute R.x = m^2 - x1 - x2
    mpz_powm_ui(R.x, m, 2, theN);
    mpz_sub(R.x, R.x, theP1.x);
    mpz_sub(R.x, R.x, theP2.x);
    mpz_mod(R.x, R.x, theN);

    // Compute R.y = m(x1 - x3) - y1
    mpz_sub(R.y, theP1.x, R.x);
    mpz_mul(R.y, R.y, m);
    mpz_sub(R.y, R.y, theP1.y);
    mpz_mod(R.y, R.y, theN);

    // Set our return results now
    mpz_set(theR->x, R.x);
    mpz_set(theR->y, R.y);

    // Clear our working point R
    mpz_clear(R.y);
    mpz_clear(R.x);
  }

  // Clear our m value used above
  mpz_clear(m);

  // Return true if illegal inversion occurred
  return anResult;
}

/**
 * Double implements the double portion of the Elliptic add method described by
 * algorithm (7.2.2).
 */
void Double(struct Point* theR, struct Point& theP, mpz_t& theN, mpz_t& theA)
{
  mpz_t m;  // Elliptic slope value m
  Point R;  // Point R to return via theR

  // Initialize our values first
  mpz_init(m);
  mpz_init(R.x);
  mpz_init(R.y);

  // Compute 2*P.y first
  mpz_mul_ui(m, theP.y, 2);

  // Compute inverse (2y)^-1 first
  mpz_invert(m, m, theN);

  // Compute (3x^2 + a)
  mpz_mul(R.x, theP.x, theP.x);
  mpz_mul_ui(R.x, R.x, 3);
  mpz_add(R.x, R.x, theA);

  // Compute m = (3x^2 + a)(2y)^-1
  mpz_mul(m, R.x, m);

  // Compute R.x = m^2 - 2x
  mpz_powm_ui(R.x, m, 2, theN);
  mpz_submul_ui(R.x, theP.x, 2);
  mpz_mod(R.x, R.x, theN);

  // Compute R.y = m(x1 - x3) - y1
  mpz_sub(R.y, theP.x, R.x);
  mpz_mul(R.y, R.y, m);
  mpz_sub(R.y, R.y, theP.y);
  mpz_mod(R.y, R.y, theN);

  // Set our return results now
  mpz_set(theR->x, R.x);
  mpz_set(theR->y, R.y);

  // Clear our values used above
  mpz_clear(R.y);
  mpz_clear(R.x);
  mpz_clear(m);
}

/**
 * Multiply implements the Elliptical multiplication method described by
 * algorithm (7.2.4). This will return true if an illegal inversion occurred
 * during one of the Elliptical add method calls.
 */
bool Multiply(struct Point* theR, mpz_t& theM, struct Point& P, mpz_t& theN,
              mpz_t& theA)
{
  bool anResult = false;  // True if illegal inversion occurred

  // If theM provided is zero, return O (point at infinity)
  if(mpz_cmp_ui(theM, 0) == 0)
  {
    // Return O since theN provided was 0
    mpz_set_ui(theR->x, 0);
    mpz_set_ui(theR->y, 1);
  }
  else
  {
    mpz_t i; // number of multiplies to perform
    mpz_t t; // Used for gcd tests which makes Multiply return faster
    Point A; // The original number provided
    Point B; // Another number to be added, starts at infinity

    // Initialize our counter
    mpz_init_set(i, theM);
    mpz_init(t);
    mpz_init_set(A.x, P.x);
    mpz_init_set(A.y, P.y);
    mpz_init_set_ui(B.x, 0);
    mpz_init_set_ui(B.y, 1);

    // Loop while no illegal inversions have occurred and theM is > 0
    while(false == anResult && mpz_cmp_ui(i, 0) > 0)
    {
      // When our counter is odd, use the Add method
      if(mpz_odd_p(i))
      {
        // Subtract one from our counter
        mpz_sub_ui(i, i, 1);

        // Compute difference between B.x and A.x
        mpz_sub(t, B.x, A.x);
        mpz_mod(t, t, theN);

        // Perform GCD test
        mpz_gcd(t, t, theN);

        // Continue our loop only if t == 1 or t == n
        anResult = !(mpz_cmp_ui(t, 1) == 0 || mpz_cmp(t, theN) == 0);

        // If A is at infinity (point O) don't do anything
        if(mpz_cmp_ui(A.x, 0) == 0 && mpz_cmp_ui(A.y, 1) == 0)
        {
          // do nothing
        }
        // If B is at infinity (point O) then B = A*O = A
        else if(mpz_cmp_ui(B.x, 0) == 0 && mpz_cmp_ui(B.y, 1) == 0)
        {
          mpz_set(B.x, A.x);
          mpz_set(B.y, A.y);
        }
        // Otherwise, just do the addition
        else
        {
          // Check for illegal inversions during Add
          anResult = Add(&B, A, B, theN);
        }
      }
      // Our counter is even, use the Double method instead of Add
      else
      {
        // Divide our counter by two for each double we perform
        mpz_tdiv_q_ui(i, i, 2);

        // Compute the value 2*A.y
        mpz_mul_ui(t, A.y, 2);
        mpz_mod(t, t, theN);

        // Perform GCD test
        mpz_gcd(t, t, theN);

        // Continue our loop only if t == 1 or t == n
        anResult = !(mpz_cmp_ui(t, 1) == 0 || mpz_cmp(t, theN) == 0);

        // Use the double method
        Double(&A, A, theN, theA);
      }
    }

    // B has our results, make sure they are not bigger than n
    mpz_mod(B.x, B.x, theN);
    mpz_mod(B.y, B.y, theN);

    // Set our results to whatever B is set to now
    mpz_set(theR->x, B.x);
    mpz_set(theR->y, B.y);

    // Clear our values used above
    mpz_clear(B.y);
    mpz_clear(B.x);
    mpz_clear(A.y);
    mpz_clear(A.x);
    mpz_clear(t);
    mpz_clear(i);
  }

  // Return true if an illegal inversions occurred or a divisor was found
  return anResult;
}

/**
 * EvaluatePoint will compute the multiple U = [m/q]P. Based on these results
 * N will be either composite or Q << N will need to be proven prime to prove
 * that N is prime. This will return -1 if N is found to be composite, 0 if
 * another point should be tested, or 1 if N is found to be prime if Q is
 * proven to be prime during the next iteration.
 */
int EvaluatePoint(Point* theU, Point* theV, Point& P, mpz_t& theN, mpz_t& theM,
                  mpz_t& theQ, mpz_t& theA, mpz_t& theB)
{
  int anResult = 0;  // Assume that we will need another point
  bool anComposite = false;  // True if Multiply returns an illegal inversion
  mpz_t t;  // Store the result of m/q for computing U

  // Initialize our temporary values
  mpz_init(t);

  // First compute t = m/q
  mpz_tdiv_q(t, theM, theQ);

  // Now compute U = [m/q]P
  anComposite = Multiply(theU, t, P, theN, theA);

  // Clear our temporary value
  mpz_clear(t);

  // Did we have an illegal inversion?
  if(anComposite)
  {
    anResult = -1; // Illegal inversion occurred, N is composite
  }
  // Make sure U != O (infinity point) before calculating V
  else if(mpz_cmp_ui(theU->x,0) == 0 && mpz_cmp_ui(theU->y,1) == 0)
  {
    anResult = 0; // Try another point, this point didn't tell us anything
  }
  else
  {
    // Now compute V = [q]U
    anComposite = Multiply(theV, theQ, *theU, theN, theA);

    // If V gives us an illegal inversion, then we now N is prime if Q is
    // proven prime
    if(anComposite)
    {
      anResult = 1; // N is prime if Q is prime
    }
    // If V == O (infinity point) then N is prime if Q is prime
    //else if(mpz_cmp_ui(theV->x,0) == 0 && mpz_cmp_ui(theV->y,1) == 0)
    //{
    //  anResult = -1; // N is composite
    //}
    else
    {
      anResult = 0; // Try another point, this point didn't tell us anything
    }
  }

  // Return N is composite, try another point, or N is prime if Q is prime
  return anResult;
}

/**
 * AtkinMorain implements the original Atkin-Morain ECPP algorithm for proving
 * N is prime.  It calls various functions above to perform the ECPP algorithm.
 */
bool AtkinMorain(mpz_t& theN)
{
  bool anResult = false; // True if theNumber is proven prime, false otherwise
  bool anDone = false;  // True if theN was found to be prime or composite
  unsigned char anIndexD = 0;  // Index to selected discriminant
  unsigned int anIterations = 0;  // Number of iterations performed so far

  mpz_t n;  // The n to be tested
  mpz_t m;  // Curve order m
  mpz_t q;  // Factor q that if proven prime means that N is prime
  mpz_t u;  // Solution u
  mpz_t v;  // Solution v
  mpz_t g;  // Quadratic nonresidue g
  mpz_t a;  // Root a
  mpz_t b;  // Root b
  mpz_t Q;  // Q in ChoosePoint step
  Point P;  // Point P used in step 4: ChoosePoint
  Point U;  // Point U used in step 5: EvaluatePoint
  Point V;  // Point V used in step 5: EvaluatePoint
  mpz_t t;  // Temporary variable for testing y^2 mod n != Q

  // Initialize n, m, q, u, and v values
  mpz_init_set(n, theN);
  mpz_init(m);
  mpz_init(q);
  mpz_init(u);
  mpz_init(v);
  mpz_init(g);
  mpz_init(a);
  mpz_init(b);
  mpz_init(P.x);
  mpz_init(P.y);
  mpz_init(U.x);
  mpz_init(U.y);
  mpz_init(V.x);
  mpz_init(V.y);
  mpz_init(t);

  // Loop through each discriminant in the gD array or until our anDone value
  // is set to true.
  while(!anDone && anIndexD+1 < MAX_DISCRIMINANTS)
  {
    // Step 0: Use Miller-Rabin to test if theN is composite since there is
    // no guarrentee that ECPP will successfully find a u and v in Step 1, but
    // Miller-Rabin guarrentees to find all composites quickly.
    int anMillerRabin = mpz_probab_prime_p(n, 10);

    // Did Miller-Rabin prove n is composite?
    if(0 == anMillerRabin)
    {
      printf("Miller-Rabin says n is composite!\n");
      // N is composite
      anResult = false;

      // We are done
      anDone = true;

      // Exit the discriminant loop
      break;
    }
    // Did Miller-Rabin prove n is prime?
    else if(2 == anMillerRabin)
    {
      // N is proven prime
      anResult = true;
      
      // We are done
      anDone = true;
      
      // Exit the discriminant loop
      break;
    }
    else
    {
      // Continue with Atkin-Morain algorithm to prove n is prime or composite
    }

    // Step 1: ChooseDiscriminant will attempt to find a discriminant D that
    // satisfies steps 1a, 1b, and step 2 by incrementing through each
    // discriminant in our gD array, note that the gD array has a dummy entry
    // 0 at the beginning so we increment first before testing
    anIndexD++;

    // Step 1a: Find a discriminant that yields a Jacobi(D,N) == 1
    if(1 != mpz_jacobi(gD[anIndexD],n))
      continue; // Jacobi returned -1 or 0, try another discriminant

    // Step 1b: Find a u and v that satisfies 4n = u^2 + ABS(D)v^2 using the
    // modified Cornacchia algorithm (2.3.13)
    if(!ModifiedCornacchia(&u, &v, n, gD[anIndexD]))
      continue; // No solution found, u and v not set, try another discriminant

    // Step 2/3: FactorOrders attempts to find a possible order m that factors
    // as m = kq where k > 1 and q is a probable prime > (n^0.25 + 1)^2. If
    // this can't be done after K_max iterations than return FALSE and choose
    // a new discriminant D and curve m.
    if(!FactorOrders(&m, &q, u, v, n, gD[anIndexD]))
      continue; // No factor q or curve m was found, try another discriminant
    gmp_printf("Steps 1-3: d=%Zd, u=%Zd, v=%Zd, m=%Zd, q=%Zd\n", 
      gD[anIndexD], u, v, m, q);

    // Step 4a: CalculateNonresidue will find a random quadratic nonresidue
    // g mod p and if D=-3 a noncube g^3 mod p for use in step 4b
    CalculateNonresidue(&g, n, gD[anIndexD]);
    gmp_printf("G=%Zd\n", g);

    // Now that we have selected a curve m with factor q to be proven, obtain
    // curve parameters and test up to MAX_POINTS to see if N is composite
    unsigned int points = 0;
    do
    {
      // Step 4b: ObtainCurveParameters will attempt to obtain the curve
      // parameters a and b for an elliptic curve that would have order m if n
      // is indeed prime.

      // Loop through this curves valid a and b values (k iterates over them)
      unsigned int k = 0;
      while(!anDone && ObtainCurveParameters(&a, &b, n, gD[anIndexD], g, k))
      {
        gmp_printf("Step 4: a=%Zd, b=%Zd\n", a, b);

        // Step 5: ChoosePoint will try and find a point (x,y) on the curve
        // using the a and b values provided from above.
        anDone = ChoosePoint(&P,n,a,b);
        gmp_printf("Step 5: P(%Zd,%Zd)\n", P.x, P.y);

        // Did we find that N is composite while choosing a point?
        if(anDone)
        {
          printf("ChoosePoint returned true!\n");
          // N is composite
          anResult = false;

          // Exit the ObtainCurveParameters and Points loops
          break;
        }

        // Step 6: EvaluatePoint will compute the multiple U = [m/q]P. Based
        // on these results N will be either composite or Q << N will need to
        // be proven prime to prove that N is prime.
        int anTest = EvaluatePoint(&U, &V, P, n, m, q, a, b);

        // Did EvaluatePoint determine N was composite?
        if(anTest < 0)
        {
          printf("EvaluatePoint returned < 0!\n");
          // N is composite
          anResult = false;

          // We are done
          anDone = true;

          // Exit the ObtainCurveParameters and Points loops
          break;
        }
        // Did EvaluatePoint determine N was prime if Q is prime?
        else if(anTest > 0)
        {
          // N is prime if q can be proven prime
          anResult = true;

          // Print certificate information
          gmp_printf("n[%d]=%Zd\n", anIterations++, n);
          gmp_printf("d=%Zd\n", gD[anIndexD]);
          gmp_printf("u=%Zd\n", u);
          gmp_printf("v=%Zd\n", v);
          gmp_printf("m=%Zd q=%Zd\n", m, q);
          gmp_printf("a=%Zd b=%Zd\n", a, b);
          gmp_printf("P(%Zd,%Zd)\n", P.x, P.y);
          gmp_printf("U(%Zd,%Zd)\n", U.x, U.y);
          gmp_printf("V(%Zd,%Zd)\n", V.x, V.y);

          // Set n = q and start over
          mpz_set(n, q);

          // Reset our discriminant choice back to 0 and loop again
          anIndexD = 0;

          // Exit the points loop
          points = MAX_POINTS;

          // Exit the ObtainCurveParameters loop
          break;
        }
        // Does EvaluatePoint want us to try another point?
        else
        {
          // Try another a, b, and P
        }

        // Increment our k value
        k++;
      } // while(!anDone && ObtainCurveParameters(&a,&b,...))

      // Increment points for every k (a and b) tested
      points += k;
    } while(!anDone && points < MAX_POINTS);
  }
  
  // Warn the user that we ran out of discriminants before finding answer
  if(false == anDone && anIndexD+1 == MAX_DISCRIMINANTS)
  {
    printf("Atkin-Morain proof inconclusive, ran out of discriminants!\n");
  }

  // Clear our values used above
  mpz_clear(V.y);
  mpz_clear(V.x);
  mpz_clear(U.y);
  mpz_clear(U.x);
  mpz_clear(P.y);
  mpz_clear(P.x);
  mpz_clear(b);
  mpz_clear(a);
  mpz_clear(g);
  mpz_clear(v);
  mpz_clear(u);
  mpz_clear(q);
  mpz_clear(m);
  mpz_clear(n);

  // Return true if value is proven prime, false otherwise
  return anResult;
}

int main(int argc, char* argv[])
{
  mpz_t anNumber;  // Number to be tested for Primality
  
  // Initialize our random generator first
  gmp_randinit_default(gRandomState);

  // Seed our random generator using the clock
  gmp_randseed_ui(gRandomState, time(0));

  // Initialize our fixed list of discriminants
  InitDiscriminants();

  // TESTTESTTEST
  if(0)
  {
    // Initialize our number
    mpz_init_set_ui(anNumber, 1073748191);
    
    do
    {
      if(!AtkinMorain(anNumber))
      {
        gmp_printf("FAIL [%Zd]\n", anNumber);
      }
      else
      {
        //gmp_printf("PASS [%Zd]\n", anNumber);
      }
    
      mpz_nextprime(anNumber, anNumber);
    } while(true);
  }
  else
  {
    // Initialize our number
    mpz_init(anNumber);

    // Get the number to test
    gmp_scanf("%Zd", &anNumber);
  
    // Use Atkin-Morain to determine if the number is prime
    if(!AtkinMorain(anNumber))
    {
      printf("composite\n");
    }
    else
    {
      printf("proven prime\n");
    }
  }

  // Clear the number before exiting the program
  mpz_clear(anNumber);

  // Return 0
  return 0;
}
