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

// Global constants
const unsigned int MAX_DISCRIMINANTS = 28;

// Global Variables
gmp_randstate_t gRandomState;  ///< Holds random generator state and algorithm type
mpz_t gDiscriminants[MAX_DISCRIMINANTS];

// Point structure
struct Point
{
  mpz_t x;
  mpz_t y;
};

/**
 * InitDiscriminants will initialize the gDiscriminants array the appropriate
 * values.
 */
void InitDiscriminants(void)
{
  int32_t anD[MAX_DISCRIMINANTS] = {0,-3,-4,-7,-8,-11,-19,-43,-67,-163,-15,-20,-24,-35,-40,-51,-52,-88,-91,-115,-123,-148,-187,-232,-235,-267,-403,-427};

  // Call mpz_init on each array element and assign its value
  for(unsigned int i=0;i<MAX_DISCRIMINANTS;i++)
  {
    // Initialize and set each gDiscriminants array element
    mpz_init_set_si(gDiscriminants[i], anD[i]);

    // Display discriminants available
    //gmp_printf("%d=%Zd\n", i, gDiscriminants[i]);
  }
}


/**
 * SquareMod returns the solution x to x^2 === a (mod p) which is used by the
 * ModifiedCornacchia to find the initial square root.
 */
bool SquareMod(mpz_t* theX, mpz_t theA, mpz_t& theP)
{
  bool anResult = false; // Was a valid X found?
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
      mpz_powm(*theX, anBase, anExp2, theP);

      // Clear our values we don't need anymore
      mpz_clear(anBase);
      mpz_clear(anExp2);

      // Multiply the result by 2*a*x mod p
      mpz_mul_ui(*theX, *theX, 2);
      mpz_mul(*theX, *theX, theA);
      mpz_mod(*theX, *theX, theP);
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
    // TODO: Implement hardest case!
    printf("TODO: SquareMod needs case 2 of algorithm 2.3.8\n");
  }

  // Clear our temporary variables before returning
  mpz_clear(anMod8);
  mpz_clear(anMod4);

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
    if(dt != xt)
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
 * LenstraECM will attempt to find the largest non-trivial prime number that
 * will factor theN provided. This is called by FindFactor to find a probable
 * prime theQ that if proven to be prime will prove theN to be prime.
 */
void LenstraECM(mpz_t* theQ, mpz_t& theN)
{
  unsigned long b1 = 10000; // Limit to 10000 iterations
  Point P;  // Point (x,y) on curve E
  mpz_t a;  // Random value a for curve E
  mpz_t b;  // Random value b for curve E
  mpz_t g;  // The factor found if any
  mpz_t t;  // Temporary value for computing b

  // Initialize our values
  mpz_init(P.x);
  mpz_init(P.y);
  mpz_init(a);
  mpz_init(b);
  mpz_init(g);
  mpz_init(t);

  // Loop through b1 iterations trying to find a non-trivial factor of theN
  do {
    // Pick a random x from 0 to N-1
    mpz_urandomm(P.x, gRandomState, theN);
    // Pick a random y from 0 to N-1
    mpz_urandomm(P.y, gRandomState, theN);
    // Pick a random a from 0 to N-1
    mpz_urandomm(a, gRandomState, theN);

    // Compute b = (y^2 - x^3 - ax) mod n
    mpz_pow_ui(t, P.x, 3); // x^3
    mpz_submul(t, a, P.x); // x^3 - ax
    mpz_pow_ui(b, P.y, 2); // y^2
    mpz_sub(b, b, t);      // (y^2) - (x^3 - ax)
    mpz_mod(b, b, theN);   // (y^2 - x^3 - ax) mod n

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
  } while(true);

  // TODO: Finish LenstraECM implementation!

  // Clear our values used above
  mpz_clear(t);
  mpz_clear(g);
  mpz_clear(b);
  mpz_clear(a);
  mpz_clear(P.y);
  mpz_clear(P.x);
}

/**
 * FindFactor will attempt to reduce theM provided to a potential prime theQ
 * that is also less than theT by factoring all easy primes out of theM.  If
 * theM is a composite then FindFactor will return false meaning theQ wasn't found because m is a composite.
 */
bool FindFactor(mpz_t* theQ, mpz_t& theM, mpz_t& theT)
{
  bool anResult = true; // A suitable theQ was found
  unsigned long count = 1000000; // Number of Prime numbers to remove from theQ
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
  do {
    // Remove as many of this prime as possible
    mpz_remove(*theQ, *theQ, prime);

    // Find the next prime to try to remove from theQ
    mpz_nextprime(prime, prime);

    // Is theQ prime now? then exit our loop
    if(mpz_probab_prime_p(theM, 10))
      break;
  } while(mpz_cmp(*theQ, theT) >= 0 && count-- > 0);

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
bool FactorOrders(mpz_t* theM, mpz_t* theQ, mpz_t& theU, mpz_t& theV, mpz_t& theN, mpz_t& theD)
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
  else if(mpz_cmp_si(theD, -4) >= -4)
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
          // nothing was found, anResult should already be false
          //printf("FactorOrder failed!\n");
        }
      }

      // Clear our values used above
      mpz_clear(m4);
      mpz_clear(m3);
    }
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
 * ObtainCurveParameters will attempt to obtain the curve parameters a and b
 * for an elliptic curve that would have order m if n is indeed prime.
 */
bool ObtainCurveParameters(mpz_t* theA, mpz_t* theB, mpz_t& theN, mpz_t& theD)
{
  bool anResult = true; // True if curve parameters a and b were found

  // Special case for D=-3 and D=-4
  if(mpz_cmp_si(theD, -3) == 0)
  {
    // a=0, b=-1 mod n
    mpz_set_ui(*theA,0);
    mpz_set(*theB,theN);
    mpz_sub_ui(*theB, *theB, 1);
  }
  else if(mpz_cmp_si(theD, -4) == 0)
  {
    // a=-1 mod n, b=0
    mpz_set(*theA,theN);
    mpz_sub_ui(*theA, *theA, 1);
    mpz_set_ui(*theB,0);
  }
  else
  {
    // TODO: Implement algorithm 7.5.9 (manually compute Hilbert) OR
    //                 algorithm 7.5.10 (Hilbert lookup table)
    printf("TODO: ObtainCurveParameters: Implement algorithm 7.5.9 or 7.5.10\n");
    anResult = false;
  }

  // Return true if we found our curve parameters a and b above
  return anResult;
}

/**
 * Add implements the Elliptic add method described by Algorithm 7.2.2. This
 * will return true if an illegal inversion occurred and will not set theR
 * unless the inversion succeeded.
 */
bool Add(struct Point* theR, struct Point& theP1, struct Point& theP2, mpz_t& theN)
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
  if(anResult)
  {
    Point R;  // Point R to return via theR

    // Initialize our working point R
    mpz_init(R.x);
    mpz_init(R.y);
    
    // Compute (y2 - y1)
    mpz_mul(R.x, theP2.y, theP1.y);

    // Compute m = (y2-y1)(x2-x1)^-1
    mpz_mul(m, R.x, m);
    
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
 * Algorithm 7.2.2.
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
 * Algorithm 7.2.4. This will return true if an illegal inversion occurred
 * during one of the Elliptical add method calls.
 */
bool Multiply(struct Point* theR, mpz_t theM, struct Point& P)
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
    // Loop while no illegal inversions have occurred and theM is > 0
    while(false == anResult && mpz_cmp_ui(theM, 0) > 0)
    {
      if(mpz_odd_p(theM))
      {
        
      }
      else
      {
      }
    }
    // TODO: Perform the multiply using Add and Subtract methods
  }
}

/**
 * AtkinMorain implements the original Atkin-Morain ECPP algorithm for proving
 * N is prime.  It calls various functions above to perform the ECPP algorithm.
 */
bool AtkinMorain(mpz_t& theN)
{
  bool anResult = false; // True if theNumber is proven prime, false otherwise
  bool anFound = false;  // True if m, u, v, and q were found in step 1
  unsigned char anIndexD = 0;  // Index to selected discriminant
  mpz_t m;  // Curve order m
  mpz_t q;  // Factor q that if proven prime means that N is prime
  mpz_t u;  // Solution u
  mpz_t v;  // Solution v
  mpz_t a;  // Root a
  mpz_t b;  // Root b
  mpz_t Q;  // Q in ChoosePoint step
  Point P;  // Point P used in step 4: ChoosePoint
  Point U;  // Point U used in step 5: EvaluatePoint
  Point V;  // Point V used in step 5: EvaluatePoint
  mpz_t t;  // Temporary variable for testing y^2 mod n != Q

  // Step 0: Use Miller-Rabbin to test if theN is composite since there is no
  // guarrentee that ECPP will successfully find a u and v in Step 1, but
  // Miller-Rabbin guarrentees to find all composites quickly.
  if(mpz_probab_prime_p(theN, 10) == 0)
    return false;

  // Initialize m, q, u, and v values
  mpz_init(m);
  mpz_init(q);
  mpz_init(u);
  mpz_init(v);
  mpz_init(a);
  mpz_init(b);
  mpz_init(Q);
  mpz_init(P.x);
  mpz_init(P.y);
  mpz_init(U.x);
  mpz_init(U.y);
  mpz_init(V.x);
  mpz_init(V.y);
  mpz_init(t);

  // Step 1: ChooseDiscriminant will attempt to find a discriminant D that
  // satisfies the following:
  //   Jacobi(D/N) = 1
  //   4n = u^2 + ABS(D)v^2 where u and v can be found using D
  //   m curve order can be found

  // Find a valid discriminant, u, v, and m values
  do
  {
    // Increment D - note first pass will skip dummy 0 entry in anD array above
    anIndexD++;

    // Find Jacobi(D,N) that equals 1
    if(1 != mpz_jacobi(gDiscriminants[anIndexD],theN))
      continue; // Jacobi returned -1 or 0

    gmp_printf("d=%Zd\n", gDiscriminants[anIndexD]);

    // Try to retrieve u and v using modified Cornacchia algorithm (2.3.13)
    if(false == ModifiedCornacchia(&u, &v, theN, gDiscriminants[anIndexD]))
    {
      continue; // No solution found, u and v not set
    }
    else
    {
      gmp_printf("u=%Zd v=%Zd\n", u, v);
    }

    // Step 2: FactorOrders attempts to find a possible order m that factors
    // as m = kq where k > 1 and q is a probable prime > (n^0.25 + 1)^2. If
    // this can't be done after K_max iterations than return FALSE and choose
    // a new discriminant D and curve m.
    if(false == FactorOrders(&m, &q, u, v, theN, gDiscriminants[anIndexD]))
    {
      continue; // No factor q or curve m was found, try another D
    }
    else
    {
      gmp_printf("m=%Zd q=%Zd\n", m, q);
    }

    // Step 3: ObtainCurveParameters will attempt to obtain the curve
    // parameters a and b for an elliptic curve that would have order m if n
    // is indeed prime.
    if(false == ObtainCurveParameters(&a, &b, theN, gDiscriminants[anIndexD]))
    {
      continue; // No curve parameter a and b was obtained
    }
    else
    {
      gmp_printf("a=%Zd b=%Zd\n", a, b);

      // D, m, q, a and b were all found
      anFound = true;
      // Exit loop
      break;
    }
  } while(anIndexD + 1 < MAX_DISCRIMINANTS);

  // Try to find a point using the curve found above that proves N prime or
  // composite
  bool anDone = false;
  while(true == anFound && false == anDone)
  {
    // Step 4: ChoosePoint will try and find a point (x,y) on the curve
    // using the a and b values provided from above.
    do
    {
      // Step 4a: Choose a random x such that Q = (x^3 + ax + b) mod n and
      // Jacobi(Q/N) != -1.
      do
      {
        // Pick a random x from 0 to N-1
        mpz_urandomm(P.x, gRandomState, theN);

        // Compute Q = x^3
        mpz_powm_ui(Q, P.x, 3, theN);
        
        // Compute Q = Q + ax
        mpz_addmul(Q, a, P.x);
        
        // Compute Q = Q + b
        mpz_add(Q, Q, b);
        
        // Compute Q = Q mod n
        mpz_mod(Q, Q, theN);
      } while(-1 == mpz_jacobi(Q,theN));

      // Step 4b: Apply Algorithm 2.3.8 or 2.3.9 (with a = Q and p = n) to find
      // an integer y that would satisfy y2 = Q (mod n) if n were prime
      anFound = SquareMod(&P.y, Q, theN);
    } while(false == anFound);

    // Step 4c: If y^2 mod n != Q then N is composite
    mpz_powm_ui(t, P.y, 2, theN);
    if(mpz_cmp(Q, t) != 0)
    {
      // N is composite
      anResult = false;

      // We are done
      anDone = true;

      // Exit this loop
      break;
    }
    else
    {
      gmp_printf("P(%Zd,%Zd)\n", P.x, P.y);
    }

    // Step 5: EvaluatePoint will compute the multiple U = [m/q]P. Based on
    // these results N will be either composite or Q << N will need to be
    // proven prime to prove that N is prime.
    //anDone = EvaluatePoint(anResult, x, y);
    break;
  }

  // Clear our values used above
  mpz_clear(t);
  mpz_clear(V.y);
  mpz_clear(V.x);
  mpz_clear(U.y);
  mpz_clear(U.x);
  mpz_clear(P.y);
  mpz_clear(P.x);
  mpz_clear(Q);
  mpz_clear(b);
  mpz_clear(a);
  mpz_clear(v);
  mpz_clear(u);
  mpz_clear(q);
  mpz_clear(m);

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

  InitDiscriminants();

  mpz_init(anNumber);

  // Get the number to test
  //cin >> anNumber;

  // Try test number of 2^89 - 1 which is prime
  mpz_ui_pow_ui(anNumber, 2, 89);
  mpz_sub_ui(anNumber, anNumber, 1);

  // Show the number to be tested
  gmp_printf("n=%Zd\n", anNumber);

  if(false == AtkinMorain(anNumber))
  {
    printf("N is composite\n");
  }
  else
  {
    printf("N is prime\n");
  }

  // Clear the number before exiting the program
  mpz_clear(anNumber);

  // Return 0
  return 0;
}
