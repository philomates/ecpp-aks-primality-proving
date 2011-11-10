/*
  Author:  Pate Williams (c) 1997 & James Wanless (c) 2008-11 & Chris Studholme (c) 2004
 & Elisavet Konstantinou & Yiannis Stamatiu & Christos Zaroliagis (c) 2003
 & Scott Contini (c) 1996 & Paul Zimmermann (c) 2005 & Timo Poikola (c) 2010

  Multiple precision Atkin primality test.
  See "A Course in Computational Algebraic Number
  Theory" by Henri Cohen Algorithm 9.2.4 page 474.

  *   This file is part of the GMP-ECPP Math C++ Program.
  *
  *   The GMP-ECPP program is free software; you can redistribute it and/or
  *   modify it under the terms of the GNU General Public
  *   License as published by the Free Software Foundation; either
  *   version 3 of the License, or (at your option) any later version.
  *
  *   The GMP-ECPP program is distributed in the hope that it will be useful,
  *   but WITHOUT ANY WARRANTY; without even the implied warranty of
  *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  *   Lesser General Public License for more details.
  *
  *   You should have received a copy of the GNU General Public License
  *   along with this program.  If not, see <http://www.gnu.org/licenses/>.


*/

#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <sys/time.h>
#include <time.h>

#include <iostream>
#include <gmpxx.h>
#include <gmp.h>

using namespace std;

#define SIEVE_LIMIT 1073741824

#define LOG2 0.301029995f

#define POLY_SIZE 8192l

#define PRECISION 10000l

#define ERROR_SHIFT 1000l
#define ITER 1000l
#define BMAX 2000l
#define DMAX 20l

#define HAVE_GETTIMEOFDAY   1

gmp_randstate_t rstate;


int quiet = 0;
int Quiet = 0;
int verbose = 0;
int one = 0;
int weber = 0;
long seed = 0;
long precision = 0;
long Bmax = 0;
long Dmax = 0;
long error_shift = 0;

bool staticBmax = false;
bool staticDmax = false;

struct complex {mpf_class x, y;};
struct point {mpz_class x, y;};

mpf_class PIf(3, PRECISION);
mpf_class Ef(2, PRECISION);
mpf_class NATLOGONEPOINTNINEf(1, PRECISION);
mpf_class ONEPOINTNINE(1.9, PRECISION);
mpf_class ZERO(0, PRECISION);
mpf_class ONE(1, PRECISION);
mpf_class TWO(2, PRECISION);
mpf_class THREE(3, PRECISION);
mpf_class FOUR(4, PRECISION);
mpf_class FIVE(5, PRECISION);
mpf_class EIGHT(8, PRECISION);
mpf_class TWELVE(12, PRECISION);
mpf_class EIGHTEEN(18, PRECISION);
mpf_class FIFTYSEVEN(57, PRECISION);
mpf_class TWOHUNDREDTHIRTYNINE(239, PRECISION);

mpz_class rand2()
/* returns GMP pseudo-random number of 32-bits */
{
  mpz_t temp;
  mpz_init(temp);
  mpz_urandomb(temp, rstate, 32L);
  mpz_class temp_class(temp);
  mpz_clear(temp);
  return temp_class;
}

mpz_class gcd(mpz_class a, mpz_class b)
/* returns greatest common divisor of a and b */
{
  mpz_t temp;
  mpz_init(temp);
  mpz_gcd(temp, a.get_mpz_t(), b.get_mpz_t());
  mpz_class temp_class(temp);
  mpz_clear(temp);
  return temp_class;
}

mpz_class inverse(mpz_class a, mpz_class b)
/* returns inverse of a modulo b or 0 if it does not exist */
{
  mpz_t temp;
  mpz_init(temp);
  if (!mpz_invert(temp, a.get_mpz_t(), b.get_mpz_t()))
    mpz_set_si(temp, 0);
  mpz_class temp_class(temp);
  mpz_clear(temp);
  return temp_class;
}

mpz_class modpos(mpz_class a, mpz_class b)
/* returns a modulo b, strictly non-negative */
{
  mpz_class temp_class;
  temp_class = a % b;
  if (temp_class < 0)
    temp_class += b;
  return temp_class;
}

mpz_class exp_mod(mpz_class x, mpz_class b, mpz_class n)
/* returns x ^ b mod n */
{
  mpz_t temp;
  mpz_init (temp);
  mpz_powm(temp, x.get_mpz_t(), b.get_mpz_t(), n.get_mpz_t());
  mpz_class temp_class(temp);
  mpz_clear (temp);
  return temp_class;
}

mpz_class exp_pow(mpz_class a, long e)
/* returns a ^ e */
{
  mpz_t temp;
  mpz_init(temp);
  mpz_pow_ui(temp, a.get_mpz_t(), e);
  mpz_class temp_class(temp);
  mpz_clear(temp);
  return temp_class;
}

int JACOBI(mpz_class a, mpz_class b)
{
    return mpz_jacobi(a.get_mpz_t(), b.get_mpz_t());
}

int Rabin_Miller(mpz_class n)
/* given an integer n >= 3 returns 0 composite, 1 probably prime */
{
  return mpz_probab_prime_p(n.get_mpz_t(), 10);
}

/*
 Given a positive integer n, this algorithm
 determines whether n is a square or not, and if it is
 outputs the square root of n.
 */
int square_test(mpz_class n, mpz_class *root)
{
  int issquare = 0;

    *root = sqrt(n);
    if ((*root * *root) == n) issquare = 1;

  return issquare;
}

/*
 Given a positive integer n, this algorithm
 determines whether n is a cube or not, and if it is
 outputs the cube root of n.
 */
int cube_test(mpz_class n, mpz_class *root)
{
  int iscube = 0;

  mpz_t temp;
  mpz_init(temp);
  iscube = mpz_root(temp, n.get_mpz_t(), 3l);
  mpz_class temp_class(temp);
  mpz_clear(temp);

  *root = temp_class;

  return iscube;
}

mpz_class nextp(mpz_class n)
/* returns next prime after n */
{
  mpz_t temp;
  mpz_init(temp);
  mpz_nextprime(temp, n.get_mpz_t());
  mpz_class temp_class(temp);
  mpz_clear(temp);
  return temp_class;
}

int addition_1(mpz_class n, struct point P1, struct point P2, struct point *P3)
/* affine coords */
/* returns 1 if P1 = -P2 therefore P1 + P2 = O, 0 otherwise */
/* P1 != P2 */
{
  mpz_class delta_x;
  mpz_class delta_y;
  mpz_class m;

  delta_x = modpos(P2.x - P1.x, n);
  delta_y = modpos(P2.y - P1.y, n);

  if (P1.x == P2.x && modpos(P1.y + P2.y, n) == 0) {
    P3->x = 0, P3->y = 1;
    return 1;
  }

  /* calculate m = (y2 - y1)(x2 - x1) ^ -1 mod n */
  m = modpos(delta_y * inverse(delta_x, n), n);

  /* calculate x3 = m ^ 2 - (x1 + x2) mod n */
  P3->x = modpos(m * m - (P1.x + P2.x), n);

  /* calculate y3 = m(x1 - x3) - y1 mod n */
  P3->y = modpos(m * (P1.x - P3->x) - P1.y, n);

  return 0;
}

void addition_2(mpz_class a, mpz_class n, struct point P1, struct point *P3)
/* affine coords */
/* P1 == P2 */
{
  mpz_class m;

  /* calculate m = (3x1 ^ 2 + a)(2y1) ^ -1 mod n */
  m = modpos((3 * P1.x * P1.x + a) * inverse(2 * P1.y, n), n);

  /* calculate x3 = m ^ 2 - 2x1 mod n */
  P3->x = modpos(m * m - 2 * P1.x, n);

  /* calculate y3 = m(x1 - x3) - y1 mod n */
  P3->y = modpos(m * (P1.x - P3->x) - P1.y, n);
}

int multiply(mpz_class a, mpz_class k, mpz_class n, struct point P, struct point *R, mpz_class *d)
/* binary ladder */
/* returns -1 if O encountered, 0 if divisor not found, 1 otherwise */
{
  int value = 1;
  struct point A, B, C;

  /*  A = P */
  A = P;
  /* B = O = (0, 1) the point at infinity */
  B.x = 0, B.y = 1;

  while (value && k > 0) {
    if (k % 2 != 0) {
      *d = gcd(modpos(B.x - A.x, n), n);

      k--;
      value = (*d == 1 || *d == n);

      if (A.x == 0 && A.y == 1);
      else if (B.x == 0 && B.y == 1) B = A;
      else if (value) {
        addition_1(n, A, B, &C);
        B = C;
      }
    }
    else {
      *d = gcd(modpos(2 * A.y, n), n);

      k >>= 1;
      value = (*d == 1 || *d == n);

      if (value) {
        addition_2(a, n, A, &C);
        A = C;
      }
    }
  }

  *R = B;
  R->x = modpos(R->x, n);
  R->y = modpos(R->y, n);

  if (R->x == 0 && R->y == 1) return -1;

  return !value;
}

/*
 Algorithm 1.5.1 (Square Root Mod p). See "A Course
 in Computational Algebraic Number Theory" by Henri
 Cohen page 33. This algorithm either outputs an
 x such that x^2 == a mod p, or says that such a
 solution does not exist.
 */
mpz_class square_root_mod1(mpz_class a, mpz_class p)
{
  mpz_class b, e = 0, m, n, q = p - 1, r, t, x, y, z;
  mpz_class TWO = 2;
  mpz_class A;

  A = modpos(a, p);

  /* p - 1 = 2 ^ e * q with q odd */
  while (q % 2 == 0) q >>= 1, e++;

  /* find generator */
  do
    do n = modpos(rand2(), p); while (n == 0);
  while (JACOBI(n, p) != -1);

  z = exp_mod(n, q, p);
  y = z, r = e;
  x = exp_mod(a, (q - 1) / 2, p);

  b = modpos((A * modpos(x, p)) * x, p);

  x = modpos((A * x), p);

  do {
    if (b == 1) return x;

    m = 1;
    while (exp_mod(b, TWO ^ m, p) != 1) m++;
    if (m == r) return 0;

    t = exp_mod(y, TWO ^ (r - m - 1), p);
    y = modpos(t * t, p);

    r = m;
    x = modpos(x * t, p);

    b = modpos(b * y, p);

  } while (true);
}

/* Function: sqrtmod
 Description: Calculates the squareroot of a mod prime p.
 Uses the RESSOL algorithm of Shanks.  For description of algorithm,
 see:  An Introduction to the Theory of Numbers, 5th edition.
 Authors: Niven, Zuckerman, and Montgomery.  Pages 110-114.
 Always return the smallest squareroot
 */
/******************************************************************
 * computes x so that x^2 == a mod p for prime p, and puts x in s.
 * if no such x exists or if p is not prime, then sets s=0.
 ******************************************************************/
void
mpz_sqrtmod (mpz_t s, mpz_t a, mpz_t p)
{

  mpz_t q, nn, x, y, t, b, z;
  register long n, i, m, r;

  mpz_init (q);
  mpz_init (y);

  if (mpz_fdiv_ui (p, 4) == 3)
  {
    mpz_add_ui (q, p, 1);
    mpz_div_2exp (q, q, 2);
    mpz_powm (s, a, q, p);
    mpz_mul (y, s, s);
    mpz_mod (y, y, p);
    if (mpz_cmp (y, a))
      mpz_set_ui (s, 0);
    mpz_clear (q);
    mpz_clear (y);
    return;
  }

  mpz_init (nn);
  mpz_init (x);
  mpz_init (t);
  mpz_init (b);
  mpz_init (z);

  mpz_sub_ui (q, p, 1);
  if (mpz_cmp_ui (q, 0) == 0)
    r = -1;
  else
  {
    r = mpz_scan1 (q, 0);
    mpz_div_2exp (q, q, r);
  }
  n = 3;
  if ((mpz_fdiv_ui (p, 4) & n) == 3)
    m = -mpz_kronecker_si (p, n);
  else
    m = mpz_kronecker_si (p, n);

  while (m != -1) {
  loop:
    n += 2;
    if (!(n % 3))
      n += 2;
    for (i = 5; i * i <= n; )  {
      if (!(n % i))
        goto loop;
      i += 2;
      if (!(n % i))
        goto loop;
      i += 4;
    }
    if ((mpz_fdiv_ui (p, 4) & n) == 3)
      m = -mpz_kronecker_si (p, n);
    else
      m = mpz_kronecker_si (p, n);
  }

  mpz_set_si (nn, n);

  mpz_powm (z, nn, q, p);

  mpz_set (y, z);
  mpz_sub_ui (t, q, 1);
  mpz_div_2exp (t, t, 1); /* t = (q-1)/2 */
  mpz_powm (x, a, t, p);
  mpz_mul (t, x, x);
  mpz_mod (t, t, p);
  mpz_mul (b, t, a);
  mpz_mod (b, b, p);      /* b = a*x^2 mod p */
  mpz_mul (x, a, x);
  mpz_mod (x, x, p);

  while (1) {
    if (mpz_cmp_ui (b, 1) == 0)
    {
      mpz_set (s, x);
      goto end;
    }

    mpz_mul (t, b, b);
    mpz_mod (t, t, p);
    m = 1;
    while (mpz_cmp_ui (t, 1) != 0)
    {
      if (m == r)
      {
        mpz_set_ui (s, 0);
        goto end;
      }
      mpz_mul (t, t, t);
      mpz_mod (t, t, p);
      m ++;
    }

    mpz_set (t, y);
    if (r <= m)
    {
      mpz_set_ui (s, 0);
      goto end;
    }
    for (i = r - m - 1; i; i--)
    {
      mpz_mul (t, t, t);
      mpz_mod (t, t, p);
    }
    mpz_mul (y, t, t);
    mpz_mod (y, y, p);
    r = m;
    mpz_mul (x, x, t);
    mpz_mod (x, x, p);
    mpz_mul (b, b, y);
    mpz_mod (b, b, p);
  }
end:
  mpz_clear (q);
  mpz_clear (nn);
  mpz_clear (x);
  mpz_clear (y);
  mpz_clear (t);
  mpz_clear (b);
  mpz_clear (z);
}

mpz_class square_root_mod3(mpz_class a, mpz_class p)
/* returns square root of a modulo p in the easy cases, 0 otherwise */
{
  if (modpos(p, 4) == 3)
    return exp_mod(a, (p + 1) / 4, p);

  else if (modpos(p, 8) == 5) {
    if (modpos(exp_mod(a, (p - 1) / 4, p), p) == 1)
      return exp_mod(a, (p + 3) / 8, p);
    else
      return modpos(2 * a * exp_mod(4 * a, (p - 5) / 8, p), p);
  }

  else
    return 0;
}

mpz_class square_root_mod(mpz_class a, mpz_class p)
/* returns square root of a modulo p if it exists, 0 otherwise */
// version 3
{
  mpz_class temp_class3;
  if ((temp_class3 = square_root_mod3(a, p)) > 0)
    return temp_class3;

  mpz_class temp_class1;
  if ((temp_class1 = square_root_mod1(a, p)) > 0)
    return temp_class1;

  mpz_t temp;
  mpz_init (temp);
  mpz_sqrtmod(temp, a.get_mpz_t(), p.get_mpz_t());
  mpz_clear (temp);
  mpz_class temp_class2(temp);
  return temp_class2;
}

long max(long x, long y)
{
  if (x > y)
    return x;
  else
    return y;
}

/*
void cprint(struct complex u)
{
  cout << "(";
  print_mpf(u.x);
  cout << ", ";
  print_mpf(u.y);
  cout << ")\n";
  fflush(stdout);
}
*/

void print_mpf (mpf_class x)
{
  mpf_set_default_prec(precision);

  gmp_printf("%.*Ff", (int)(LOG2 * (x.get_prec())), x.get_mpf_t());

  return;
}

mpf_class sinf2(mpf_class x)
{
  mpz_class i, fact_2n_1 = 1;
  mpf_class x_2n_1(x, precision), sin_x(0, precision), old_sin_x(1, precision);

  bool sign = true;

  if (x < 0)
    return -sinf2(-x);

  if (x >= 2 * PIf)
    return sinf2(x - 2 * PIf * floor(x / (2 * PIf)));

  if (x >= PIf)
    return -sinf2(2 * PIf - x);

  if (x >= PIf / 2)
    return sinf2(PIf - x);


  for (i = 1; old_sin_x != sin_x; i += 2) {
    old_sin_x = sin_x;

    if (sign)
      sin_x += x_2n_1 / fact_2n_1;
    else
      sin_x -= x_2n_1 / fact_2n_1;
    x_2n_1 *= x * x;
    fact_2n_1 *= (i + 1) * (i + 2);
    sign = !sign;
  }

    return sin_x;
}

mpf_class cosf2(mpf_class x)
{
  mpz_class i, fact_2n = 1;
  mpf_class x_2n(1, precision), cos_x(0, precision), old_cos_x(1, precision);

  bool sign = true;

  if (x < 0)
    return cosf2(-x);

  if (x >= 2 * PIf)
    return cosf2(x - 2 * PIf * floor(x / (2 * PIf)));

  if (x >= PIf)
    return cosf2(2 * PIf - x);

  if (x >= PIf / 2)
    return -cosf2(PIf - x);


  for (i = 0; old_cos_x != cos_x; i += 2) {
    old_cos_x = cos_x;

    if (sign)
      cos_x += x_2n / fact_2n;
    else
      cos_x -= x_2n / fact_2n;
    x_2n *= x * x;
    fact_2n *= (i + 1) * (i + 2);
    sign = !sign;
  }

    return cos_x;
}

mpf_class atanf2(mpf_class x)
{
  mpz_class i, twon_1 = 1;
  mpf_class x_2n_1(x, precision), atan_x(0, precision), old_atan_x(1, precision);

  bool sign = true;

  if (x == 1)
    return (PIf / FOUR);

  if (x == -1)
    return (-PIf / FOUR);

  if (x > 1)
    return (PIf / TWO - atanf2(ONE / x));

  if (x < -1)
    return (-PIf / TWO - atanf2(ONE / x));

  for (i = 1; old_atan_x != atan_x; i += 2) {

    if (i > 1000000)
      break;

    old_atan_x = atan_x;

    if (sign)
      atan_x += x_2n_1 / twon_1;
    else
      atan_x -= x_2n_1 / twon_1;
    x_2n_1 *= x * x;
    twon_1 += 2;
    sign = !sign;
  }

    return atan_x;
}

mpf_class atan2f2(mpf_class x, mpf_class y)
{

  mpf_class r(0, precision);

  if (x == 0 && y == 0)
    return 0;
  if (abs(y) >= abs(x)) {
    r = atanf2(x / y);
    if (y < 0) {
      if (x >= 0) return (r += PIf);
      return (r -= PIf);
    }
  }
  else {
    r = -atanf2(y / x);
    if (x < 0) return (r -= PIf / TWO);
    return (r += PIf / TWO);
  }

  return r;
}

mpf_class logf2(mpf_class x)
{
  mpz_class i, n = 1;

  bool sign = true;

  if (x <= 0)
    return (ZERO);

  if (x == 1)
    return (ZERO);

  if (x > ONEPOINTNINE)
    return (logf2(x / ONEPOINTNINE) + NATLOGONEPOINTNINEf);

  if (x < 1)
    return (logf2(x * Ef) - ONE);

  x--;

  mpf_class x_n(x, precision), log_x(0, precision), old_log_x(1, precision);


  for (i = 1; old_log_x != log_x; i += 1) {
    old_log_x = log_x;

    if (sign)
      log_x += x_n / n;
    else
      log_x -= x_n / n;
    x_n *= x;
    n++;
    sign = !sign;
  }

    return log_x;
}

mpf_class expf2(mpf_class x)
{
  mpz_class i, n = 1;
  mpf_class x_n(x, precision), exp_x(0, precision), old_exp_x(1, precision);

  if (x > 1)
    return (expf2(x - ONE) * Ef);

  if (x < 0)
    return (expf2(x + ONE) / Ef);

  for (i = 1; old_exp_x != exp_x; i += 1) {
    old_exp_x = exp_x;

    exp_x += x_n / n;
    x_n *= x;
    n *= (i + 1);
  }

    return ONE + exp_x;
}

mpf_class sinhf2(mpf_class x)
{
  return ((expf2(x) - expf2(-x)) / TWO);
}

mpf_class coshf2(mpf_class x)
{
  return ((expf2(x) + expf2(-x)) / TWO);
}

mpf_class cabs2(struct complex u)
{
  return sqrt(u.x * u.x + u.y * u.y);
}

struct complex cadd(struct complex u, struct complex v)
{
  struct complex w;

  w.x = u.x + v.x;
  w.y = u.y + v.y;
  return w;
}

struct complex csub(struct complex u, struct complex v)
{
  struct complex w;

  w.x = u.x - v.x;
  w.y = u.y - v.y;
  return w;
}

struct complex cdiv(struct complex u, struct complex v)
{
  mpf_class temp1, temp2;
  struct complex w;

  if (abs(v.x) <= abs(v.y)) {
    temp1 = v.x / v.y;
    temp2 = v.y + temp1 * v.x;
    w.x = (temp1 * u.x + u.y) / temp2;
    w.y = (temp1 * u.y - u.x) / temp2;
  }
  else {
    temp1 = v.y / v.x;
    temp2 = v.x + temp1 * v.y;
    w.x = (u.x + temp1 * u.y) / temp2;
    w.y = (u.y - temp1 * u.x) / temp2;
  }
  return w;
}

struct complex cmul(struct complex u, struct complex v)
{
  struct complex w;

  w.x = u.x * v.x - u.y * v.y;
  w.y = u.x * v.y + u.y * v.x;
  return w;
}

struct complex cconj(struct complex u)
{
  struct complex w;

  w.x = u.x;
  w.y = - u.y;
  return w;
}

struct complex csqrt2(struct complex u)
{
  mpf_class r = sqrt(cabs2(u));
  mpf_class theta = atan2f2(u.y, u.x), phi = 0.5 * theta;
  struct complex w;

  w.x = r * cosf2(phi);
  w.y = r * sinf2(phi);
  return w;
}

struct complex cexp2(struct complex u)
{
  mpf_class e = expf2(u.x);
  struct complex w;

  w.x = e * cosf2(u.y);
  w.y = e * sinf2(u.y);
  return w;
}

struct complex csin2(struct complex u)
{
  struct complex w;

  w.x = sinf2(u.x) * coshf2(u.y);
  w.y = cosf2(u.x) * sinhf2(u.y);
  return w;
}

struct complex ccos2(struct complex u)
{
  struct complex w;

  w.x = cosf2(u.x) * coshf2(u.y);
  w.y = - sinf2(u.x) * sinhf2(u.y);
  return w;
}

struct complex clog2(struct complex u)
{
  mpf_class r = cabs2(u);
  struct complex w;

  w.x = logf2(r);
  w.y = atan2f2(u.y, u.x);
  return w;
}

struct complex cpow2(struct complex u, struct complex v)
{
  return cexp2(cmul(v, clog2(u)));
}

mpz_class round2(mpf_class x)
{
  if (x >= 0.0) return x + 0.5;
  return x - 0.5;
}

void cpoly_mul(long m, long n, struct complex *a,
               struct complex *b, struct complex *c,
               long *p)
{
  long i, j, k;
  struct complex ai, bj, sum, term;

  *p = m + n;
  for (k = 0; k <= *p; k++) {
    sum.x = sum.y = 0.0;
    for (i = 0; i <= k; i++) {
      j = k - i;
      if (i > m) ai.x = ai.y = 0.0;
      else ai = a[i];
      if (j > n) bj.x = bj.y = 0.0;
      else bj = b[j];
      term = cmul(ai, bj);
      sum = cadd(sum, term);
    }
    c[k] = sum;
  }
}

struct complex Sum(struct complex q)
{
  mpf_class sign = - 1.0;
  mpf_class error_term = 1;
  long e1, e2, n = 1;
  struct complex expon1 = {0.0, 0.0}, expon2 = {0.0, 0.0};
  struct complex sum = {0.0, 0.0}, term1, term2;

  for (long i = 0; i < error_shift; i++)
    error_term /= 2;

  /*
   cout << "error_term = ";
   print_mpf(error_term);
   cout << "\n";
   fflush(stdout);
   */

  do {
    e1 = n * (3 * n - 1) / 2;
    e2 = n * (3 * n + 1) / 2;
    expon1.x = e1;
    expon2.x = e2;
    term1 = cpow2(q, expon1);
    term2 = cpow2(q, expon2);
    term1 = cadd(term1, term2);
    term1.x *= sign;
    term1.y *= sign;
    sum = cadd(sum, term1);
    sign = - sign;
    n++;
  } while (cabs2(term1) > error_term);

  sum.x += 1.0;

  return sum;
}

struct complex Delta(struct complex q)
{
  struct complex s = {24.0, 0.0};
  struct complex Delta;

  Delta = cmul(q, cpow2(Sum(q), s));

  return Delta;
}

struct complex Theta(long D, long A, long B)
{
  struct complex s, Theta;

  s.x = -sqrt((double)D) * PIf / A;
  s.y = -B * PIf / A;

  Theta = cexp2(s);

//  cout << "Theta.x=" << Theta.x << ", Theta.y=" << Theta.y << "\n";

  return Theta;
}

struct complex j(struct complex tau)
{
  struct complex c = {0.0, 2.0 * PIf}, f, f1, q, q2,
  twofiftysix = {256.0, 0.0}, one = {1.0, 0.0}, three = {3.0, 0.0};

  q = cexp2(cmul(c, tau));
  q2 = cmul(q, q);
  f = cdiv(Delta(q2), Delta(q));

  f1 = cmul(twofiftysix, f);
  f1 = cadd(one, f1);
  f1 = cpow2(f1, three);

  return cdiv(f1, f);
}

/* Weber functions as listed "Explicit Construction of the Hilbert Class Fields of Imaginary
 Quadratic Fields by Integer Lattice Reduction" by Kaltofen-Yui pp. 180-1 and
 "Elliptic Curves and Primality Proving" by Atkin-Morain p. 7 */
struct complex F(long j, long D, long A, long B)
{
  struct complex theta, minustwentyfourth = {-1.0/24.0, 0.0}, twelth = {1.0/12.0, 0.0};
  struct complex theta12, theta24, minustheta, theta2;
  struct complex result1, result2, result;
  struct complex root_two;
  struct complex two = {2.0, 0.0};

  root_two = csqrt2(two);

//  cout << "j=" << j << ", D=" << D << ", A=" << A << ", B=" << B << "\n";

  theta = Theta(D, A, B);

  theta24 = cpow2(theta, minustwentyfourth);
  theta12 = cpow2(theta, twelth);
  theta2 = cmul(theta, theta);

  minustheta.x = -theta.x;
  minustheta.y = -theta.y;

  if (j == 0)
    result1 = Sum(minustheta);
  if (j == 1)
    result1 = Sum(theta);
  if (j == 2)
    result1 = Sum(cmul(cmul(theta, theta), cmul(theta, theta)));

  result2 = cdiv(result1, Sum(theta2));

  if (j == 0 || j == 1)
    result = cmul(theta24, result2);
  if (j == 2)
    result = cmul(root_two, cmul(theta12, result2));

  return result;
}


/* the following functions implement the routines described in IEEEP1363 Standard */
long G(long D)
{
  if (D % 3 == 0)
    return 3;
  else
    return 1;
}

long I(long D)
{
  long t1, t2;
  long result = 0;

  t1 = D % 8;
  t2 = D % 3;

  if (t1 == 1 || t1 == 2 || t1 == 6 || t1 == 7)
    result = 3;

  if (t1 == 3 && t2 != 0)
    result = 0;

  if (t1 == 3 && t2 == 0)
    result = 2;

  if (t1 == 5)
    result = 6;

  return result;
}

long J1(long A, long C)
{
  long t1, t2;
  long result = 0;

  t1 = A * C;
  t2 = t1 % 2;

  if (t2 == 1)
    result = 0;

  t2 = C % 2;
  if (t2 == 0)
    result = 1;

  t2 = A % 2;
  if (t2 == 0)
    result = 2;

  return result;
}

long K(long D)
{
  long t1;
  long result = 0;

  t1 = D % 8;

  if (t1 == 1 || t1 == 2 || t1 == 6)
    result = 2;
  if (t1 == 3 || t1 == 7)
    result = 1;
  if (t1 == 5)
    result = 4;

  return result;
}

long L(long D, long A, long C)
{
  long t1, t2, t3;
  long result = 0;

  t1 = D % 8;
  t2 = (A * C) % 2;

  t3 = C % 2;

  if (t2 == 1 || (t1 == 5 && t3 == 0))
    result = A * A * C + A - C;
  if ((t1 == 1 || t1 == 2 || t1 == 3 || t1 == 6 || t1 == 7) && t3 == 0)
    result = A + C + C - A * C * C;

  t3 = A % 2;

  if (t1 == 3 && t3 == 0)
    result = A - C + 5 * A * C * C;

  if ((t1 == 1 || t1 == 2 || t1 == 5 || t1 == 6 || t1 == 7) && t3 == 0)
    result = A - C - A * C * C;

  return result;
}

long M(long A, long C)
{
  long t1, t2, t3;
  long result = 0;

  t1 = A % 2;

  if (t1 == 1) {
    t2 = (A * A - 1) / 8;
    t3 = t2 % 2;

    if (t3 == 0)
      result = 1;

    if (t3 == 1)
      result = -1;
  }

  if (t1 == 0) {
    t2 = (C * C - 1) / 8;
    t3 = t2 % 2;

    if (t3 == 0)
      result = 1;

    if (t3 == 1)
      result = -1;
  }

  return result;
}


long N(long D, long A, long C)
{

  long t1, t2;
  long result = 0;

  t1 = D % 8;
  t2 = (A * C) % 2;

  if (t1 == 5 || (t1 == 3 && t2 == 1) || (t1 == 7 && t2 == 0))
    result = 1;
  if ((t1 == 1 || t1 == 2 || t1 == 6) || (t1 == 7 && t2 == 1))
    result = M(A, C);
  if (t1 == 3 && t2 == 0)
    result = -M(A, C);

  return result;
}


/* Weber switch as listed "Explicit Construction of the Hilbert Class Fields of Imaginary
 Quadratic Fields by Integer Lattice Reduction" by Kaltofen-Yui p. 178 and
 "Elliptic Curves and Primality Proving" by Atkin-Morain p. 17  and
 Weber "Lehrbuch der Algebra III" Section 127 */
struct complex u(long d, long a, long b, long c)
{
  struct complex result;
  struct complex f, temp1, temp2, temp3, temp4, temp5;
  long t1;

  long n, l, i, k, j, g;

  n = N(d, a, c);
  l = L(d, a, c);
  i = I(d);
  k = K(d);
  j = J1(a, c);
  g = G(d);
//  cout << "n=" << n << ", l=" << l << ", i=" << i << ", k=" << k << ", j=" << j << ", g=" << g << "\n";

  f = F(j, d, a, b);
//  cout << "f.x=" << f.x << ",f.y=" << f.y << "\n";

  temp1.x = k;
  temp1.y = 0.0;
  temp2 = cpow2(f, temp1);

  temp3.x = 0;
  temp3.y = -PIf/24.0;
  t1 = k * b * l;
  t1 = t1%48;
  if (t1 < 0) t1 += 48;
  temp3.y = temp3.y * t1;

  temp1 = cexp2(temp3);

  temp1.x *= n;
  temp1.y *= n;

  temp3.x = 2.0;
  temp3.y = 0.0;
  temp4.x = -i/6.0;
  temp4.y = 0.0;
  temp5 = cpow2(temp3, temp4);

  temp3 = cmul(temp1, temp5);
  temp4 = cmul(temp3, temp2);

  temp5.x = g;
  temp5.y = 0.0;

  result = cpow2(temp4, temp5);

/*
  if (d % 8 == -3)
    result = j(tau);
*/
  return result;
}

/*
 Algorithm 7.6.1 (Hilbert Class Polynomial). See "A Course
 in Computational Algebraic Number Theory" by Henri
 Cohen page 415. Given a negative discriminant D, this
 algorithm computes the monic polynomial of degree h(D)
 in Z[X] of which j((D+sqrt(D))/2) is a root.
 */
void Hilbert(long D, mpz_class *P, long *dP)
{
  mpz_class B;
  long i, t, a, b, c, d, dQ, dR;
  struct complex J, aa = {0.0, 0.0}, bb = {0.0, 0.0};
  struct complex one = {1.0, 0.0}, two = {2.0, 0.0};
  struct complex DD = {0.0, 0.0}, sqrtD, tau;
  struct complex PP[POLY_SIZE], Q[2], R[POLY_SIZE];

  d = D;

  DD.x = d;
  DD.y = 0.0;
  sqrtD = csqrt2(DD);

  PP[0].x = 1.0;
  PP[0].y = 0.0;
  *dP = 0;

  b = D % 2;
  if (b < 0) b += 2;
  B = sqrt(labs(D) / 3.0);
L2:
  t = (b * b - D) / 4;
  a = max(b, 1);
L3:
  if (t % a != 0) goto L4;

    aa.x = 2 * a;
    aa.y = 0.0;
    bb.x = - b;
    bb.y = 0.0;

    tau = cdiv(cadd(bb, sqrtD), aa);

    J = j(tau);

  if (a == b || a * a == t || b == 0) {
    Q[0].x = - J.x;
    Q[0].y = - J.y;
    Q[1].x = 1.0;
    Q[1].y = 0.0;
    dQ = 1;
  }
  else {
    Q[0].x = J.x * J.x + J.y * J.y;
    Q[0].y = 0.0;
    Q[1].x = - 2.0 * J.x;
    Q[1].y = 0.0;
    Q[2].x = 1.0;
    Q[2].y = 0.0;
    dQ = 2;
  }

  cpoly_mul(*dP, dQ, PP, Q, R, &dR);

  *dP = dR;
  for (i = 0; i <= dR; i++) PP[i] = R[i];
L4:
  a++;
  if (a * a <= t) goto L3;

  b += 2;
  if (b <= B) goto L2;

  for (i = 0; i <= *dP; i++) P[i] = round2(PP[i].x);
}

void Weber(long D, mpz_class *P, long *dP)
{
  long t1, t2, t3;
  mpz_class temp2;
  long i, a, b, c, A, B1, C, h1, help1, x2, temp1, dQ, dR;
  struct complex J;
  struct complex PP[POLY_SIZE], Q[2], R[POLY_SIZE];

  if (D % 4 == 0)
    D /= 4;

  help1 = (long) sqrt(labs(D) * 4.0 / 3.0);

  PP[0].x = 1.0;
  PP[0].y = 0.0;
  *dP = 0;

  b = 0;

  while (b <= help1) {
    t1 = b * b;
    t2 = -D * 4;
    t3 = t1 + t2;
    t1 = 4;

    t2 = t3 % t1;

    if (t2 == 0) {
      t1 = t3 / 4;

      if (b > 1)
        a = b;
      else
        a = 1;

      x2 = (long) sqrt((double)t1);

      while (a <= x2) {
        h1 = t1 % a;
        if (h1 == 0) {
          C = t1 / a;
          B1 = b;
          A = a;

          temp2 = gcd(gcd(A, B1), C);

          if (temp2 == 1) {
            temp1 = B1 / 2;

            J = u(-D, A, temp1, C);


            if (! ((B1 > 0) && (C > A) && (A > B1))) {
              Q[0].x = - J.x;
              Q[0].y = - J.y;
              Q[1].x = 1.0;
              Q[1].y = 0.0;
              dQ = 1;
            }
            else {
              Q[0].x = J.x * J.x + J.y * J.y;
              Q[0].y = 0.0;
              Q[1].x = - 2.0 * J.x;
              Q[1].y = 0.0;
              Q[2].x = 1.0;
              Q[2].y = 0.0;
              dQ = 2;
            }

            cpoly_mul(*dP, dQ, PP, Q, R, &dR);

            *dP = dR;
            for (i = 0; i <= dR; i++) PP[i] = R[i];
          }
        }
        a++;
      }
    }
    b++;
  }
  for (i = 0; i <= *dP; i++) P[i] = round2(PP[i].x);
}

// Hilbert root from Watson root as listed "An Improved Las Vegas Primality Test" by Kaltofen-Valente 1989, p.4
mpz_class Vegas(mpz_class u, mpz_class N, long d)
{
  mpz_class j;
  mpz_class two = 2, three = 3, four = 4, eight = 8;
  mpz_class A = 0;
  bool passthrough = false;


/*
  cout << "d = " << d << "\n";
  cout << "dmod3 = " << d % 3 << "\n";
  cout << "dmod8 = " << d % 8 << "\n";
  fflush(stdout);
*/

  switch (d % 8) {
    case -1:
      A = 4 * exp_mod(inverse(u, N), four, N);
      break;
    case -2:
    case -6:
      A = -4 * exp_mod(u, four, N);
      break;
    case -3:
      A = 16 * exp_mod(inverse(u, N), eight, N);
      break;
    case -4:
      A = -8 * exp_mod(u, two, N);
      break;
    case -5:
      A = 4 * exp_mod(inverse(u, N), two, N);
      break;
    case -7:
      A = exp_mod(inverse(u, N), eight, N);
      break;
    default:
      // D == 0 mod 32
      passthrough = true;
      break;
  }

  if (d % 3 != 0)
    A = exp_mod(A, three, N);

  j = inverse(A, N) * exp_mod(A - 16, three, N);
  j = modpos(j, N);

/*
  if (d % 8 == -3)
    j = u;
*/
  if (passthrough)
    j = u;

  return j;
}

/*
  Algorithm 1.5.3 (Modified Cornacchia). See "A Course
  in Computational Algebraic Number Theory" by Henri
  Cohen page 36. Let p be a prime number and D be a
  negative number such that D = 0 or 1 modulo 4 and
  | D | < 4 * p. This algorithm either outputs an
  integer solution (x, y) to the Diophantine equation
  x * x + | D | * y * y = 4 * p, or says that such a
  solution does not exist.
*/
int modified_Cornacchia(mpz_class D, mpz_class p,
                        mpz_class *x, mpz_class *y)
{
  int value = 0;
  mpz_class dd, xx;
  mpz_class a = 0, b = 0, c = 0, d = 0, e = 0;
  mpz_class l = 0, r = 0, x0 = 0;

/*
  if (p == 2) {
    a = D + 8;
    if (a > 0) {
      if (square_test(a, x)) {
        *y = 1;
        value = 1;
      }
    }
  }

  else {
 */
    if (JACOBI(D, p) != -1) {
      x0 = square_root_mod(D, p);

      dd = modpos(D, 2);
      xx = modpos(x0, 2);

      if (dd != xx)
        x0 = p - x0;

      a = p * 2;
      b = x0;
      c = sqrt(p);

      l = c * 2;

      while (b > l) {
        r = modpos(a, b);
        a = b;
        b = r;
      }

      c = p * 4;
      a = b * b;
      e = c - a;
      d = abs(D);

      c = e / d;
      r = modpos(e, d);

      if ((r == 0) && square_test(c, y)) {
        *x = b;
        value = 1;
      }
    }
//  }

  return value;
}

void zpoly_print(long da, mpz_class *za)
{
  long i;

  if (!Quiet) {
    for (i = da; i >= 0; i--) {
      cout << za[i];
      printf(" ");
    }
    printf("\n");
  }
}

/* multiplication (xy) modulo m */
void myzmulmod(mpz_t * res, mpz_t * x, mpz_t * y, mpz_t * m)
{
  mpz_t temp;

  mpz_init(temp);

  mpz_mul(temp, *x, *y);
  mpz_mod(*res, temp, *m);

  mpz_clear(temp);

}

/* division (x/y) modulo m */
void myzdivmod(mpz_t * res, mpz_t * x, mpz_t * y, mpz_t * m)
{
  int i;
  mpz_t h1;

  mpz_init(h1);

  i = mpz_invert(h1, *y, *m);

  if (i == 0) {
    printf("inverse is undefined!\n");
    exit(0);
  }

  else
    myzmulmod(res, x, &h1, m);

  mpz_clear(h1);
}

void zpoly_copy(long da, mpz_t * za, mpz_t * zb, long *db)
{

  long i;



  *db = da;

  for (i = 0; i <= da; i++)
    mpz_set(zb[i], za[i]);

}

void zpoly_mul(long m, long n, mpz_t * za, mpz_t * zb, mpz_t * zc, long *p)
{

  long i, j, k;

  mpz_t zd, zai, zbk, zsum, zterm;

  mpz_init(zd);
  mpz_init(zai);
  mpz_init(zbk);
  mpz_init(zsum);
  mpz_init(zterm);

  *p = m + n;

  for (k = 0; k <= *p; k++) {

    mpz_set_ui(zsum, 0);

    for (i = 0; i <= k; i++) {

      j = k - i;

      if (i > m)
        mpz_set_ui(zai, 0);

      else
        mpz_set(zai, za[i]);

      if (j > n)
        mpz_set_ui(zbk, 0);

      else
        mpz_set(zbk, zb[j]);

      mpz_mul(zterm, zai, zbk);

      mpz_set(zd, zsum);

      mpz_add(zsum, zterm, zd);

    }

    mpz_set(zc[k], zsum);

  }

  mpz_clear(zd);
  mpz_clear(zai);
  mpz_clear(zbk);
  mpz_clear(zsum);
  mpz_clear(zterm);
}



void zpoly_div(long m, long n, mpz_t * zu, mpz_t * zv, mpz_t * zq, mpz_t * zr,
         long *p, long *s)
{

  long j, jk, k, nk;

  mpz_t za, zb, zvn;

  mpz_init(za);
  mpz_init(zb);
  mpz_init(zvn);

  mpz_set(zvn, zv[n]);

  for (j = 0; j <= m; j++)
    mpz_set(zr[j], zu[j]);

  if (m < n) {
    *p = 0, *s = m;
    mpz_set_ui(zq[0], 0);
  } else {
    *p = m - n, *s = n - 1;
    for (k = *p; k >= 0; k--) {
      nk = n + k;
      mpz_pow_ui(za, zvn, k);
      mpz_mul(zq[k], zr[nk], za);

      for (j = nk - 1; j >= 0; j--) {
        jk = j - k;
        if (jk >= 0) {
          mpz_mul(za, zvn, zr[j]);
          mpz_mul(zb, zr[nk], zv[jk]);
          mpz_sub(zr[j], za, zb);
        } else {
          mpz_set(za, zr[j]);
          mpz_mul(zr[j], zvn, za);
        }
      }
    }

    while (*p > 0 && mpz_cmp_ui(zq[*p], 0l) == 0)
      *p = *p - 1;

    while (*s > 0 && mpz_cmp_ui(zr[*s], 0l) == 0)
      *s = *s - 1;
  }

  mpz_clear(za);
  mpz_clear(zb);
  mpz_clear(zvn);
}

void zpoly_mod(mpz_t * zp, mpz_t * za, long *da)
{
  long i;
  mpz_t zb;
  mpz_init(zb);

  for (i = 0; i <= *da; i++) {
    mpz_mod(zb, za[i], *zp);
    mpz_set(za[i], zb);
  }

  while (*da > 0 && mpz_cmp_ui(za[*da], 0l) == 0)
    *da = *da - 1;

  mpz_clear(zb);

}

void zpoly_pow(long degreeA, long degreem, mpz_t * zn, mpz_t * zp, mpz_t * zA,
         mpz_t * zm, mpz_t * zs, long *ds)
{

  long dP, dq, dx = degreeA, i;

  mpz_t za, zb, zP[POLY_SIZE], zq[POLY_SIZE], zx[POLY_SIZE], zy[POLY_SIZE];
  mpz_init(za);
  mpz_init(zb);

  for (i = 0; i < POLY_SIZE; i++) {
    mpz_init(zP[i]);
    mpz_init(zq[i]);
    mpz_init(zx[i]);
    mpz_init(zy[i]);
  }

  *ds = 0;

  mpz_set(za, *zn);

  mpz_set_ui(zs[0], 1l);

  for (i = 0; i <= dx; i++)
    mpz_set(zx[i], zA[i]);

  while (mpz_cmp_ui(za, 0l) > 0) {

    if (mpz_tdiv_ui(za, 2l) == 1) {
      // s = (s * x) % m;
      zpoly_mul(*ds, dx, zs, zx, zP, &dP);
      zpoly_div(dP, degreem, zP, zm, zq, zs, &dq, ds);
      zpoly_mod(zp, zs, ds);
    }

    mpz_set(zb, za);
    mpz_tdiv_q_2exp(za, zb, 1l);

    if (mpz_cmp_ui(za, 0l) > 0) {
      // x = (x * x) % m;
      for (i = 0; i <= dx; i++)
        mpz_set(zy[i], zx[i]);

      zpoly_mul(dx, dx, zx, zy, zP, &dP);
      zpoly_div(dP, degreem, zP, zm, zq, zx, &dq, &dx);
      zpoly_mod(zp, zx, &dx);
    }
  }

  mpz_clear(za);
  mpz_clear(zb);

  for (i = 0; i < POLY_SIZE; i++) {
    mpz_clear(zP[i]);
    mpz_clear(zq[i]);
    mpz_clear(zx[i]);
    mpz_clear(zy[i]);
  }
}

void zpoly_sub(long da, long db, mpz_t * za, mpz_t * zb, mpz_t * zc, long *dc)
{
  long i;
  mpz_t zz;
  mpz_init(zz);
  mpz_set_ui(zz, 0);
  if (da >= db) {
    for (i = 0; i <= db; i++)
      mpz_sub(zc[i], za[i], zb[i]);

    for (i = db + 1; i <= da; i++)
      mpz_set(zc[i], za[i]);

    *dc = da;
  } else {
    for (i = 0; i <= da; i++)
      mpz_sub(zc[i], za[i], zb[i]);

    for (i = da + 1; i <= db; i++)
      mpz_sub(zc[i], zz, zb[i]);

    *dc = db;
  }

  mpz_clear(zz);
}

void zpoly_gcd(long degreeA, long degreeB, mpz_t * zp, mpz_t * zA, mpz_t * zB,
         mpz_t * za, long *da)
{
  int nonzero = 0, zero;
  long db, dq, dr, i;

  mpz_t zc;
  mpz_t zb[POLY_SIZE], zq[POLY_SIZE], zr[POLY_SIZE];
  mpz_init(zc);

  for (i = 0; i < POLY_SIZE; i++) {
    mpz_init(zb[i]);
    mpz_init(zq[i]);
    mpz_init(zr[i]);
  }

  if (degreeA > degreeB) {
    *da = degreeA;
    db = degreeB;

    for (i = 0; i <= *da; i++)
      mpz_set(za[i], zA[i]);

    for (i = 0; i <= db; i++)
      mpz_set(zb[i], zB[i]);
  } else {
    *da = degreeB;
    db = degreeA;

    for (i = 0; i <= *da; i++)
      mpz_set(za[i], zB[i]);

    for (i = 0; i <= db; i++)
      mpz_set(zb[i], zA[i]);
  }

  for (i = 0; i <= db && !nonzero; i++)
    nonzero = mpz_cmp_ui(zb[i], 0l) != 0;

  while (nonzero) {
    zpoly_div(*da, db, za, zb, zq, zr, &dq, &dr);

    for (i = 0; i <= dr; i++) {
      mpz_set(zc, zr[i]);
      mpz_mod(zr[i], zc, *zp);
    }

    zero = 1;
    for (i = dr; i >= 0 && zero; i--) {
      zero = mpz_cmp_ui(zr[i], 0l) == 0;
      if (zero && dr > 0)
        dr--;
    }

    for (i = 0; i <= db; i++)
      mpz_set(za[i], zb[i]);

    *da = db;

    for (i = 0; i <= dr; i++)
      mpz_set(zb[i], zr[i]);

    db = dr;

    nonzero = 0;
    for (i = 0; i <= db && !nonzero; i++)
      nonzero = mpz_cmp_ui(zb[i], 0l) != 0;
  }

  mpz_clear(zc);

  for (i = 0; i < POLY_SIZE; i++) {
    mpz_clear(zb[i]);
    mpz_clear(zq[i]);
    mpz_clear(zr[i]);
  }
}

/* find the roots of the polynomial zA modulo zp */
void Recurse(long degreeA, mpz_t * zp, mpz_t * zA, mpz_t * zroot,
       long *rootSize)
{
  long dd, degreeB, dq, dr, du = 1, i, flag = 0;
  mpz_t zD, za, zb, zc, ze;
  mpz_t zn, x0;
  mpz_t zB[POLY_SIZE], zd[POLY_SIZE];
  mpz_t zq[POLY_SIZE], zr[POLY_SIZE];
  mpz_t zu[2];
  mpz_init(zD);
  mpz_init(za);
  mpz_init(zb);
  mpz_init(zc);
  mpz_init(ze);
  mpz_init(zn);
  mpz_init(x0);

  for (i = 0; i < POLY_SIZE; i++) {
    mpz_init(zB[i]);
    mpz_init(zd[i]);
    mpz_init(zq[i]);
    mpz_init(zr[i]);
  }

  mpz_init(zu[0]);
  mpz_init(zu[1]);
  mpz_set_ui(x0, (long) 102);
  if (mpz_cmp_ui(zA[degreeA], (long) 1) != 0) {
    for (i = 0; i < (degreeA + 1); i++) {
      myzdivmod(&zA[i], &zA[i], &zA[degreeA], zp);
    }
  }
  while (degreeA != 1) {
    i = 0;
    do {
      mpz_sub_ui(za, *zp, 1l);
      mpz_tdiv_q_2exp(zn, za, 1l);
      mpz_add_ui(x0, x0, (long) 1);
      mpz_set(zu[0], x0);
      mpz_set_ui(zu[1], 1l);
      zpoly_mod(zp, zu, &du);
      zpoly_pow(du, degreeA, &zn, zp, zu, zA, zd, &dd);
      zpoly_mod(zp, zd, &dd);
      mpz_sub_ui(zd[0], zd[0], 1l);
      zpoly_gcd(dd, degreeA, zp, zd, zA, zB, &degreeB);
      zpoly_mod(zp, zB, &degreeB);
      if ((mpz_cmp_ui(x0, (long) 200) == 1)
        || *rootSize > (degreeA - 1)) {
        flag = 1;
        goto L1;
      }
    } while ((degreeB == 0 || degreeB == degreeA));

    if (degreeB >= 1 && flag != 1) {
      Recurse(degreeB, zp, zB, zroot, rootSize);

      zpoly_div(degreeA, degreeB, zA, zB, zq, zr, &dq, &dr);
      zpoly_mod(zp, zq, &dq);
      zpoly_copy(dq, zq, zA, &degreeA);

      Recurse(degreeA, zp, zA, zroot, rootSize);
    }
  }

  if (degreeA == 1) {
    mpz_invert(za, zA[1], *zp);
    mpz_mul(zb, zA[0], za);
    mpz_neg(zb, zb);
    mpz_mod(zroot[*rootSize], zb, *zp);

    *rootSize = *rootSize + 1;
  }

L1:
  mpz_clear(zD);
  mpz_clear(za);
  mpz_clear(zb);
  mpz_clear(zc);
  mpz_clear(ze);
  mpz_clear(zn);
  mpz_clear(x0);

  for (i = 0; i < POLY_SIZE; i++) {
    mpz_clear(zB[i]);
    mpz_clear(zd[i]);
    mpz_clear(zq[i]);
    mpz_clear(zr[i]);
  }

  mpz_clear(zu[0]);
  mpz_clear(zu[1]);
}

/* find a single root of polynomial zA modulo zp */
/*
void myRecurse(long degreeA, mpz_t * zp, mpz_t * zA, mpz_t * zroot,
         long *rootSize)
{

  long dd, degreeB, dq, dr, du = 1, i, flag = 0;

  mpz_t zn, x0, za;

  mpz_t zB[POLY_SIZE], zd[POLY_SIZE];

  mpz_t zq[POLY_SIZE], zr[POLY_SIZE];

  mpz_t zu[2];


  mpz_init(zn);
  mpz_init(x0);
  mpz_init(za);

  for (i = 0; i < POLY_SIZE; i++) {
    mpz_init(zB[i]);
    mpz_init(zd[i]);
    mpz_init(zq[i]);
    mpz_init(zr[i]);
  }

  mpz_init(zu[0]);
  mpz_init(zu[1]);


  mpz_set_ui(x0, (long) 102);


  if (mpz_cmp_ui(zA[degreeA], (long) 1) != 0) {
    for (i = 0; i < (degreeA + 1); i++) {
      myzdivmod(&zA[i], &zA[i], &zA[degreeA], zp);

    }
  }

  while (degreeA != 1) {
    i = 0;
    do {

      mpz_sub_ui(za, *zp, 1l);
      mpz_tdiv_q_2exp(zn, za, 1l);

      mpz_add_ui(x0, x0, (long) 1);

      mpz_set(zu[0], x0);
      mpz_set_ui(zu[1], 1l);

      i++;

      zpoly_mod(zp, zu, &du);

      zpoly_pow(du, degreeA, &zn, zp, zu, zA, zd, &dd);
      zpoly_mod(zp, zd, &dd);

      mpz_sub_ui(zd[0], zd[0], 1l);

      zpoly_gcd(dd, degreeA, zp, zd, zA, zB, &degreeB);

      zpoly_mod(zp, zB, &degreeB);


      if (i > 50 || *rootSize == 1) {
        flag = 1;

        goto L1;
      }


    } while ((degreeB == 0 || degreeB == degreeA));


    if (degreeB >= 1 && flag != 1) {
      if (degreeB > degreeA / 2) {
        zpoly_div(degreeA, degreeB, zA, zB, zq, zr, &dq, &dr);
        zpoly_mod(zp, zq, &dq);
        zpoly_copy(dq, zq, zA, &degreeA);

        Recurse(degreeA, zp, zA, zroot, rootSize);
      }

      else {
        zpoly_copy(degreeB, zB, zA, &degreeA);

        Recurse(degreeA, zp, zA, zroot, rootSize);

      }
    }

  }


  if (degreeA == 1) {
    mpz_invert(za, zA[1], *zp);
    mpz_mul(x0, zA[0], za);
    mpz_neg(x0, x0);
    mpz_mod(zroot[*rootSize], x0, *zp);

    *rootSize = *rootSize + 1;

  }

L1:

  mpz_clear(zn);
  mpz_clear(x0);
  mpz_clear(za);

  for (i = 0; i < POLY_SIZE; i++) {
    mpz_clear(zB[i]);
    mpz_clear(zd[i]);
    mpz_clear(zq[i]);
    mpz_clear(zr[i]);
  }

  mpz_clear(zu[0]);
  mpz_clear(zu[1]);
}
*/

/*
 Algorithm 1.6.1 (Roots Mod p). See "A Course
 in Computational Algebraic Number Theory" by
 Henri Cohen page 37.
 */
void FindRootsModuloAPrime(long degreeP, mpz_class p, mpz_class *P, mpz_class *root, long *rootSize)
{
  long i, j;
  mpz_class t;
  mpz_t zp, zA[POLY_SIZE], zroot[POLY_SIZE];

  mpz_init(zp);

  for (i = 0; i < POLY_SIZE; i++) {
    mpz_init(zA[i]);
    mpz_init(zroot[i]);
  }

  mpz_set(zp, p.get_mpz_t());

  for (i = 0; i < POLY_SIZE; i++) {
    mpz_set(zA[i], P[i].get_mpz_t());
    mpz_set(zroot[i], root[i].get_mpz_t());
  }

  // added JGW 2011-03-13
  *rootSize = 0;

  Recurse(degreeP, &zp, zA, zroot, rootSize);

  for (i = 0; i < POLY_SIZE; i++) {
    mpz_class temp_class(zA[i]);
    P[i] = temp_class;
    mpz_class temp_class2(zroot[i]);
    root[i] = temp_class2;
  }

  /* sort the roots using selection sort */
  for (i = 0; i < *rootSize - 1; i++) {
    for (j = i + 1; j < *rootSize; j++) {
      if (root[i] > root[j]) {
        t = root[i];
        root[i] = root[j];
        root[j] = t;
      }
    }
  }

  mpz_clear(zp);
  for (i = 0; i < POLY_SIZE; i++) {
    mpz_clear(zA[i]);
    mpz_clear(zroot[i]);
  }
}

/*
  Algorithm 10.3.3 (Lenstra's ECM). See "A Course
  in Computational Algebraic Number Theory" by
  Henri Cohen page 488.
*/
int LenstrasECM(mpz_class *N, mpz_class *g, mpz_class Bmax)
{
  int found = 0;
  mpz_class B = 1000l;
  mpz_class l, q, q1, newq;
  struct point x, y;
  mpz_class a, d = 0;
  long curve = 0;

  do {
    do {
      a = modpos(rand2(), *N);
      x.x = 0;
      x.y = 1;

      curve++;
//    cout << "B=" << B << ", curve#" << curve << ", a=" << a << "          \r";
//    fflush(stdout);

      newq = 2;
      for (; newq < B && found != 1;) {
        q = newq;
        q1 = q;
        l = B / q;
        while (q1 <= l)
          q1 *= q;

        found = multiply(a, q1, *N, x, &y, &d);

        x.x = y.x;
        x.y = y.y;
        //  cout << "X=" << x.zx << "," << x.zy << "\n";

        newq = nextp(q);
      }

      *g = gcd(d, *N);

    } while (curve < 20 && (*g == *N || *g == 1));

    B = B * 5;
    curve = 0;

  } while (B < Bmax && (*g == *N || *g == 1));

//  cout << "\n";
  if (*g == *N || *g == 1)
    return 0;
  else {
    *N /= *g;
    return 1;
  }
}

// return true if t < *q < m and *q prime
bool check_for_factor(mpz_class *q, mpz_class m, mpz_class t)
{
  *q = m;

  if (Rabin_Miller(*q))
    return false;

// trial division
  while (*q % 2 == 0)
    *q /= 2;
  while (*q % 3 == 0)
    *q /= 3;

  if (*q < t)
    return false;
  if (*q < m && Rabin_Miller(*q))
    return true;

  mpz_class d = 1;
  mpz_class oldq;

  do {
    oldq = *q;
    LenstrasECM(q, &d, Bmax);
    if (*q < t)
      return false;
    if (*q < m && Rabin_Miller(*q))
      return true;
  } while (*q < oldq);

  return false;
}

bool find_curve(mpz_class *a, mpz_class *b, long D, mpz_class N, mpz_class *root, long *rootSize)
{
  mpz_t zp, zA[POLY_SIZE];
  long i;

  mpz_init(zp);

  for (i = 0; i < POLY_SIZE; i++) {
    mpz_init(zA[i]);
  }

  if (D == -3) {
    *a = 0;
    *b = -1;
    *rootSize = 1;
    return true;
  }
  else if (D == -4) {
    *a = -1;
    *b = 0;
    *rootSize = 1;
    return true;
  }

  mpz_class T[POLY_SIZE];
  long dT;

/*
  cout << "D = " << D << "\n";
  cout << "Dmod3 = " << D % 3 << "\n";
  cout << "Dmod8 = " << D % 8 << "\n";
  fflush(stdout);
*/
  if (weber && D % 32 != 0)
    Weber(D, T, &dT);
  else
    Hilbert(D, T, &dT);

  if (!Quiet) {
    if (weber && D % 32 != 0)
      cout << "D = " << D << ", dW = " << dT << ", W = ";
    else
      cout << "D = " << D << ", dT = " << dT << ", T = ";
    fflush(stdout);
    zpoly_print(dT, T);
    fflush(stdout);
  }

  if (!(weber && D % 32 != 0)) {
    mpz_class SmallestCoeff, Root3;
    SmallestCoeff = T[0];
    if (cube_test(SmallestCoeff, &Root3) == 0) {
      if (!Quiet) {
      cout << "HCP loss of precision\n";
      fflush(stdout);
      }
      return false;
    }
  }

  mpz_set(zp, N.get_mpz_t());

  for (i = 0; i < POLY_SIZE; i++) {
    mpz_set(zA[i], T[i].get_mpz_t());
  }

  zpoly_mod(&zp, zA, &dT);

  for (i = 0; i < POLY_SIZE; i++) {
    mpz_class temp_class(zA[i]);
    T[i] = temp_class;
  }

  FindRootsModuloAPrime(dT, N, T, root, rootSize);

/*
  if (rootSize > 0) {
    printf("the roots are as follows:\n");
    for (long i = 0; i < rootSize; i++) {
      cout << root[i];
      printf(" ");
    }
    printf("\n");
    fflush(stdout);
  }
*/

  mpz_clear(zp);
  for (i = 0; i < POLY_SIZE; i++) {
    mpz_clear(zA[i]);
  }

  if (*rootSize > 0)
    return true;
  else
    return false;
}

int Atkin(mpz_class N)
/* returns 0 if N is composite, 1 if prime */
{
  int value1, value2;
  long p;
  mpz_class Ni = N, a, b, d, i = 0, m, q, t, x, y;
  struct point P, P1, P2;
  long n = 0, k = 0;
  long D;
  bool found = false, found2 = false;
  mpz_class u, v;
  mpz_class g;
  long points_tried = 0;
  mpz_class root[POLY_SIZE];
  long rootSize = 0;
  mpz_class j;

  do {
//  cout << "Step 2\n";
    if (Ni <= SIEVE_LIMIT) {
      p = 2;
      while (p <= sqrt(Ni)) {
        if (Ni % p == 0) {
          cout << "1 factor = " << p << "\n";
        fflush(stdout);
          return 0;
        }
        p++;
      }
      return 1;
    }
    if (!Rabin_Miller(Ni)) return 0;

    while (Bmax <= 2000000000 && Dmax <= 312500) {
      if (!Quiet) {
        cout << "Bmax = " << Bmax << "\n";
        cout << "Dmax = " << Dmax << "\n";
        fflush (stdout);
      }

      n = 1;
      while ((D = -n++) > -Dmax) {
        if (D%4 != 0 && D%4 != -3 && D%4 != 1)
          continue;

        found = false;
        found2 = false;
        if (verbose) {
          cout << n << '\r';
          fflush(stdout);
        }

//      cout << "Step 3\n";
        if (JACOBI(D + Ni, Ni) != 1)
          continue;
        if (!modified_Cornacchia(D, Ni, &u, &v))
          continue;

//      cout << "Step 4\n";
        t = sqrt(sqrt(Ni)) + 1;
        t = t * t;
        if (check_for_factor(&q, m = Ni + 1 + u, t))
          found = true;
        else if (check_for_factor(&q, m = Ni + 1 - u, t))
          found = true;
        else if (D == -4) {
          if (check_for_factor(&q, m = Ni + 1 + 2*v, t))
            found = true;
          else if (check_for_factor(&q, m = Ni + 1 - 2*v, t))
            found = true;
        }
        else if (D == -3) {
          if (check_for_factor(&q, m = Ni + 1 + (u + 3*v)/2, t))
            found = true;
          else if (check_for_factor(&q, m = Ni + 1 - (u + 3*v)/2, t))
            found = true;
          else if (check_for_factor(&q, m = Ni + 1 + (u - 3*v)/2, t))
            found = true;
          else if (check_for_factor(&q, m = Ni + 1 - (u - 3*v)/2, t))
            found = true;
        }

        if (!found)
          continue;

//      cout << "Step 6\n";
        rootSize = 0;
        if (!find_curve(&a, &b, D, Ni, root, &rootSize))
          continue;

        for (long roots_tried = 0; roots_tried < rootSize; roots_tried++) {
          if (D == -3 || D == -4) {
            a = modpos(a, Ni);
            b = modpos(b, Ni);
          }
          else {
            j = modpos(root[roots_tried], Ni);

            if (weber && D % 32 != 0) {
              if (!Quiet) {
                cout << "u = " << j << "\n";
                fflush(stdout);
              }
              if (D % 4 == 0)
                j = Vegas(j, Ni, D / 4);
              else
                j = Vegas(j, Ni, D);
            }

            if (!Quiet) {
              cout << "j = " << j << "\n";
              fflush(stdout);
            }

            mpz_class c;
            c = modpos(j * inverse(j - 1728, Ni), Ni);
            a = modpos(-3 * c, Ni);
            b = modpos(2 * c, Ni);
          }

//        cout << "Step 7\n";
          do {
            do g = modpos(rand2(), Ni); while (g == 0);
            if (JACOBI(g, Ni) != -1)
              continue;
            if (D == -3)  // Cohen
//          if (Ni % 3 == 1) // Studholme
              if (exp_mod(g, (Ni - 1)/3, Ni) == 1)
                continue;
            break;
          } while (true);

//        cout << "Step 8\n";
          points_tried = 0;
          do {
            do {
              do {
                do x = modpos(rand2(), Ni); while (x == 0);
                y = modpos(((x * x) % Ni * x) % Ni + a * x + b, Ni);
              } while (JACOBI(y, Ni) == -1);
              y = square_root_mod(y, Ni);
            } while (y == 0);

//          cout << "Step 9\n";
//          cout << "Trying point [" << x << "," << y << "]\n";
//          fflush(stdout);
            P.x = x, P.y = y;
            points_tried++;
            k = 0;

            do {
//            cout << "Step 12\n";
              value2 = multiply(a, m/q, Ni, P, &P2, &d);
              if (value2 == 1) {
                cout << "3 factor = " << d << "\n";
                fflush(stdout);
                return 0;
              }

              value1 = multiply(a, q, Ni, P2, &P1, &d);
              if (value1 == 1) {
                cout << "2 factor = " << d << "\n";
                fflush(stdout);
                return 0;
              }
              if ((value1 == -1) && (value2 == 0)) {
//              cout << "Step 13\n";
                found2 = true;
                break;
              }

//            cout << "Step 10\n";
              ++k;
              if (D == -3) {
                if (k >= 6)
                  break;
                b *= g;
              }
              else if (D == -4) {
                if (k >= 4)
                  break;
                a *= g;
              }
              else {
                if (k >= 2)
                  break;
                a *= (g * g);
                b *= (g * g * g);
              }

              a = modpos(a, Ni);
              b = modpos(b, Ni);
            } while (true);
          } while (!found2 && points_tried < 100);
          if (found2)
            break;
        }

        if (found2)
          break;
      }

      if (found2) {
        if (!staticBmax)
          Bmax = BMAX;
        if (!staticDmax)
          Dmax = DMAX;
        break;
      }
      else {
        if (!staticBmax)
          Bmax *= 10;
        if (!staticDmax)
          Dmax *= 5;
      }
    }

    if (Dmax > 312500) {
      if (!Quiet) {
        cout << "ProvePrime: ran out of discriminants\n";
        fflush(stdout);
      }
      return 0;
    }
    if (Bmax > 2000000000) {
      if (!Quiet) {
        cout << "ProvePrime: exceeded maximum factoring bounds\n";
        fflush(stdout);
      }
      return 0;
    }

    cout << "N[" << i << "] = " << Ni << "\n";
    cout << "a = " << a << "\n";
    cout << "b = " << b << "\n";
    cout << "m = " << m << "\n";
    cout << "q = " << q << "\n";
    cout << "P = (" << P.x << ", " << P.y << ")\n";
    cout << "P1 = (" << P1.x << ", " << P1.y << ")\n";
    cout << "P2 = (" << P2.x << ", " << P2.y << ")\n";

    fflush(stdout);

    i++;
    Ni = q;

  } while (!one);

  return 0;
}

/* The name of this program. */
const char* program_name;

void print_usage (FILE* stream, int exit_code)
{
  fprintf(stream, "Usage: %s options [ < inputfile ]\n", program_name);
  fprintf(stream,
      " -h  --help    Display this usage information.\n"
      " -q  --quiet   Print reduced messages.\n"
      " -Q  --Quiet   Print very reduced messages.\n"
      " -v  --verbose Print verbose messages (including progress marker).\n"
      " -o  --one   Run just one iteration before termination.\n"
      " -w  --weber   Use Weber polynomials (rather than Hilbert).\n"
      " -s  --seed    Set (pseudo-)random seed.\n"
      " -p  --precision Set GMP FP precision (default 10000).\n"
      " -B  --Bmax    Set max ECM factoring bound (default 2000 autoincrement).\n"
      " -D  --Dmax    Set max discriminant bound (default 20 autoincrement).\n"
      " -e  --error_shift Set precision in HCP generation (default 1000).\n");
  exit(exit_code);
}

int main(int argc, char* argv[])
{
  int next_option;

  const char* const short_options = "hqQvows:p:B:D:e:";
  const struct option long_options[] = {
    {"help", 0, NULL, 'h'},
    {"quiet", 0, NULL, 'q'},
    {"Quiet", 0, NULL, 'Q'},
    {"verbose", 0, NULL, 'v'},
    {"one", 0, NULL, 'o'},
    {"weber", 0, NULL, 'w'},
    {"seed", 1, NULL, 's'},
    {"precision", 1, NULL, 'p'},
    {"Bmax", 1, NULL, 'B'},
    {"Dmax", 1, NULL, 'D'},
    {"error_shift", 1, NULL, 'e'},
    {NULL,0,NULL,0}
  };

  program_name = argv[0];

  do {
    next_option = getopt_long(argc, argv, short_options, long_options, NULL);
    switch (next_option)
    {
      case 'h':
        print_usage(stdout, 0);

      case 'q':
        quiet = 1;
        break;

      case 'Q':
        Quiet = 1;
        break;

      case 'v':
        verbose = 1;
        break;

      case 'o':
        one = 1;
        break;

      case 'w':
        weber = 1;
        break;

      case 's':
        seed = atol(optarg);
        break;

      case 'p':
        precision = atol(optarg);
        break;

      case 'B':
        Bmax = atol(optarg);
        staticBmax = true;
        break;

      case 'D':
        Dmax = atol(optarg);
        staticDmax = true;
        break;

      case 'e':
        error_shift = atol(optarg);
        break;

      case '?':
        print_usage(stderr, 1);

      case -1:
        break;

      default:
        abort();
    }
  } while (next_option != -1);

  gmp_randinit_default (rstate);

  {
#if HAVE_GETTIMEOFDAY
    struct timeval tv;
    gettimeofday (&tv, NULL);
    if (!seed)
        seed = tv.tv_sec + tv.tv_usec;

#else
    time_t t;
    time (&t);
    if (!seed)
      seed = t;
#endif
  }

  gmp_randseed_ui (rstate, seed);

  if (!staticBmax)
    Bmax = BMAX;
  if (!staticDmax)
    Dmax = DMAX;
  if (!error_shift)
    error_shift = ERROR_SHIFT;
  if (!precision)
    precision = PRECISION;

  if (!Quiet) {
    cout << "random seed = " << seed << "\n";
    cout << "error_shift = " << error_shift << "\n";
    cout << "precision = " << precision << "\n";
    fflush(stdout);
  }

  mpz_class N;

  mpf_set_default_prec(precision);

  PIf.set_prec(precision);
  Ef.set_prec(precision);
  NATLOGONEPOINTNINEf.set_prec(precision);

  ONEPOINTNINE.set_prec(precision);

  ZERO.set_prec(precision);
  ONE.set_prec(precision);
  TWO.set_prec(precision);
  THREE.set_prec(precision);
  FOUR.set_prec(precision);
  FIVE.set_prec(precision);
  EIGHT.set_prec(precision);
  TWELVE.set_prec(precision);
  EIGHTEEN.set_prec(precision);
  FIFTYSEVEN.set_prec(precision);
  TWOHUNDREDTHIRTYNINE.set_prec(precision);

  PIf = FOUR * (atanf2(ONE/TWO) + atanf2(ONE/THREE)); // Hutton
/*
  PIf = FOUR *
   (TWELVE * atanf2(ONE/EIGHTEEN)
    + EIGHT * atanf2(ONE/FIFTYSEVEN)
     - FIVE * atanf2(ONE/TWOHUNDREDTHIRTYNINE));  // Gauss
*/

  Ef = expf2(ONE);
  NATLOGONEPOINTNINEf = logf2(ONEPOINTNINE);

  if (!Quiet && !quiet) {
    cout << "PI = ";
    print_mpf(PIf);
    cout << "\n";
    cout << "********************\n";
    cout << "E = ";
    print_mpf(Ef);
    cout << "\n";
    cout << "********************\n";
    cout << "NATLOGONEPOINTNINE = ";
    print_mpf(NATLOGONEPOINTNINEf);
    cout << "\n";
    cout << "********************\n";
    fflush(stdout);
  }

  if (verbose) {
    long D;
    long dCP;
    mpz_class CP[POLY_SIZE];

    if (weber)
      cout << "\nFirst few Weber Class Polynomials :\n";
    else
      cout << "\nFirst few Hilbert Class Polynomials :\n";

    for (int n = 1; n <= 20; n++) {
      D = -n;
      if (D%4 != 0 && D%4 != -3 && D%4 != 1)
        continue;

      if (weber && D % 32 != 0)
        Weber(D, CP, &dCP);
      else
        Hilbert(D, CP, &dCP);
      if (weber && D % 32 != 0)
        cout << "D = " << D << ", dW = " << dCP << ", W = ";
      else
        cout << "D = " << D << ", dT = " << dCP << ", T = ";
      fflush(stdout);
      zpoly_print(dCP, CP);
      fflush(stdout);
    }
    cout << "\n";
    fflush(stdout);
  }

  if (!Quiet && !quiet) {
    cout << "number to be tested:\n";
    fflush(stdout);
  }
  cin >> N;
  if (Atkin(N))
    cout << ("proven prime\n");
  else {
    if (!one) {
      cout << ("composite\n");
    }
  }
  fflush(stdout);

  return 0;
}
