/****************************************************************

file fixedpoint.h

Elkies package for the three cubes problem,
Version 1.0

Copyright (C) 2007 A.-S. Elsenhans and J. Jahnel

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

Address of the authors

	A.-S. Elsenhans/ J. Jahnel
	Math. Institut der Universitaet
	Bunsenstrasze 3--5
	D-37073 Goettingen
        Germany

WWW     http://www.uni-math.gwdg.de/jahnel

*****************************************************************/

#include<gmp.h>

typedef unsigned long int UDItype;
typedef long int           DItype;
typedef UDItype ulong;

#define MAX(h, i) ((h) > (i) ? (h) : (i))
#define MIN(h, i) ((h) < (i) ? (h) : (i))

/* g++ confuses the usual macros add_ssaaaa and subf_ddmmss, because of the unusual
waythe casting is used (?) */
#define addf_ssaaaa(sh, sl, ah, al, bh, bl)                       \
  __asm__ ("addq %5,%1\n\tadcq %3,%0"                             \
           : "=r" (sh),           "=&r" (sl)                      \
           : "0"  ((UDItype)(ah)), "g" ((UDItype)(bh)),           \
             "%1" ((UDItype)(al)), "g" ((UDItype)(bl)))

#define subf_ddmmss(sh, sl, ah, al, bh, bl)                       \
  __asm__ ("subq %5,%1\n\tsbbq %3,%0"                             \
           : "=r" (sh), "=&r" (sl)                                \
           : "0" ((UDItype)(ah)), "g" ((UDItype)(bh)),            \
             "1" ((UDItype)(al)), "g" ((UDItype)(bl)))

#define umul_ppmm(w1, w0, u, v)                                   \
  __asm__ ("mulq %3"                                              \
	   : "=a" (w0), "=d" (w1)                                 \
	   : "%0" ((UDItype)(u)), "rm" ((UDItype)(v)))

#define smul_ppmm(w1, w0, u, v)                                   \
  __asm__ ("imulq %3"                                             \
	   : "=a" (w0), "=d" (w1)                                 \
	   : "%0" ((DItype)(u)), "rm" ((DItype)(v)))


/* We can encode two different varieties of fixed-point numbers:
   1. A number between 0 and 1 with precision 2^-128, with the format
          fixed[1] * 2**(-64) + fixed[0] * 2**(-128).
      Datatype "mpx_t".
   2. A signed number between -2^63 and 2^63 with precision 2^-64, with the format
          (signed long) fixed[1] + fixed[0] * 2**(-64).
      Datatype "mpxg_t". */
typedef unsigned long int mpx_t [2];


/* Naive conversion from mpxg_t to double. */
#define mpxg_get_d_simple(fixed)                                  \
  (ldexp (((double) (fixed)[0]), -64)                             \
+ ((double) ((signed long) (fixed)[1])))


/* Convert mpxg_t to double.
   If the value is close to 0, the accuracy is very good
   Even for those just below 0, encoded as
         (-1) + fixed[0] * 2**(-64) kodiert
   */
inline double mpxg_get_d (mpx_t fixed) {
 /* fixed[0] >= 2**63 <==> broken(?) part >= 0.5. */
 if ((signed long) fixed[0] < 0) {
  /* Artificiallly complicated converstion into double.
     The first summand is [erg] - 1, which may be small in magnitude
     and can be converted very precisely into double. */
  return (- ldexp (((double) -fixed[0]), -64)
                 + ((double) (1 + (signed long) fixed[1])));
 }
 /* Naive conversion into double.
 fixed[1] is taken to be signed. */
 return (   ldexp (((double) fixed[0]), -64) 
                 + ((double) ((signed long) fixed[1])));
}


/* Naive conversion from mpx_t to double. */
#define mpx_get_d(fixed)                                          \
   (ldexp (((double) (fixed)[1]), -64)                            \
 + ldexp (((double) (fixed)[0]), -128))


/* Prerequisite: 0 < floa < 1. */
inline void mpx_set_d (mpx_t erg, double floa) {
 floa = ldexp (floa, 64);
 erg[1] = (ulong) floor (floa);
 floa -= erg[1];
 floa = ldexp (floa, 64);
 erg[0] = (ulong) floor (floa);
 /* Cannot use lround due to problems with numbers > 2^63! */
}


/* Prerequisite: 0 < floa < 1. Precision of fl >=128 bit. 
   Warning: Extremely inefficient! */
inline void mpx_set_mpf (mpx_t erg, mpf_t fl) {
 /* Efficient code looks like:
 mp_ptr  ptr;
 long    pos;
 ptr = PTR (fl);
 pos = ABSIZ (fl) - 1;
 erg[1] = ptr[pos];
 erg[0] = ptr[pos - 1]; */

 mpf_t  tmp;

 mpf_init (tmp);
 mpf_mul_2exp (tmp, fl, 64);
 erg[1] = mpf_get_ui (tmp);
 mpf_sub_ui (tmp, tmp, erg[1]);
 mpf_mul_2exp (tmp, tmp, 64);
 erg[0] = mpf_get_ui (tmp);
 mpf_clear (tmp);
}


/* For output with gmp_printf(). Can certainly be optimized. */
inline void mpf_set_mpx (mpf_t erg, mpx_t fixed) {
 mpf_set_ui (erg, fixed[0]);
 mpf_div_2exp (erg, erg, 64);
 mpf_add_ui (erg, erg, fixed[1]);
 mpf_div_2exp (erg, erg, 64);
}


/* Produce of two mpx_t-numbers (128-Bit fixed point number between 0 and 1).
   Not for mpxg_t types. */
inline void mpx_mul (mpx_t erg, mpx_t fak1, mpx_t fak2) {
 ulong  argh, argl;

 /* Main part of the product */
 umul_ppmm (erg[1], erg[0], fak1[1], fak2[1]);
 /* Cross-multiply, ignoring the third limb. */
 umul_ppmm (argh, argl, fak1[1], fak2[0]);
 addf_ssaaaa (erg[1], erg[0], erg[1], erg[0], 0L, argh);
 /* Cross-multiply again, stil ignoring the third limb. */
 umul_ppmm (argh, argl, fak1[0], fak2[1]);
 addf_ssaaaa (erg[1], erg[0], erg[1], erg[0], 0L, argh);
 /* Leaves out the prpduct of the low-value limbs entirely. */
}


/* Unsigned fixed-point multiplied by long.
   fix is unsigned fixed-point i.e. mpx_t type.
   The result erg is a signed fixed-point of type mpxg_t. */
inline void mpx_mul_si (mpx_t erg, mpx_t fix, long ganz) {
 ulong  argh, argl;

 /* Treat completely as unsigned */
 umul_ppmm (erg[1], erg[0], fix[1], ganz);
 umul_ppmm (argh, argl, fix[0], ganz);
 /* Ignore argl. */
 addf_ssaaaa (erg[1], erg[0], erg[1], erg[0], 0, argh);

 /* Correction when < 0.
    We computed with an error of exactly 2^64 */
 if (ganz < 0)
  subf_ddmmss (erg[1], erg[0], erg[1], erg[0], fix[1], fix[0]);
}


/* Sum two mpx_t into an mpxg_t */
#define mpx_add(summe, summand1, summand2)                        \
 addf_ssaaaa (summe[1], summe[0], summand1[1], summand1[0],       \
                                  summand2[1], summand2[0])


/* Subtract two  mpx_t into an mpxg_t. */
#define mpx_sub(differenz, minuend, subtrahend)                   \
 subf_ddmmss (differenz[1], differenz[0], minuend[1], minuend[0], \
                                    subtrahend[1], subtrahend[0])


/* Sum of and mpx_t and a long.
   We pay attention to make sure overflow is handled correctly.*/
#define mpxb_add_si(summe, summand1, summand2)                    \
 addf_ssaaaa (summe[1], summe[0], summand1[1], summand1[0],       \
                                                    summand2, 0)
