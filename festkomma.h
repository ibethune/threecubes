/****************************************************************

file festkomma.h

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


#define MAX(h, i) ((h) > (i) ? (h) : (i))
#define MIN(h, i) ((h) < (i) ? (h) : (i))

/* Der g++-Compiler vertraegt die normalen Makros add_ssaaaa und subf_ddmmss
   nicht, weil da ein seltsamer Cast beim Schreiben gemacht werden soll. */
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


/* Kann zwei Sorten von Festkommazahlen kodieren. 
   1. Eine Zahl zwischen 0 und 1 mit Genauigkeit 2*(-128) nach der Regel
          fixed[1] * 2**(-64) + fixed[0] * 2**(-128).
      Datentyp "mpx_t".
   2. Eine vorzeichenbehaftete Zahl zwischen (-2**63) und 2**63 mit Genauigkeit
      2*(-64) nach der Regel
          (signed long) fixed[1] + fixed[0] * 2**(-64).
      Datentyp "mpxg_t". */
typedef unsigned long int mpx_t [2];


/* Die naive Umrechnung von mpxg_t in double. */
#define mpxg_get_d_simple(fixed)                                  \
  (ldexp (((double) (fixed)[0]), -64)                             \
+ ((double) ((signed long) (fixed)[1])))


/* Umrechnung von mpxg_t in double.
   Wollen bei Werten nahe 0 besonders grosze Genauigkeit.
   Auch bei solchen knapp unter 0, die als
         (-1) + fixed[0] * 2**(-64) kodiert
   sind. */
inline double mpxg_get_d (mpx_t fixed) {
 /* fixed[0] >= 2**63 <==> gebrochener Teil >= 0.5. */
 if ((signed long) fixed[0] < 0) {
  /* Kuenstlich komplizierte Umrechnung in double.
     Der erste Summand ist [erg] - 1, was eventuell betragsmaeszig klein ist
     und sehr genau in double umgerechnet werden kann. */
  return (- ldexp (((double) -fixed[0]), -64)
                 + ((double) (1 + (signed long) fixed[1])));
 }
 /* Die naive Umrechnung in double.
 fixed[1] soll als signed verstanden werden. */
 return (   ldexp (((double) fixed[0]), -64) 
                 + ((double) ((signed long) fixed[1])));
}


/* Die naive Umrechnung von mpx_t in double. */
#define mpx_get_d(fixed)                                          \
   (ldexp (((double) (fixed)[1]), -64)                            \
 + ldexp (((double) (fixed)[0]), -128))


/* Voraussetzung: 0 < floa < 1. */
inline void mpx_set_d (mpx_t erg, double floa) {
 floa = ldexp (floa, 64);
 erg[1] = (ulong) floor (floa);
 floa -= erg[1];
 floa = ldexp (floa, 64);
 erg[0] = (ulong) floor (floa);
 /* Hier darf kein lround stehen. Gibt Probleme bei Zahlen >2**63! */
}


/* Voraussetzungen: 0 < floa < 1. Precision von fl >=128 Bit. 
   Achtung: Extrem ineffizient! */
inline void mpx_set_mpf (mpx_t erg, mpf_t fl) {
 /* Effizienter Code saehe so aus:
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


/* Zur Ausgabe mit gmp_printf gedacht. Hat sicher Optimierungspotential. */
inline void mpf_set_mpx (mpf_t erg, mpx_t fixed) {
 mpf_set_ui (erg, fixed[0]);
 mpf_div_2exp (erg, erg, 64);
 mpf_add_ui (erg, erg, fixed[1]);
 mpf_div_2exp (erg, erg, 64);
}


/* Produkt zweier mpx_t-Zahlen (128-Bit Festkommazahlen zwischen 0 und 1).
   Nicht fuer Datentyp mpxg_t. */
inline void mpx_mul (mpx_t erg, mpx_t fak1, mpx_t fak2) {
 ulong  argh, argl;

 /* Hauptteil des Produkts. */
 umul_ppmm (erg[1], erg[0], fak1[1], fak2[1]);
 /* Multiplikation ueber Kreuz. Ignorieren das dritte Limb. */
 umul_ppmm (argh, argl, fak1[1], fak2[0]);
 addf_ssaaaa (erg[1], erg[0], erg[1], erg[0], 0L, argh);
 /* Multiplikation ueber Kreuz andersrum. Ignorieren wieder das dritte Limb. */
 umul_ppmm (argh, argl, fak1[0], fak2[1]);
 addf_ssaaaa (erg[1], erg[0], erg[1], erg[0], 0L, argh);
 /* Lassen das Produkt der niederwertigen Limbs komplett weg. */
}


/* Vorzeichenloses Festkomma * long.
   fix ist vorzeichenloses Festkomma, also vom Datentyp mpx_t.
   Das Ergebnis erg ist vorzeichenbehaftetes Festkomma vom Datentyp mpxg_t. */
inline void mpx_mul_si (mpx_t erg, mpx_t fix, long ganz) {
 ulong  argh, argl;

 /* Behandle ganz als vorzeichenlos. */
 umul_ppmm (erg[1], erg[0], fix[1], ganz);
 umul_ppmm (argh, argl, fix[0], ganz);
 /* Ignorieren argl. */
 addf_ssaaaa (erg[1], erg[0], erg[1], erg[0], 0, argh);

 /* Korrektur falls ganz < 0.
    Wir haben bei ganz mit einem Fehler von exakt 2**64 gerechnet. */
 if (ganz < 0)
  subf_ddmmss (erg[1], erg[0], erg[1], erg[0], fix[1], fix[0]);
}


/* Klappt bei mpx_t wie bei mpxg_t. */
#define mpx_add(summe, summand1, summand2)                        \
 addf_ssaaaa (summe[1], summe[0], summand1[1], summand1[0],       \
                                  summand2[1], summand2[0])


/* Klappt bei mpx_t wie bei mpxg_t. */
#define mpx_sub(differenz, minuend, subtrahend)                   \
 subf_ddmmss (differenz[1], differenz[0], minuend[1], minuend[0], \
                                    subtrahend[1], subtrahend[0])


/* Klappt bei summand1 vom Datentyp mpx_t oder mpxg_t. summand2 ist eine long.
   Auf Overflow ist auch hier, wie immer bei add und sub, bei der Anwendung
   zu achten. */
#define mpxb_add_si(summe, summand1, summand2)                    \
 addf_ssaaaa (summe[1], summe[0], summand1[1], summand1[0],       \
                                                    summand2, 0)
