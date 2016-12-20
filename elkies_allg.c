/****************************************************************

file elkies_allg.c

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


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<limits.h>
#include<gmp.h>
#include<math.h>
#include"festkomma.h"

#define        SUCHWEITE         1.0e14
#define        UNTERSCHRANKE     1.0e11

/* Tripel mit |x**3 + y**3 - z**3| < AUSGABE_SCHRANKE werden ausgegeben. */
#define        AUSGABE_SCHRANKE  1000


/* Die globalen Variablen halbe_schrittweite, halbe_fliesenbreite und
   fliesenversatz werden in berechne_fliesenformat errechnet.
   Diese Funktion wird alle FLIESE_NEU Durchlaeufe der aeuszeren Schleife
   aufgerufen.
   Der FAKTOR regelt dabei die Groesze der Fliese.

   Beachte, dabei:
   halbe_schrittweite ist die halbe Ausdehnung einer Fliese in x-Richtung.
   halbe_fliesenbreite ist die halbe Ausdehnung einer Fliese in y-Richtung.
   Es muss gelten:
      halbe_fliesenbreite = 1.001 * y''(x) * halbe_schrittweite**2 / 4
   mit
      y(x) := (1 - x**3)**(1/3).
   Der Inhalt der Fliese ist damit 1.001 * y''(x) * halbe_schrittweite**3.

   Dieser Wert wird als 1.001 * (FAKTOR/SUCHWEITE)**3 angesetzt.
   Aus FAKTOR, SUCHWEITE und der aktuellen zweiten Ableitung sind somit
   halbe_schrittweite und halbe_fliesenbreite berechenbar. */
#define        FLIESE_NEU        1000000
#define        FAKTOR            5.5
double  halbe_schrittweite, halbe_fliesenbreite;
mpx_t                   x_0, fliesenversatz, Ax;
mpx_t                 y_diff, y_inv, y_inv_diff;
char                      ausg[1000], file[100];
long                                      zeile;


/* Zaehler fuer die Gesamtzahl der bearbeiteten Fliesen. */
long                                    fliesen;

/* Diese mpf-Variablen sollten an sich lokale Variablen, hauptsaechlich in den
   Funktionen berechne_y_wert und berechne_drei_linearf bzw.
   berechne_drei_linearf_mpf, sein.
   Mit globalen Variablen laeuft der Code aber schneller. */
mpf_t                                tmp1, tmp2;


/* Nur ein Prototyp. */
void init (long v[3][3], mpx_t y_0_type, mpx_t x_0_type,
  double halbe_schrittweite, mpx_t fliesenversatz,
  double halbe_fliesenbreite, double suchweite);


void out () {
 FILE *fp;

 fp = fopen (file, "a");
 fprintf (fp, "%6ld \t ", zeile);
 zeile++;
 fprintf (fp, ausg);
 fclose (fp);
}


/* pr := m1 * m2. Produkt zweier 3x3-Matrizen mit ganzzahligen Eintraegen.
   Wir rechnen in long und achten auf overflow ueber \pm2**63.
   pr = m1 und pr = m2 sind erlaubt. */
inline long matrix_prod (long pr[3][3],
                         long m1[3][3], long m2[3][3]) {
 long               i, j;
 long         prod[3][3];
 long  prodh, argh, argl;
 long                err;

 err = 0;
 for (i = 0; i < 3; i++) {
  for (j = 0; j < 3; j++) {
   smul_ppmm (prodh, prod[i][j], m1[i][0], m2[0][j]);
   smul_ppmm (argh, argl, m1[i][1], m2[1][j]);
   addf_ssaaaa (prodh, prod[i][j], prodh, prod[i][j], argh, argl);
   smul_ppmm (argh, argl, m1[i][2], m2[2][j]);
   addf_ssaaaa (prodh, prod[i][j], prodh, prod[i][j], argh, argl);
   if (((prod[i][j] >= 0) && (prodh != 0))
     || ((prod[i][j] < 0) && (prodh != -1))) {
    /* sprintf (ausg, "Overflow beim Matrizenprodukt. "); out ();
       sprintf (ausg, " prod[i][j] = %ld, prodh = %ld.\n",
                                              prod[i][j], prodh); out ();
    sprintf (ausg, "m1 = [%ld %ld %ld], [%ld %ld %ld], [%ld %ld %ld]\n",
         m1[0][0], m1[0][1], m1[0][2],
         m1[1][0], m1[1][1], m1[1][2],
         m1[2][0], m1[2][1], m1[2][2]); out ();
    sprintf (ausg, "m2 = [%ld %ld %ld], [%ld %ld %ld], [%ld %ld %ld]\n",
         m2[0][0], m2[0][1], m2[0][2],
         m2[1][0], m2[1][1], m2[1][2],
         m2[2][0], m2[2][1], m2[2][2]); out (); */
    err = 1; /* Overflow! */
   }
  }
 }
 /* Umkopieren. Noetig fuer den Fall, dass pr an die Adresse von m1 oder m2
    geschrieben werden soll. */
 for (i = 0; i < 3; i++) {
  for (j = 0; j < 3; j++) {
   pr[i][j] = prod[i][j];
  }
 }
 return (err);
}


/* Berechnung von halbe_schrittweite und halbe_fliesenbreite an der aktuellen
   Stelle x_0 \in [0,1]. */
void berechne_fliesenformat (mpx_t schrittweite) {
 double  x, zweite_abl, zusatzsummand;
 long                              fl;

 /* Output nur alle 100*FLIESE_NEU Fliesen. */
 fl = fliesen % (100 * FLIESE_NEU);

 x = mpx_get_d (x_0);

 if (fl == 0) {
  sprintf (ausg, "\n"); out ();
  sprintf (ausg, "x_0 = %.15f.\n", x); out ();
 }

 zweite_abl = 2*x / pow (1 - x*x*x, 5.0/3.0);
 halbe_schrittweite = pow (zweite_abl, -1.0/3.0) * FAKTOR / SUCHWEITE;

 zusatzsummand = AUSGABE_SCHRANKE
                              / (UNTERSCHRANKE * UNTERSCHRANKE * UNTERSCHRANKE);
 halbe_fliesenbreite =  1.001 * zweite_abl *
                                halbe_schrittweite*halbe_schrittweite / 4
                      + zusatzsummand;
                     /* Zusatzsummand soll sicherstellen, dass Loesungen mit
                        Hoehe >UNTERSCHRANKE garantiert gefunden werden. */
 if (fl == 0) {
  sprintf (ausg, "halbe_schrittweite = %.9e, halbe_fliesenbreite = %.9e.\n",
                               halbe_schrittweite, halbe_fliesenbreite);
  out ();
 }

 mpx_set_d (schrittweite, 2 * halbe_schrittweite);
 mpx_set_d (fliesenversatz,
                  zweite_abl * halbe_schrittweite*halbe_schrittweite / 4);
}


/***************************************************************************
 *
 * Eigentlicher Code
 *
 **************************************************************************/


/* Wird einmal pro FLIESE_NEU Fliesen aufgerufen.
   Berechnet 1/y_0 (- 1) aufwendig und exakt. */
void y_inv_init (mpx_t y_0) {
 /* Initialisierung von y_inv. */
 mpf_set_ui (tmp1, 1);
 mpf_set_mpx (tmp2, y_0);
 mpf_div (tmp1, tmp1, tmp2);
 mpf_sub_ui (tmp1, tmp1, 1);
 mpx_set_mpf (y_inv, tmp1);
 /* gmp_sprintf (ausg, "y_inv = %.*Ff.\n", 40, tmp1); */
}


/* Wird einmal pro FLIESE_NEU Fliesen aufgerufen.
   Initialisiert y_diff als y'(x_0) * schrittweite.
   Berechnet y_inv_diff als y_diff / y_0**2. */
void y_diff_init (mpx_t y_0, mpx_t schrittweite) {
 double  diff, y;

 /* Initialisierung von y_diff. */
 mpx_mul (y_diff, Ax, schrittweite);
 /* Im Falle y_0 < x_0 ist in Wirklichkeit Ax > 1. Korrigieren hier nach. */
 if ((y_0[1] < x_0[1]) || ((y_0[1] == x_0[1]) && (y_0[0] <= x_0[0])))
  mpx_add (y_diff, y_diff, schrittweite);
 /* mpf_set_mpx (tmp1, y_diff);
 gmp_sprintf (ausg, "y_diff = %.*Ff.\n\n", 40, tmp1); */

 /* Initialisierung von y_inv_diff. */
 y = mpx_get_d (y_0);
 diff = mpx_get_d (y_diff) / (y*y);
 mpx_set_d (y_inv_diff, diff);
 /* mpf_set_mpx (tmp1, y_inv_diff);
 gmp_sprintf (ausg, "y_inv_diff = %.*Ff.\n", 40, tmp1); */
}


/* y_0 := (1 - x_0**3)**(1/3). */
inline void berechne_y_wert (mpx_t y_0) {
 double        nenner, diff;
 mpx_t  tmpx1, tmpx2, tmpx3;

 mpx_add (y_0, y_0, y_diff);
 /* y_0 -= schrittweite * A. Lineare Verbesserung des Startwerts.
    Minus, weil wir rueckwaerts durch das Intervall laufen.
    y_diff wird einmal pro 1000000 Fliesen ausgerechnet.
    A ist immer negativ, also (-schrittweite)*A > 0.
    Beachte, y_diff > 0 nach Initialisierung. */

 mpx_mul (tmpx1, x_0, x_0);
 mpx_mul (tmpx2, tmpx1, x_0);
 tmpx3[0] = ~tmpx2[0]; tmpx3[1] = ~tmpx2[1];
 /* y3 = 1 - x_0*x_0*x_0; */
 /* Wir machen einen Fehler von exakt 2**(-128). */

 mpx_mul (tmpx1, y_0, y_0);
 nenner = 3.0 * mpx_get_d (tmpx1);

 mpx_mul (tmpx2, y_0, tmpx1);
 mpx_sub (tmpx1, tmpx2, tmpx3);
 diff = mpx_get_d (tmpx1);
 if (diff > 0.5) {
  sprintf (ausg, "Underflow.\n"); out ();
  exit (0);
 }
 diff /= nenner;

 mpx_set_d (tmpx1, diff);
 mpx_sub (y_0, y_0, tmpx1);
 /* Eine Newton-Iteration. */
 /* y_0 = (2*y_0 + y3/(y_0*y_0)) / 3
        = (2*y_0*y_0*y_0 + y3) / (3*y_0*y_0)
        = y_0 + (y3 - y_0*y_0*y_0) / (3*y_0*y_0);
    Differenz kann in double berechnet werden. */
 /* Eine Newton-Iteration. */

 /* mpf_set_mpx (tmp1, y_0);
 gmp_sprintf (ausg, "y_0 = %.*Ff.\n", 40, tmp1); out (); */
}


/* Schreibt y'(x_0) in die globale Variable Ax. */
inline void berechne_y_strich (mpx_t y_0) {
 mpx_t  tmpx1, tmpx2;

 /* Berechne zunaechst 1/y_0. */

 /* Lineare Verbesserung des Startwerts. */
 mpx_sub (y_inv, y_inv, y_inv_diff);

 /* Jetzt Newton-Iteration inv_neu = inv + (1 - y_0 * inv) * inv.
 Da y_inv den Wert 1/y - 1 annaehert, muessen wir
      y_inv = y_inv + (1 - y_0 * (1 + y_inv)) * (1 + y_inv)
            = y_inv + (1 - y_0 - y_0 * y_inv) + (1 - y_0 - y_0 * y_inv) * y_inv
 rechnen. */
 mpx_mul (tmpx1, y_0, y_inv);
 mpx_add (tmpx2, tmpx1, y_0); /* tmpx2 ist y_0 + y_0 * y_inv. */
 tmpx2[0] = ~tmpx2[0]; tmpx2[1] = ~tmpx2[1];
 /* tmpx2 ist jetzt 1 - y_0 - y_0 * y_inv. Fehler von 2**(-128) ignoriert. */
 mpx_mul (tmpx1, tmpx2, y_inv);
 mpx_add (tmpx1, tmpx1, tmpx2);
 /* tmpx1 ist (1 - y_0 - y_0 * y_inv) + (1 - y_0 - y_0 * y_inv) * y_inv. */
 mpx_add (y_inv, y_inv, tmpx1);
 /* mpf_set_mpx (tmp1, y_inv);
 gmp_sprintf (ausg, "y_inv = %.*Ff.\n", 40, tmp1); */

 /* Berechne jetzt 
         y'(x_0) = Ax := x_0**2 / y_0**2
                 = x_0**2 * (1 + y_inv)**2
                 = x_0**2 + x_0**2 * (y_inv**2 + 2*y_inv). */
 mpx_mul (tmpx1, y_inv, y_inv);
 mpx_add (tmpx1, tmpx1, y_inv);
 mpx_add (tmpx1, tmpx1, y_inv); /* tmpx1 ist jetzt y_inv**2 + 2*y_inv. */
 mpx_mul (tmpx2, x_0, x_0); /* tmpx2 ist x_0**2. */
 mpx_mul (Ax, tmpx2, tmpx1); /* Ax ist jetzt x_0**2 * (y_inv**2 + 2*y_inv). */
 mpx_add (Ax, Ax, tmpx2);
 /* mpf_set_mpx (tmp1, Ax);
 gmp_sprintf (ausg, "A = %.*Ff.\n", 40, tmp1); */

 /* Ax = - x_0*x_0 / (y_0*y_0); */
 /* Ableitung von
        y(x) := (1 - x**3)**(1/3)
    nach dem Satz ueber implizite Funktionen. */
 /* A = (y_1 - y_0) / (x_1 - x_0); */
 /* A = - x_0*x_0 / pow (1 - x_0*x_0*x_0, 2.0/3.0); */
}


inline void berechne_drei_linearf
  (double l[3][3], long v[3][3], mpx_t y_0, double d, double N) {
 long                         i;
 mpx_t  tmpx1, tmpx2, tmpx3, Bx;


 /* mpf_set_mpx (tmp1, x_0);
 gmp_sprintf (ausg, "\nx_0 = %.*Ff.\n", 40, tmp1); out (); */
 mpx_mul (tmpx1, Ax, x_0);
 /* Im Falle y_0 < x_0 ist in Wirklichkeit Ax > 1. Korrigieren hier nach. */
 if ((y_0[1] < x_0[1]) || ((y_0[1] == x_0[1]) && (y_0[0] <= x_0[0])))
  mpx_add (tmpx1, tmpx1, x_0);

 /* mpf_set_mpx (tmp1, fliesenversatz);
 gmp_sprintf (ausg, "fliesenversatz = %.*Ff.\n", 40, tmp1); out (); */
 mpx_sub (tmpx2, y_0, fliesenversatz);
 mpx_add (Bx, tmpx1, tmpx2); /* Ax ist in Wirklichkeit negativ. */
 /* Bx = y_0 - fliesenversatz - Ax * x_0; */
 /* Es gilt 1 <= Bx <= 1.6. Speichern in Wirklichkeit Bx - 1. */
 /* mpf_set_mpx (tmp1, Ax); mpf_set_mpx (tmp2, Bx);
 gmp_sprintf (ausg, "A = %.*Ff.\nB = %.*Ff.\n", 40, tmp1, 40, tmp2); out (); */

 for (i = 0; i < 3; i++) {
  l[0][i] = v[i][2] / N;
  /* sprintf (ausg, "l[0][%ld] = %f.\n", i, l[0][i]); out (); */
  /* Hier keine Probleme mit der Genauigkeit. */

  mpx_mul_si (tmpx1, x_0, -v[i][2]);
  mpxb_add_si (tmpx2, tmpx1, v[i][0]);
  /* tmpx = v[i][0] - v[i][2] * x_0; */
  l[1][i] = mpxg_get_d_simple (tmpx2);
  l[1][i] /= (halbe_schrittweite * N);
  /* sprintf (ausg, "l[1][%ld] = %f.\n", i, l[1][i]); out (); */
  /* Es ist schrittweite*N \approx 100.
     Wir brauchen also tmp2 auf >2 (7?) Nachkommastellen.
     Hierfuer sollte double ausreichen.
     tmp2 ist auf 128 Bit, also >=23 Nachkommastellen genau. */

  mpx_mul_si (tmpx1, Ax, v[i][0]);
  /* Im Falle y_0 < x_0 ist in Wirklichkeit Ax > 1. Korrigieren hier nach. */
  if ((y_0[1] < x_0[1]) || ((y_0[1] == x_0[1]) && (y_0[0] <= x_0[0])))
   mpxb_add_si (tmpx1, tmpx1, v[i][0]);

  mpx_mul_si (tmpx2, Bx, -v[i][2]);
  mpx_add (tmpx3, tmpx1, tmpx2); /* Ax ist in Wirklichkeit negativ. */
  mpxb_add_si (tmpx1, tmpx3, v[i][1] - v[i][2]); /* Bx \in [1..2]. */
  /* tmpx = v[i][1] - Ax * v[i][0] - Bx * v[i][2]; */
  l[2][i] = mpxg_get_d (tmpx1);
  l[2][i] /= (d*N);
  /* sprintf (ausg, "l[2][%ld] = %f.\n", i, l[2][i]); out (); */
  /* l[2][i] reagiert am sensibelsten auf zu wenig Precision.
     Es ist d*N \approx 1.0e-15.
     Wir brauchen also tmp2 auf > 15 (20?) Dezimalstellen hinter dem Komma.
     Also Ax und Bx auf 35 Dezimalstellen, entsprechend 128 Bit. */
 }
}


/* Makros fuer lll. */
#define scal_prod(prod, l, vec1, vec2)                            \
do {                                                              \
 long          t1, t2;                                            \
 double  l1[3], l2[3];                                            \
                                                                  \
 for (t1 = 0; t1 < 3; t1++)                                       \
  l1[t1] = l2[t1] = 0;                                            \
 for (t1 = 0; t1 < 3; t1++) {                                     \
  for (t2 = 0; t2 < 3; t2++) {                                    \
   l1[t1] += l[t1][t2] * vec1[t2];                                \
   l2[t1] += l[t1][t2] * vec2[t2];                                \
  }                                                               \
 }                                                                \
 prod = 0;                                                        \
 for (t1 = 0; t1 < 3; t1++)                                       \
  prod += l1[t1]*l2[t1];                                          \
} while (0);



#define gram()                                                    \
 do {                                                             \
  long  t1;                                                       \
                                                                  \
  /* sprintf (ausg, "Gram!\n"); out (); */                        \
  for (t1 = 0; t1 < 3; t1++)                                      \
   vec_gram[0][t1] = vec[0][t1];      /* double = long */         \
  scal_prod (B[0], l, vec_gram[0], vec_gram[0]);                  \
  scal_prod (mu[1][0], l, vec[1], vec_gram[0]); mu[1][0] /= B[0]; \
                                                                  \
  for (t1 = 0; t1 < 3; t1++)                                      \
   vec_gram[1][t1] = vec[1][t1] - mu[1][0] * vec_gram[0][t1];     \
  scal_prod (B[1], l, vec_gram[1], vec_gram[1]);                  \
  scal_prod (mu[2][0], l, vec[2], vec_gram[0]); mu[2][0] /= B[0]; \
  scal_prod (mu[2][1], l, vec[2], vec_gram[1]); mu[2][1] /= B[1]; \
                                                                  \
  for (t1 = 0; t1 < 3; t1++)                                      \
   vec_gram[2][t1] = vec[2][t1] - mu[2][0] * vec_gram[0][t1]      \
                                - mu[2][1] * vec_gram[1][t1];     \
  scal_prod (B[2], l, vec_gram[2], vec_gram[2]);                  \
} while (0)


#define red_k_k1()                                                \
 do {                                                             \
  /* sprintf (ausg, "Mache RED (%ld, %ld).\n", k, k-1); out (); */\
  tm = floor (mu[k][k-1] + 0.5);                                  \
  /* if (fabs (tm) > (1LU << 63) - 1) {                           \
   sprintf (ausg, "q ist zu lang.\n"); out ();                    \
   exit (0);                                                      \
  } */                                                            \
  q = lround (tm);                                                \
  if (q != 0) {                                                   \
   for (i = 0; i < 3; i++)                                        \
    vec[k][i] -= q * vec[k-1][i];                                 \
   mu[k][k-1] -= tm;                                              \
   if (k == 2)                                                    \
    mu[2][0] -= tm * mu[1][0];                                    \
  }                                                               \
} while (0)


#define red_2_0()                                                 \
 do {                                                             \
  /* sprintf (ausg, "Mache RED (2, 0).\n"); out (); */            \
  tm = floor (mu[2][0] + 0.5);                                    \
  q = lround (tm);                                                \
  /* if (fabs (tm) > (1LU << 63) - 1) {                           \
   sprintf (ausg, "q ist zu lang.\n"); out ();                    \
   exit (0);                                                      \
  } */                                                            \
  if (q != 0)                                                     \
   for (i = 0; i < 3; i++)                                        \
    vec[2][i] -= q * vec[0][i];                                   \
} while (0)


#define lll_swap()                                                \
do {                                                              \
 /* sprintf (ausg, "Mache SWAP (%ld).\n", k); out (); */          \
 for (i = 0; i < 3; i++) {                                        \
  tmp = vec[k][i];                                                \
  vec[k][i] = vec[k-1][i];                                        \
  vec[k-1][i] = tmp;                                              \
 }                                                                \
                                                                  \
 muc = mu[k][k-1];                                                \
 Bc = B[k] + muc*muc * B[k-1];                                    \
 mu[k][k-1] = muc * B[k-1] / Bc;                                  \
 B[k] = B[k-1]*B[k] / Bc;                                         \
 B[k-1] = Bc;                                                     \
                                                                  \
 if (k == 1) {                                                    \
  t = mu[2][1];                                                   \
  mu[2][1] = mu[2][0] - muc*t;                                    \
  mu[2][0] = t + mu[1][0]*mu[2][1];                               \
 }                                                                \
 else { /* k == 2*/                                               \
  tm = mu[2][0];                                                  \
  mu[2][0] = mu[1][0];                                            \
  mu[1][0] = tm;                                                  \
 }                                                                \
 k = 1;                                                           \
} while (0)


/* LLL -- double-Version.
   Die Funktion erwartet die drei Linearformen als in l gegeben.
   Sie gibt in v einen Satz von fuer l kurzen Vektoren zurueck.

   Der Ablauf des Algorithmus ist (bis auf einen Fehler in gram) von H. Cohen
   uebernommen.

   Wir wenden lll in der Weise an, dasz l schon bezueglich einer bestimmten
   Basis ausgerechnet ist. Die Rueckgabematrix vec ist dann eine
   Transformationsmatrix zu einer neuen Basis aus kurzen Vektoren.

   ACHTUNG: Diese Funktion arbeitet nur, wenn l aus nicht zu groszen Zahlen
   besteht. Wir laufen an der Kurve entlang und hoffen, dass die Vektoren, die
   bei letzten Fliese kurz waren, jetzt nicht sehr lang sind. */
inline void lll (double l[3][3], long vec[3][3]) {
 long                    i, k;
 double  B[3], vec_gram[3][3];
 double              mu[3][3];
 long                  q, tmp;
 double        muc, Bc, t, tm;

 /* Standardbasis. Der erste Kandidat fuer die Transformationsmatrix. */
 vec[0][0] = 1;  vec[0][1] = 0;  vec[0][2] = 0;
 vec[1][0] = 0;  vec[1][1] = 1;  vec[1][2] = 0;
 vec[2][0] = 0;  vec[2][1] = 0;  vec[2][2] = 1;

 /* LLL. */
 k = 1;
 gram ();

 while (k < 3) {
  /* sprintf (ausg, "%ld\n", k); out (); */
  red_k_k1 ();

  /* Test LLL condition */
  if (B[k] < ((0.99 - mu[k][k-1]*mu[k][k-1]) * B[k-1])) {
  /* if (B[k] < 0.5 * B[k-1]) { */
   /* sprintf (ausg, "Swap mit k = %ld \n", k); out (); */
   lll_swap ();
   /* out (); */
  }
  else {
   if (k == 2)
    red_2_0 ();
   k++;
  }
 }
}


#define lf_neu(lf, l, e)                                          \
do {                                                              \
 /* long  t1, t2;                                                 \
 for (t1 = 0; t1 < 3; t1++)                                       \
  for (t2 = 0; t2 < 3; t2++)                                      \
   lf[t1][t2]                                                     \
   = l[t1][0]*e[t2][0] + l[t1][1]*e[t2][1] + l[t1][2]*e[t2][2]; */\
                                                                  \
 lf[0][0] = l[0][0]*e[0][0] + l[0][1]*e[0][1] + l[0][2]*e[0][2];  \
 lf[0][1] = l[0][0]*e[1][0] + l[0][1]*e[1][1] + l[0][2]*e[1][2];  \
 lf[0][2] = l[0][0]*e[2][0] + l[0][1]*e[2][1] + l[0][2]*e[2][2];  \
                                                                  \
 lf[1][0] = l[1][0]*e[0][0] + l[1][1]*e[0][1] + l[1][2]*e[0][2];  \
 lf[1][1] = l[1][0]*e[1][0] + l[1][1]*e[1][1] + l[1][2]*e[1][2];  \
 lf[1][2] = l[1][0]*e[2][0] + l[1][1]*e[2][1] + l[1][2]*e[2][2];  \
                                                                  \
 lf[2][0] = l[2][0]*e[0][0] + l[2][1]*e[0][1] + l[2][2]*e[0][2];  \
 lf[2][1] = l[2][0]*e[1][0] + l[2][1]*e[1][1] + l[2][2]*e[1][2];  \
 lf[2][2] = l[2][0]*e[2][0] + l[2][1]*e[2][1] + l[2][2]*e[2][2];  \
} while (0);


/******************************************************************************
 *
 * Der Code fuer die Suche nach Gitterpunkten in der Pyramide.
 *
 *****************************************************************************/

/* res = m * v, 3x3-Matrix mal Spaltenvektor. */
#define MMUL(res, m, v)                                           \
do {                                                              \
 /* long  k;                                                      \
                                                                  \
 for (k = 0; k < 3; k++)                                          \
  res[k] = m[k][0]*v[0] + m[k][1]*v[1] + m[k][2]*v[2]; */         \
                                                                  \
 res[0] = m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2];             \
 res[1] = m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2];             \
 res[2] = m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2];             \
} while (0)


/* Berechne die vier Ecken der Pyramide. */
inline void pyr_ecken
 (double pt1[3], double pt2[3], double pt3[3], double pt4[3], double lf[3][3]) {
 long                             i, j;
 double  det_inv, adj[3][3], inv[3][3];

 /* Die Ecken sind gegeben durch (l1, l2, l3) (pti) = vi.
    Wichtig: v1, ..., v4 sind in der Reihenfolge des Umlaufs gelistet. */
 double      v1[3] = {1.0,  1.0,  1.0};
 double      v2[3] = {1.0, -1.0,  1.0};
 double      v3[3] = {1.0, -1.0, -1.0};
 /* double      v4[3] = {1.0,  1.0, -1.0}; */

 /* Inverse Matrix bis auf Skalierung. */
 adj[0][0] = lf[1][1]*lf[2][2] - lf[1][2]*lf[2][1];
 adj[0][1] = - (lf[0][1]*lf[2][2] - lf[0][2]*lf[2][1]);
 adj[0][2] = lf[0][1]*lf[1][2] - lf[0][2]*lf[1][1];

 adj[1][0] = - (lf[1][0]*lf[2][2] - lf[1][2]*lf[2][0]);
 adj[1][1] = lf[0][0]*lf[2][2] - lf[0][2]*lf[2][0];
 adj[1][2] = - (lf[0][0]*lf[1][2] - lf[0][2]*lf[1][0]);

 adj[2][0] = lf[1][0]*lf[2][1] - lf[1][1]*lf[2][0];
 adj[2][1] = - (lf[0][0]*lf[2][1] - lf[0][1]*lf[2][0]);
 adj[2][2] = lf[0][0]*lf[1][1] - lf[0][1]*lf[1][0];

 /* Determinante */
 det_inv = 1.0 / (lf[0][0]*adj[0][0] + lf[0][1]*adj[1][0] + lf[0][2]*adj[2][0]);
 /* sprintf (ausg, "det_inv = %.15f\n", det_inv); out (); */
 /* Inverse */
 for (i = 0; i < 3; i++)
  for (j = 0; j < 3; j++)
   inv[i][j] = adj[i][j] * det_inv;

 MMUL (pt1, inv, v1);
 MMUL (pt2, inv, v2);
 MMUL (pt3, inv, v3);
 /* MMUL (pt4, inv, v4); */
 pt4[0] = pt1[0] + pt3[0] - pt2[0];
 pt4[1] = pt1[1] + pt3[1] - pt2[1];
 pt4[2] = pt1[2] + pt3[2] - pt2[2];

 /* sprintf (ausg, "Ecken:\n"); out ();
 sprintf (ausg, "pt1 = (%f, %f, %f)\n", pt1[0], pt1[1], pt1[2]); out ();
 sprintf (ausg, "pt2 = (%f, %f, %f)\n", pt2[0], pt2[1], pt2[2]); out ();
 sprintf (ausg, "pt3 = (%f, %f, %f)\n", pt3[0], pt3[1], pt3[2]); out ();
 sprintf (ausg, "pt4 = (%f, %f, %f)\n", pt4[0], pt4[1], pt4[2]); out (); */
}


/* z laeuft in der auszersten Schleife.
   Brauchen also Maximum und Minimum der z-Koordinaten der 5 Ecken
   der Pyramide. */
inline void z_schranken (long *z_anf, long *z_end,
                         double *p1, double *p2, double *p3, double *p4) {
 double  tmp1, tmp2;

 tmp1 = MIN (MIN (MIN (p1[2], p2[2]), MIN (p3[2], p4[2])), 0.0);
 tmp2 = MAX (MAX (MAX (p1[2], p2[2]), MAX (p3[2], p4[2])), 0.0);

 *z_anf = lround (tmp1 + 0.5); /* Aufrunden. */
 *z_end = lround (tmp2 - 0.5); /* Abrunden. */
 /* Ganze Zahlen werden eventuell falsch gerundet.
 Dies ist aber kein Bug, weil wir ohnehin nur die Punkte im Innern der Pyramide
 brauchen. */

 /* sprintf (ausg, "z-Schranken: [%ld, %ld]\n", *z_anf, *z_end); out (); */
}


/* Eine Kante dargestellt als array der Laenge 6.
   y = kante[2]*z + kante[3] und
   x = kante[4]*z + kante[5]
   fuer kante[0] <= z <= kante[1]. */
inline void kante (double *kante, double *anf, double *end) {
 double  dz_inv;

 dz_inv = 1 / (anf[2] - end[2]);

 if (anf[2] < end[2]) {
  kante[0] = anf[2];
  kante[1] = end[2];
 }
 else {
  kante[0] = end[2];
  kante[1] = anf[2];
 }
 /* Ist dies eine Division durch 0, dann wird sowieso kante[0] = kante[1]
 (und dies ist hoffentlich keine ganze Zahl). */
 kante[2] = (anf[1] - end[1]) * dz_inv; /* dy/dz */
 kante[3] = anf[1] - kante[2]*anf[2];

 kante[4] = (anf[0] - end[0]) * dz_inv; /* dx/dz */
 kante[5] = anf[0] - kante[4]*anf[2];

 /* sprintf (ausg, "Kante: y = %f*z + %f, x = %f*z + %f fuer %f <= z <= %f\n",
               kante[2], kante[3], kante[4], kante[5], kante[0], kante[1]);
 out (); */
}


/* Berechne die acht Kanten der Pyramide. */
inline void pyr_kanten (double kant[8][6],
                 double *p1, double *p2, double *p3, double *p4) {
 double  s[3] = {0.0, 0.0, 0.0};

 kante (kant[0], s, p1);
 kante (kant[1], s, p2);
 kante (kant[2], s, p3);
 kante (kant[3], s, p4);
 kante (kant[4], p1, p2);
 kante (kant[5], p2, p3);
 kante (kant[6], p3, p4);
 kante (kant[7], p4, p1);
}


/* x und y laufen in den beiden inneren Schleifen.
   Wir schneiden also die Pyramide mit der z = ... - Ebene und haben ein
   Polygon. Wir brauchen das Maximum und das Minimum der x- bzw. y-Koordinaten
   der Ecken dieses Polygons.
   Ecken entstehen durch den Schnitt der Ebene mit den Kanten. */
inline void x_und_y_schranken
                  (long *x_anf, long *x_end, long *y_anf, long *y_end,
                   double kant[8][6], long z) {
 int                          i;
 double                 *kan, k;
 double  xanf, xend, yanf, yend;

 xanf = LONG_MAX;
 xend = LONG_MIN;
 yanf = LONG_MAX;
 yend = LONG_MIN;

 /* Schleife durch die acht Kanten. */
 for (i = 0; i < 8; i++) {
  kan = kant[i]; /* Aktuelle Kante. */
  /* Treffen wir diese Kante ueberhaupt? */
  if ((kan[0] <= z) && (z <= kan[1])) {
   k = kan[2]*z + kan[3];
   yanf = MIN (yanf, k);
   yend = MAX (yend, k);
   k = kan[4]*z + kan[5];
   xanf = MIN (xanf, k);
   xend = MAX (xend, k);
  }
 }
 *y_anf = lround (yanf + 0.5); /* Aufrunden. */
 *y_end = lround (yend - 0.5); /* Abrunden. */
 *x_anf = lround (xanf + 0.5); /* Aufrunden. */
 *x_end = lround (xend - 0.5); /* Abrunden. */

 /* sprintf (ausg, "x_y-Schranken fuer z = %ld: [%ld, %ld] x [%ld, %ld]\n",
                                         z, *x_anf, *x_end, *y_anf, *y_end);
 out (); */
}


/* Ausgabe der Loesung. */
void post_proc (long vec_out0, long vec_out1, long vec_out2, long f0) {
 /* double  quot, x0, x1; */

 /* quot = ((double) vec_out0) / ((double) vec_out2); */
 /* sprintf (ausg, "%.18f\n", quot); out (); */

 /* x0 = mpx_get_d (x_0);
 x0 -= halbe_schrittweite;
 x1 = x0 + 2 * halbe_schrittweite; */
 /* sprintf (ausg, "Intervall [%.18f,%.18f]\n\n", x0, x1); out (); */

 /* if ((x0 < quot) && (x1 > quot)) */ {
  sprintf (ausg, "(%ld, %ld, %ld) Loesung fuer k = %ld.\n",
                                          vec_out0, vec_out1, vec_out2, f0);
  out ();
  sprintf (ausg, "\n"); out ();
  /* sprintf (ausg, "Loesung gefunden.\n"); out (); */
  /* exit (0); */
 }
}


/* Was ist v2**3 - v0**3 - v1**3 bei den aktuellen Werten von z, y und x?
   Alle Rechnungen modulo 2**64.
   Ausgabe der Tripel (v0, v1, v2), bei denen v0**3 + v1**3 - v2**3 wirklich
   gleich \pm3 ist. */
inline ulong ein_funktionswert (long x, long v00, long v01, long v02,
                                long vec_h0, long vec_h1, long vec_h2) {
 long  vec_out0, vec_out1, vec_out2;
 ulong                           f0;

 /* sprintf (ausg, "x = %ld.\n", x); out (); */

 vec_out0 = x * v00 + vec_h0; /* Waere hier Rechnung in ulong logischer? */
 vec_out1 = x * v01 + vec_h1; /* Haben sowieso keinen Overflow. */
 vec_out2 = x * v02 + vec_h2;

 /* Rechnung modulo 2**64. Deswegen Cast nach unsigned. */
 f0 = ((ulong) vec_out2) * ((ulong) vec_out2) * ((ulong) vec_out2) -
      ((ulong) vec_out0) * ((ulong) vec_out0) * ((ulong) vec_out0) -
      ((ulong) vec_out1) * ((ulong) vec_out1) * ((ulong) vec_out1);

 /* Rechnen in kleine_vect nur den Fall z >= 0.
 Deshalb brauchen wir an dieser Stelle plus und minus.
 Ausgabe natuerlich als signed long ints. */
 if ((f0 < AUSGABE_SCHRANKE) || (-f0 < AUSGABE_SCHRANKE))
  if (f0 != 0)
   post_proc (vec_out0, vec_out1, vec_out2, f0); /* Ausgabe */

 return (f0);
}


#define XY_SCHLEIFE {                                                       \
  for (y = y_anf; y <= y_end; y++) {                                        \
   /* sprintf (ausg, "Rechne y = %ld.\n", y); out (); */                    \
   vec_h0 = y * v[1][0] + z * v[2][0];                                      \
   vec_h1 = y * v[1][1] + z * v[2][1];                                      \
   vec_h2 = y * v[1][2] + z * v[2][2];                                      \
                                                                            \
   x = x_anf;                                                               \
   /* Dieser Funktionsaufruf verursacht eine Ausgabe, falls f0 oder -f0     \
   unterhalb der AUSGABE_SCHRANKE liegt. */                                 \
   f0 = ein_funktionswert (x, v00, v01, v02, vec_h0, vec_h1, vec_h2);       \
   anz++;                                                                   \
                                                                            \
   x++;                                                                     \
   f1 = ein_funktionswert (x, v00, v01, v02, vec_h0, vec_h1, vec_h2);       \
   anz++;                                                                   \
                                                                            \
   x++;                                                                     \
   f2 = ein_funktionswert (x, v00, v01, v02, vec_h0, vec_h1, vec_h2);       \
   anz++;                                                                   \
                                                                            \
   x++;                                                                     \
                                                                            \
   /* Differenzenschema.*/                                                  \
   d1 = f1 - f0;    d = f2 - f1;                                            \
   dd =  d - d1;                                                            \
   f = f2;                                                                  \
   dd += ddd; /* dd einen Schritt im Vorlauf. */                            \
   for (; x <= x_end; x++) {                                                \
    d  += dd;                                                               \
    dd += ddd;                                                              \
    f  += d;                                                                \
    /* if (f != ein_funktionswert (x, v00, v01, v02, vec_h0, vec_h1, vec_h2)) {\
     sprintf (ausg, "%lu %lu\n", f,                                         \
                 ein_funktionswert (x, v00, v01, v02, vec_h0, vec_h1, vec_h2));\
     out ();
     sprintf (ausg, "Bug im Differenzenschema!\n"); out ();                 \
     exit (0);                                                              \
    } */                                                                    \
    anz++;                                                                  \
    if ((f < AUSGABE_SCHRANKE) || (-f < AUSGABE_SCHRANKE))                  \
     /* Dient nur der Ausgabe. Rechnen den Funktionswert ueberfluessigerweise
     nochmals aus. Ist schneller als eine extra Funktion aufzurufen. */     \
     ein_funktionswert (x, v00, v01, v02, vec_h0, vec_h1, vec_h2);          \
   }                                                                        \
  }                                                                         \
} while (0);


/* Suche nach ganzen Vektoren mit
        0 < l1 < 1, |l2| < l1 und |l3| < l1.
   Alle Rechnungen finden bezueglich der reduzierten Basis v statt.
   Sollte (vec[0], vec[1], vec[2]) aus durch 3 teilbaren Zahlen bestehen,
   rufen wir post_proc gar nicht erst auf.
   (Wegen vec[0]**3 + vec[1]**3 + vec[2]**3 = 0 (mod 3) sind entweder alle
   drei Komponenten durch 3 teilbar oder gar keine.) */
inline long kleine_vect (double lf[3][3], long v[3][3]) {
 long                       x, y, z;
 long                  x_anf, x_end;
 long                  y_anf, y_end;
 long                  z_anf, z_end;
 long        vec_h0, vec_h1, vec_h2;
 long                           anz;
 long                 v00, v01, v02;
 double  p1[3], p2[3], p3[3], p4[3];
 double                  kant[8][6];
 ulong                f, f0, f1, f2;
 ulong               d, dd, ddd, d1;

 /* sprintf (ausg, "lll-Basis:\n[%ld %ld %ld],\n[%ld %ld %ld],\n[%ld %ld %ld]\n",
           v[0][0], v[0][1], v[0][2],
           v[1][0], v[1][1], v[1][2],
           v[2][0], v[2][1], v[2][2]);
 out ();
 sprintf (ausg, "Linearformen auf dem ersten Vektor: %f %f %f\n",
                                               lf[0][0], lf[1][0], lf[2][0]);
 out (); */

 /* Berechne die Ecken der Pyramide. Die Spitze ist der Nullpunkt. */
 pyr_ecken (p1, p2, p3, p4, lf);
 /* Berechne die Schranken fuer z. */
 z_schranken (&z_anf, &z_end, p1, p2, p3, p4);

 /* Berechne die acht Kanten der Pyramide. */
 pyr_kanten (kant, p1, p2, p3, p4);

 /* Dreifachloop durch die Pyramide. */
 anz = 0;
 v00 = v[0][0]; v01 = v[0][1]; v02 = v[0][2];
 /* Differenzenschema, erster Teil.
 ddd ist unabhaengig von z, ueber die gesamte Fliese konstant. */
 ddd =  ((ulong) 6) *
       (((ulong) v02) * ((ulong) v02) * ((ulong) v02)
       -((ulong) v00) * ((ulong) v00) * ((ulong) v00)
       -((ulong) v01) * ((ulong) v01) * ((ulong) v01));

 for (z = z_anf; z <= z_end; z++) {
  /* Berechne die Schranken fuer x und y. */
  x_und_y_schranken (&x_anf, &x_end, &y_anf, &y_end, kant, z);
  XY_SCHLEIFE;
 }
 return (anz);
}


/* Init. Die aeuszere Schleife. */
void rechne_intervall (mpx_t x_0_anf, mpx_t x_0_ende) {
 mpx_t         schrittweite;
 mpx_t                  y_0;
 double  lf[3][3], ln[3][3];
 long      e[3][3], v[3][3];
 long     zaehler, anz, err;

 x_0[0] = x_0_ende[0]; x_0[1] = x_0_ende[1];
 /* Laufen rueckwaerts von x_0_ende nach x_0_anf. */

 berechne_fliesenformat (schrittweite);

 /* Erster Schleifendurchlauf. Hier rechnen wir LLL mit viel Precision.
    Auszerdem wird y_0 initialisiert. */
 init (v, y_0, x_0,
           halbe_schrittweite, fliesenversatz, halbe_fliesenbreite, SUCHWEITE);
 y_inv_init (y_0);
 y_diff[0] = 0; y_diff[1] = 0; y_inv_diff[0] = 0; y_inv_diff[1] = 0;
 /* y_diff = y_inv_diff = 0. */
 sprintf (ausg, "Init fertig.\n"); out ();

 /* Scheife. Hier reicht fuer fast alles die Genauigkeit von double. */
 zaehler = FLIESE_NEU; fliesen = -FLIESE_NEU; err = 0;
 /* while (x_0 >= x_0_anf) ... . */
 while ((x_0[1] > x_0_anf[1]) ||
       ((x_0[1] == x_0_anf[1]) && (x_0[0] >= x_0_anf[0]))) {
  /* mpf_set_mpx (tmp1, x_0);
  gmp_sprintf (ausg, "\nx_0 = %.*Ff.\n", 40, tmp1); out (); */

  berechne_y_wert (y_0);
  /* Beim ersten Mal ist y_diff = 0 => y_0 wird richtig ausgerechnet. */
  /* mpf_set_mpx (tmp1, y_0);
  gmp_sprintf (ausg, "y_0 = %.*Ff.\n\n", 40, tmp1); out (); */
  berechne_y_strich (y_0);
  /* Beim ersten Mal ist y_inv_diff = 0 => y'(x_0) wird richtig ausgerechnet. */
  /* mpf_set_mpx (tmp1, Ax);
  gmp_sprintf (ausg, "A = %.*Ff.\n\n", 40, tmp1); out (); */

  if (zaehler == FLIESE_NEU) {
   zaehler = 0; fliesen += FLIESE_NEU;
   berechne_fliesenformat (schrittweite);
   y_inv_init (y_0); y_diff_init (y_0, schrittweite);
  }

  /* Letztes v war nicht korrekt wegen Overflow. */
  if (err > 0) {
   init (v, y_0, x_0,
           halbe_schrittweite, fliesenversatz, halbe_fliesenbreite, SUCHWEITE);
   mpf_set_mpx (tmp1, x_0);
   gmp_sprintf (ausg, "Neustart bei x_0 = %.*Ff.\n", 25, tmp1); out ();
  }

  berechne_drei_linearf (lf, v, y_0, halbe_fliesenbreite, SUCHWEITE);

  lll (lf, e);
  err = matrix_prod (v, e, v); /* v ist auf jeden Fall richtig modulo 2**64. */
 
  lf_neu (ln, lf, e);
  anz = kleine_vect (ln, v);
  /* sprintf (ausg, "Finden %ld Gitterpunkte.\n", anz); out (); */

  /* Fliesen gehen von x_0 nach rechts und links. */
  mpx_sub (x_0, x_0, schrittweite); zaehler++;
 }

 sprintf (ausg, "Insgesamt %ld Fliesen behandelt.\n", fliesen + zaehler);
 out ();
}


/* Wir wollen mit etwas wie
      elkies_allg 0.4 0.000001
   starten.

   Dabei ist 0.4 der Anfangswert.
   0.000001 ist die Intervalllaenge, die der Prozessor schaffen soll. */
int main (int argc, char *argv[]) {
 mpx_t   diff, x_0_anf, x_0_ende;

 mpf_set_default_prec (128);
 mpf_init (tmp1); mpf_init (tmp2);

 mpf_set_str (tmp1, argv[1], 10); mpx_set_mpf (x_0_anf, tmp1);
 mpf_set_str (tmp2, argv[2], 10); mpx_set_mpf (diff, tmp2);

 zeile = 1;
 mpx_add (x_0_ende, x_0_anf, diff);

 /* Name der Ausgabedatei */
 mpf_set_mpx (tmp1, x_0_anf); mpf_set_mpx (tmp2, x_0_ende);
 gmp_sprintf (file, "liste_allg_%.*Ff_%.*Ff.txt", 8, tmp1, 8, tmp2);

 /* Erste Ausgabe */
 gmp_sprintf (ausg, "Starte Rechnung von %.*Ff bis %.*Ff.\n",
                                                      15, tmp1, 15, tmp2);
 out ();

 rechne_intervall (x_0_anf, x_0_ende);

 mpf_clear (tmp1); mpf_clear (tmp2);
 return 0;
}

