/****************************************************************

file elkies_allg_init.c

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

#include <gmpxx.h>
#include <math.h>
#include "festkomma.h"

mpf_class y3(0, 128);
mpf_class tmpxx(0, 128), Axx(0, 128), Bxx(0, 128);

/******************************************************************************
 *
 * INIT -- Funktionen (Rechnungen in mpf-class)
 *
 *****************************************************************************/

/* y_0 := (1 - x_0**3)**(1/3). */
void compute_y_value_mpf(mpf_class x_0, mpf_class &y_0)
{
    mpf_class tmp;

    tmp = 1 - x_0 * x_0 * x_0;
    y_0 = pow(tmp.get_d(), 1.0 / 3);
    y_0 = (2 * y_0 + tmp / (y_0 * y_0)) / 3; /* double ist nicht genau genug!
                                                Eine Newton-Iteration. */
}

/* Berechnung von l1, l2 und l3 mit hoher Genauigkeit.
   Genauer wird li[j] = L_i (v_j) ausgerechnet.
   Besonders bei l[2] tritt extreme Ausloeschung auf, weshalb die hohe
   Genauigkeit gebraucht wird.*/
void compute_three_linearf_mpf(mpf_class l[3][3], long v[3][3], mpf_class x_0,
                               mpf_class y_0, double half_step,
                               mpx_t tile_offset, double d, double N)
{
    long i;

    Axx = -x_0 * x_0 / (y_0 * y_0);
    /* Ableitung von
           y(x) := (1 - x**3)**(1/3)
       nach dem Satz ueber implizite Funktionen. */
    /* A = (y_1 - y_0) / (x_1 - x_0); */
    /* A = - x_0*x_0 / pow (1 - x_0*x_0*x_0, 2.0/3.0); */

    Bxx = y_0 - mpx_get_d(tile_offset) - Axx * x_0;
    /* gmp_printf ("x_0 = %.*Ff.\n", 40, x_0.get_mpf_t ());
    gmp_printf ("A = %.*Ff, B = %.*Ff.\n",
                40, Axx.get_mpf_t (), 40, Bxx.get_mpf_t ()); */

    for (i = 0; i < 3; i++)
    {
        l[0][i] = v[i][2] / N; /* l0 (vi)*/
        l[1][i] = (v[i][0] - v[i][2] * x_0) / (half_step * N);
        l[2][i] = (v[i][1] - Axx * v[i][0] - Bxx * v[i][2]) / (d * N);
    }
}

/* Makros fuer lll_mpf. */
#define scal_prod_mpf(prod, l, vec1, vec2)                                     \
    do                                                                         \
    {                                                                          \
        long t1, t2;                                                           \
        mpf_class l1[3], l2[3];                                                \
                                                                               \
        for (t1 = 0; t1 < 3; t1++)                                             \
            l1[t1] = l2[t1] = 0;                                               \
        for (t1 = 0; t1 < 3; t1++)                                             \
            for (t2 = 0; t2 < 3; t2++)                                         \
            {                                                                  \
                l1[t1] += l[t1][t2] * vec1[t2];                                \
                l2[t1] += l[t1][t2] * vec2[t2];                                \
            }                                                                  \
        prod = 0;                                                              \
        for (t1 = 0; t1 < 3; t1++)                                             \
            prod += l1[t1] * l2[t1];                                           \
    } while (0);

#define gram_mpf()                                                             \
    do                                                                         \
    {                                                                          \
        long t1;                                                               \
                                                                               \
        /* printf ("Gram!\n"); */                                              \
        for (t1 = 0; t1 < 3; t1++)                                             \
            vec_gram[0][t1] = vec[0][t1]; /* mpf_class = mpz_class*/           \
        scal_prod_mpf(B[0], l, vec_gram[0], vec_gram[0]);                      \
        scal_prod_mpf(mu[1][0], l, vec[1], vec_gram[0]);                       \
        mu[1][0] /= B[0];                                                      \
                                                                               \
        for (t1 = 0; t1 < 3; t1++)                                             \
            vec_gram[1][t1] = vec[1][t1] - mu[1][0] * vec_gram[0][t1];         \
        scal_prod_mpf(B[1], l, vec_gram[1], vec_gram[1]);                      \
        scal_prod_mpf(mu[2][0], l, vec[2], vec_gram[0]);                       \
        mu[2][0] /= B[0];                                                      \
        scal_prod_mpf(mu[2][1], l, vec[2], vec_gram[1]);                       \
        mu[2][1] /= B[1];                                                      \
                                                                               \
        for (t1 = 0; t1 < 3; t1++)                                             \
            vec_gram[2][t1] = vec[2][t1] - mu[2][0] * vec_gram[0][t1] -        \
                              mu[2][1] * vec_gram[1][t1];                      \
        scal_prod_mpf(B[2], l, vec_gram[2], vec_gram[2]);                      \
    } while (0)

#define red_k_k1_mpf()                                                         \
    do                                                                         \
    {                                                                          \
        /* printf ("Mache RED (%ld, %ld).\n", k, k-1); */                      \
        tm = floor(mu[k][k - 1] + 0.5);                                        \
        q = tm.get_si();                                                       \
        if (!tm.fits_slong_p())                                                \
        {                                                                      \
            printf("q ist zu lang.\n");                                        \
            exit(0);                                                           \
        }                                                                      \
        for (i = 0; i < 3; i++)                                                \
            vec[k][i] -= q * vec[k - 1][i];                                    \
        mu[k][k - 1] -= q;                                                     \
        for (i = 0; i < k - 1; i++)                                            \
            mu[k][i] -= q * mu[k - 1][i];                                      \
    } while (0)

#define red_2_0_mpf()                                                          \
    do                                                                         \
    {                                                                          \
        /* printf ("Mache RED (2, 0).\n"); */                                  \
        tm = floor(mu[2][0] + 0.5);                                            \
        q = tm.get_si();                                                       \
        if (!tm.fits_slong_p())                                                \
        {                                                                      \
            printf("q ist zu lang.\n");                                        \
            exit(0);                                                           \
        }                                                                      \
        for (i = 0; i < 3; i++)                                                \
            vec[2][i] -= q * vec[0][i];                                        \
        mu[2][0] -= q;                                                         \
    } while (0)

#define lll_swap()                                                             \
    do                                                                         \
    {                                                                          \
        /* printf ("Mache SWAP (%ld).\n", k); */                               \
        for (i = 0; i < 3; i++)                                                \
        {                                                                      \
            tmp = vec[k][i];                                                   \
            vec[k][i] = vec[k - 1][i];                                         \
            vec[k - 1][i] = tmp;                                               \
        }                                                                      \
                                                                               \
        muc = mu[k][k - 1];                                                    \
        Bc = B[k] + muc * muc * B[k - 1];                                      \
        mu[k][k - 1] = muc * B[k - 1] / Bc;                                    \
        B[k] = B[k - 1] * B[k] / Bc;                                           \
        B[k - 1] = Bc;                                                         \
                                                                               \
        if (k == 1)                                                            \
        {                                                                      \
            t = mu[2][1];                                                      \
            mu[2][1] = mu[2][0] - muc * t;                                     \
            mu[2][0] = t + mu[1][0] * mu[2][1];                                \
        }                                                                      \
        else                                                                   \
        {/* k == 2*/                                                           \
            tm = mu[2][0];                                                     \
            mu[2][0] = mu[1][0];                                               \
            mu[1][0] = tm;                                                     \
        }                                                                      \
        k = 1;                                                                 \
    } while (0)

/* LLL -- mpf-Version.
   Die Funktion erwartet die drei Linearformen als in l gegeben.
   Sie gibt in v einen Satz von fuer l kurzen Vektoren zurueck.
   Der Ablauf des Algorithmus ist (bis auf einen Fehler in gram) von H. Cohen
   uebernommen. */
inline void lll_mpf(mpf_class l[3][3], long vec[3][3])
{
    long i, k;
    mpf_class B[3], vec_gram[3][3];
    mpf_class mu[3][3];
    long q, tmp;
    mpf_class muc, Bc, t, tm;

    /* Das Standardgitter. */
    vec[0][0] = 1;
    vec[0][1] = 0;
    vec[0][2] = 0;
    vec[1][0] = 0;
    vec[1][1] = 1;
    vec[1][2] = 0;
    vec[2][0] = 0;
    vec[2][1] = 0;
    vec[2][2] = 1;

    /* LLL. */
    k = 1;
    gram_mpf();

    while (k < 3)
    {
        /* printf ("%ld\n", k); */
        red_k_k1_mpf();

        /* Test LLL condition */
        if (B[k] < ((0.99 - mu[k][k - 1] * mu[k][k - 1]) * B[k - 1]))
        {
            /* printf ("Swap mit k = %ld \n", k); */
            lll_swap();
            /* out (); */
        }
        else
        {
            if (k == 2)
                red_2_0_mpf();
            k++;
        }
    }
}

/* Rechnet einen Satz kurzer Vektoren fuer die erste Fliese aus.
   Berechnet auszerdem einen Startwert fuer y_0. */
void init(long v[3][3], mpx_t y_0_xt, mpx_t x_0_xt, double half_step,
          mpx_t tile_offset, double half_tilewidth, double upper_bound)
{
    mpf_set_default_prec(192);
    mpf_class x_0, y_0;
    mpf_class l[3][3];
    long e[3][3];
    long i, j;

    mpf_set_mpx(x_0.get_mpf_t(), x_0_xt);
    /* gmp_printf ("x_0 = %.*Ff.\n", 15, x_0.get_mpf_t ()); */

    compute_y_value_mpf(x_0, y_0);
    /* gmp_printf ("y_0 = %.*Ff.\n\n", 15, y_0.get_mpf_t ()); */
    mpx_set_mpf(y_0_xt, y_0.get_mpf_t());

    /* Standardbasis */
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            e[i][j] = 0;
    for (i = 0; i < 3; i++)
        e[i][i] = 1;

    compute_three_linearf_mpf(l, e, x_0, y_0, half_step, tile_offset,
                              half_tilewidth, upper_bound);

    lll_mpf(l, v);
    /* printf ("v[0] = (%ld, %ld, %ld)\n", v[0][0], v[0][1], v[0][2]);
    printf ("v[1] = (%ld, %ld, %ld)\n", v[1][0], v[1][1], v[1][2]);
    printf ("v[2] = (%ld, %ld, %ld)\n\n", v[2][0], v[2][1], v[2][2]); */

    mpf_set_default_prec(128);
}
