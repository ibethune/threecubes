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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <gmp.h>
#include <math.h>
#include "festkomma.h"

#define UPPER_BOUND 1.0e14
#define LOWER_BOUND 1.0e11

/* Triples with |x**3 + y**3 - z**3| < MAX_K are allowed. */
#define MAX_K 1000

/* The global variables half_step, half_tilewidth und
   tile_offset are calculated in compute_tile_params().

   This function will perform all NEW_TILES iterations of the
   outermost loop.
   The FACTOR controls the size of the tile.

   NB:
   half_step is the half-dimension of a tile in the x-direction.
   half_tilewidth is the half-dimension of a tile in the y-direction.
   It is defined as:
      half_tilewidth = 1.001 * y''(x) * half_step**2 / 4
   wthere
      y(x) := (1 - x**3)**(1/3).
   The area of the tile is thus 1.001 * y''(x) * half_step**3.

   This value is set as 1.001 * (FACTOR/UPPER_BOUND)**3.
   Thus FACTOR, UPPER_BOUND and the current second derivative control the
   half_step und half_tilewidth calculation. */
#define NEW_TILES 1000000
#define FACTOR 5.5
double half_step, half_tilewidth;
mpx_t x_0, tile_offset, Ax;
mpx_t y_diff, y_inv, y_inv_diff;
char output[1000], file[100];
long row;

/* Counter of the total number of processed tiles. */
long tiles;

/* These mpf variables should be locals, mainly in the
   functions compute_y_value and compute_three_linearf or
   compute_three_linearf_mpf.
   The code runs faster with global variables. */
mpf_t tmp1, tmp2;

/* Just a prototype */
void init(long v[3][3], mpx_t y_0_type, mpx_t x_0_type, double half_step,
          mpx_t tile_offset, double half_tilewidth, double upper_bound);

void out()
{
    FILE *fp;

    fp = fopen(file, "a");
    fprintf(fp, "%6ld \t ", row);
    row++;
    fprintf(fp, output);
    fclose(fp);
}

/* pr := m1 * m2. Product of two 3x3 matrices with integer entries.
   We compute in long integer and watch for overflows \pm 2**63.
   pr = m1 or pr = m2 are allowed. */
inline long matrix_prod(long pr[3][3], long m1[3][3], long m2[3][3])
{
    long i, j;
    long prod[3][3];
    long prodh, argh, argl;
    long err;

    err = 0;
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            smul_ppmm(prodh, prod[i][j], m1[i][0], m2[0][j]);
            smul_ppmm(argh, argl, m1[i][1], m2[1][j]);
            addf_ssaaaa(prodh, prod[i][j], prodh, prod[i][j], argh, argl);
            smul_ppmm(argh, argl, m1[i][2], m2[2][j]);
            addf_ssaaaa(prodh, prod[i][j], prodh, prod[i][j], argh, argl);
            if (((prod[i][j] >= 0) && (prodh != 0)) ||
                ((prod[i][j] < 0) && (prodh != -1)))
            {
                /* sprintf (output, "Overflow during matrix product. "); out ();
                   sprintf (output, " prod[i][j] = %ld, prodh = %ld.\n",
                                                          prod[i][j], prodh);
                out ();
                sprintf (output, "m1 = [%ld %ld %ld], [%ld %ld %ld], [%ld %ld
                %ld]\n",
                     m1[0][0], m1[0][1], m1[0][2],
                     m1[1][0], m1[1][1], m1[1][2],
                     m1[2][0], m1[2][1], m1[2][2]); out ();
                sprintf (output, "m2 = [%ld %ld %ld], [%ld %ld %ld], [%ld %ld
                %ld]\n",
                     m2[0][0], m2[0][1], m2[0][2],
                     m2[1][0], m2[1][1], m2[1][2],
                     m2[2][0], m2[2][1], m2[2][2]); out (); */
                err = 1; /* Overflow! */
            }
        }
    }
    /* Copyback - needed for the case that pr is aliased to m1 or m2 */
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            pr[i][j] = prod[i][j];
        }
    }
    return (err);
}

/* Compute half_step and half_tilewidth at the current
   point x_0 \in [0,1]. */
void compute_tile_params(mpx_t step)
{
    double x, second_deriv, padding;
    long do_output;

    /* Output only every 100*NEW_TILES tiles. */
    do_output = tiles % (100 * NEW_TILES);

    x = mpx_get_d(x_0);

    if (do_output == 0)
    {
        sprintf(output, "\n");
        out();
        sprintf(output, "x_0 = %.15f.\n", x);
        out();
    }

    second_deriv = 2 * x / pow(1 - x * x * x, 5.0 / 3.0);
    half_step = pow(second_deriv, -1.0 / 3.0) * FACTOR / UPPER_BOUND;

    padding = MAX_K / (LOWER_BOUND * LOWER_BOUND * LOWER_BOUND);
    half_tilewidth =
        1.001 * second_deriv * half_step * half_step / 4 + padding;
    /* padding chosen so that solutions of
       size >LOWER_BOUND are guaranteed. */

    if (do_output == 0)
    {
        sprintf(output, "half_step = %.9e, half_tilewidth = %.9e.\n", half_step,
                half_tilewidth);
        out();
    }

    mpx_set_d(step, 2 * half_step);
    mpx_set_d(tile_offset, second_deriv * half_step * half_step / 4);
}

/***************************************************************************
 *
 * Actual computational routines...
 *
 **************************************************************************/

/* Called once every NEW_TILES tiles. Wird einmal pro NEW_TILES Fliesen aufgerufen.
   Computes 1/y_0 -1 - expensive and exact. */
void y_inv_init(mpx_t y_0)
{
    /* Initialisation of y_inv. */
    mpf_set_ui(tmp1, 1);
    mpf_set_mpx(tmp2, y_0);
    mpf_div(tmp1, tmp1, tmp2);
    mpf_sub_ui(tmp1, tmp1, 1);
    mpx_set_mpf(y_inv, tmp1);
    /* gmp_sprintf (output, "y_inv = %.*Ff.\n", 40, tmp1); */
}

/* Called once every NEW_TILES tiles.
   Initialises y_diff to y'(x_0) * step.
   Computes y_inv_diff as y_diff / y_0**2. */
void y_diff_init(mpx_t y_0, mpx_t step)
{
    double diff, y;

    /* Initialisation of y_diff. */
    mpx_mul(y_diff, Ax, step);
    /* If y_0 < x_0 then Ax > 1. Correct it here. */
    if ((y_0[1] < x_0[1]) || ((y_0[1] == x_0[1]) && (y_0[0] <= x_0[0])))
        mpx_add(y_diff, y_diff, step);
    /* mpf_set_mpx (tmp1, y_diff);
    gmp_sprintf (output, "y_diff = %.*Ff.\n\n", 40, tmp1); */

    /* Initialisation of y_inv_diff. */
    y = mpx_get_d(y_0);
    diff = mpx_get_d(y_diff) / (y * y);
    mpx_set_d(y_inv_diff, diff);
    /* mpf_set_mpx (tmp1, y_inv_diff);
    gmp_sprintf (output, "y_inv_diff = %.*Ff.\n", 40, tmp1); */
}

/* y_0 := (1 - x_0**3)**(1/3). */
inline void compute_y_value(mpx_t y_0)
{
    double nenner, diff;
    mpx_t tmpx1, tmpx2, tmpx3;

    mpx_add(y_0, y_0, y_diff);
    /* y_0 -= step * A. 
       Linear improvements of the startin value.
       Minus, because we run backwards through the interval.
       y_diff s calculated once per 1000000 (FLEISEN_NEU) tiles.
       A is always negative, so (-step)*A > 0.
       NB: y_diff > 0 after initialisation. */

    mpx_mul(tmpx1, x_0, x_0);
    mpx_mul(tmpx2, tmpx1, x_0);
    tmpx3[0] = ~tmpx2[0];
    tmpx3[1] = ~tmpx2[1];
    /* y3 = 1 - x_0*x_0*x_0; */
    /* We make an error here of exactly 2**(-128). */

    mpx_mul(tmpx1, y_0, y_0);
    nenner = 3.0 * mpx_get_d(tmpx1);

    mpx_mul(tmpx2, y_0, tmpx1);
    mpx_sub(tmpx1, tmpx2, tmpx3);
    diff = mpx_get_d(tmpx1);
    if (diff > 0.5)
    {
        sprintf(output, "Underflow.\n");
        out();
        exit(0);
    }
    diff /= nenner;

    mpx_set_d(tmpx1, diff);
    mpx_sub(y_0, y_0, tmpx1);
    /* One Newton Iteration. */
    /* y_0 = (2*y_0 + y3/(y_0*y_0)) / 3
           = (2*y_0*y_0*y_0 + y3) / (3*y_0*y_0)
           = y_0 + (y3 - y_0*y_0*y_0) / (3*y_0*y_0);
       Difference can be computed in double precision. */
    /* One Newton Iteration. */

    /* mpf_set_mpx (tmp1, y_0);
    gmp_sprintf (output, "y_0 = %.*Ff.\n", 40, tmp1); out (); */
}

/* Store y'(x_0) in global variable Ax. */
inline void compute_y_deriv(mpx_t y_0)
{
    mpx_t tmpx1, tmpx2;

    /* First compute 1/y_0. */

    /* Linear improvement of starting value */
    mpx_sub(y_inv, y_inv, y_inv_diff);

    /* Now Newton Iteration inv_new = inv + (1 - y_0 * inv) * inv.
    Since y_inv is 1/y - 1, we must compute:
         y_inv = y_inv + (1 - y_0 * (1 + y_inv)) * (1 + y_inv)
               = y_inv + (1 - y_0 - y_0 * y_inv) + (1 - y_0 - y_0 * y_inv) *
                 y_inv */
    mpx_mul(tmpx1, y_0, y_inv);
    mpx_add(tmpx2, tmpx1, y_0); /* tmpx2 ist y_0 + y_0 * y_inv. */
    tmpx2[0] = ~tmpx2[0];
    tmpx2[1] = ~tmpx2[1];
    /* tmpx2 is now 1 - y_0 - y_0 * y_inv. Error of 2**(-128) is ignored. */
    mpx_mul(tmpx1, tmpx2, y_inv);
    mpx_add(tmpx1, tmpx1, tmpx2);
    /* tmpx1 is (1 - y_0 - y_0 * y_inv) + (1 - y_0 - y_0 * y_inv) * y_inv. */
    mpx_add(y_inv, y_inv, tmpx1);
    /* mpf_set_mpx (tmp1, y_inv);
    gmp_sprintf (output, "y_inv = %.*Ff.\n", 40, tmp1); */

    /* Now compute
            y'(x_0) = Ax := x_0**2 / y_0**2
                    = x_0**2 * (1 + y_inv)**2
                    = x_0**2 + x_0**2 * (y_inv**2 + 2*y_inv). */
    mpx_mul(tmpx1, y_inv, y_inv);
    mpx_add(tmpx1, tmpx1, y_inv);
    mpx_add(tmpx1, tmpx1, y_inv); /* tmpx1 is now y_inv**2 + 2*y_inv. */
    mpx_mul(tmpx2, x_0, x_0);     /* tmpx2 is x_0**2. */
    mpx_mul(Ax, tmpx2, tmpx1); /* Ax is now x_0**2 * (y_inv**2 + 2*y_inv). */
    mpx_add(Ax, Ax, tmpx2);
    /* mpf_set_mpx (tmp1, Ax);
    gmp_sprintf (output, "A = %.*Ff.\n", 40, tmp1); */

    /* Ax = - x_0*x_0 / (y_0*y_0); */
    /* Derivation of
           y(x) := (1 - x**3)**(1/3)
       according to the theorem of implicit functions. */
    /* A = (y_1 - y_0) / (x_1 - x_0); */
    /* A = - x_0*x_0 / pow (1 - x_0*x_0*x_0, 2.0/3.0); */
}

inline void compute_three_linearf(double l[3][3], long v[3][3], mpx_t y_0,
                                  double d, double N)
{
    long i;
    mpx_t tmpx1, tmpx2, tmpx3, Bx;

    /* mpf_set_mpx (tmp1, x_0);
    gmp_sprintf (output, "\nx_0 = %.*Ff.\n", 40, tmp1); out (); */
    mpx_mul(tmpx1, Ax, x_0);
    /* If y_0 < x_0 then Ax > 1. Correct it here. */
    if ((y_0[1] < x_0[1]) || ((y_0[1] == x_0[1]) && (y_0[0] <= x_0[0])))
        mpx_add(tmpx1, tmpx1, x_0);

    /* mpf_set_mpx (tmp1, tile_offset);
    gmp_sprintf (output, "tile_offset = %.*Ff.\n", 40, tmp1); out (); */
    mpx_sub(tmpx2, y_0, tile_offset);
    mpx_add(Bx, tmpx1, tmpx2); /* Ax is always negative. */
    /* Bx = y_0 - tile_offset - Ax * x_0; */
    /* 1 <= Bx <= 1.6. Actually store Bx - 1. */
    /* mpf_set_mpx (tmp1, Ax); mpf_set_mpx (tmp2, Bx);
    gmp_sprintf (output, "A = %.*Ff.\nB = %.*Ff.\n", 40, tmp1, 40, tmp2);
    out(); */

    for (i = 0; i < 3; i++)
    {
        l[0][i] = v[i][2] / N;
        /* sprintf (output, "l[0][%ld] = %f.\n", i, l[0][i]); out (); */
        /* No problem with accuracy here */

        mpx_mul_si(tmpx1, x_0, -v[i][2]);
        mpxb_add_si(tmpx2, tmpx1, v[i][0]);
        /* tmpx = v[i][0] - v[i][2] * x_0; */
        l[1][i] = mpxg_get_d_simple(tmpx2);
        l[1][i] /= (half_step * N);
        /* sprintf (output, "l[1][%ld] = %f.\n", i, l[1][i]); out (); */
        /* step*N \approx 100.
           So we need tmp2 to >2 (7?) decimal places.
           Double precision should suffice.
           tmp2 is 128 Bit, also >=23 decimal places of precision. */

        mpx_mul_si(tmpx1, Ax, v[i][0]);
        /* If y_0 < x_0 then Ax > 1. Correct it here. */
        if ((y_0[1] < x_0[1]) || ((y_0[1] == x_0[1]) && (y_0[0] <= x_0[0])))
            mpxb_add_si(tmpx1, tmpx1, v[i][0]);

        mpx_mul_si(tmpx2, Bx, -v[i][2]);
        mpx_add(tmpx3, tmpx1, tmpx2); /* Ax is always negative. */
        mpxb_add_si(tmpx1, tmpx3, v[i][1] - v[i][2]); /* Bx \in [1..2]. */
        /* tmpx = v[i][1] - Ax * v[i][0] - Bx * v[i][2]; */
        l[2][i] = mpxg_get_d(tmpx1);
        l[2][i] /= (d * N);
        /* sprintf (output, "l[2][%ld] = %f.\n", i, l[2][i]); out (); */
        /* l[2][i] is very sensitive to lack of precision.
           d*N \approx 1.0e-15.
           So we need tmp2 to > 15 (20?) decimal places.
           Also Ax and Bx to 35 decimal places, corresponding to 128 Bit. */
    }
}

/* Macros for LLL. */
#define scal_prod(prod, l, vec1, vec2)                                         \
    do                                                                         \
    {                                                                          \
        long t1, t2;                                                           \
        double l1[3], l2[3];                                                   \
                                                                               \
        for (t1 = 0; t1 < 3; t1++)                                             \
            l1[t1] = l2[t1] = 0;                                               \
        for (t1 = 0; t1 < 3; t1++)                                             \
        {                                                                      \
            for (t2 = 0; t2 < 3; t2++)                                         \
            {                                                                  \
                l1[t1] += l[t1][t2] * vec1[t2];                                \
                l2[t1] += l[t1][t2] * vec2[t2];                                \
            }                                                                  \
        }                                                                      \
        prod = 0;                                                              \
        for (t1 = 0; t1 < 3; t1++)                                             \
            prod += l1[t1] * l2[t1];                                           \
    } while (0);

#define gram()                                                                 \
    do                                                                         \
    {                                                                          \
        long t1;                                                               \
                                                                               \
        /* sprintf (output, "Gram!\n"); out (); */                             \
        for (t1 = 0; t1 < 3; t1++)                                             \
            vec_gram[0][t1] = vec[0][t1]; /* double = long */                  \
        scal_prod(B[0], l, vec_gram[0], vec_gram[0]);                          \
        scal_prod(mu[1][0], l, vec[1], vec_gram[0]);                           \
        mu[1][0] /= B[0];                                                      \
                                                                               \
        for (t1 = 0; t1 < 3; t1++)                                             \
            vec_gram[1][t1] = vec[1][t1] - mu[1][0] * vec_gram[0][t1];         \
        scal_prod(B[1], l, vec_gram[1], vec_gram[1]);                          \
        scal_prod(mu[2][0], l, vec[2], vec_gram[0]);                           \
        mu[2][0] /= B[0];                                                      \
        scal_prod(mu[2][1], l, vec[2], vec_gram[1]);                           \
        mu[2][1] /= B[1];                                                      \
                                                                               \
        for (t1 = 0; t1 < 3; t1++)                                             \
            vec_gram[2][t1] = vec[2][t1] - mu[2][0] * vec_gram[0][t1] -        \
                              mu[2][1] * vec_gram[1][t1];                      \
        scal_prod(B[2], l, vec_gram[2], vec_gram[2]);                          \
    } while (0)

#define red_k_k1()                                                             \
    do                                                                         \
    {                                                                          \
        /* sprintf (output, "Do RED (%ld, %ld).\n", k, k-1); out (); */     \
        tm = floor(mu[k][k - 1] + 0.5);                                        \
        /* if (fabs (tm) > (1LU << 63) - 1) {                                  \
         sprintf (output, "q is too long.\n"); out ();                         \
         exit (0);                                                             \
        } */                                                                   \
        q = lround(tm);                                                        \
        if (q != 0)                                                            \
        {                                                                      \
            for (i = 0; i < 3; i++)                                            \
                vec[k][i] -= q * vec[k - 1][i];                                \
            mu[k][k - 1] -= tm;                                                \
            if (k == 2)                                                        \
                mu[2][0] -= tm * mu[1][0];                                     \
        }                                                                      \
    } while (0)

#define red_2_0()                                                              \
    do                                                                         \
    {                                                                          \
        /* sprintf (output, "Do RED (2, 0).\n"); out (); */                 \
        tm = floor(mu[2][0] + 0.5);                                            \
        q = lround(tm);                                                        \
        /* if (fabs (tm) > (1LU << 63) - 1) {                                  \
         sprintf (output, "q is too long.\n"); out ();                         \
         exit (0);                                                             \
        } */                                                                   \
        if (q != 0)                                                            \
            for (i = 0; i < 3; i++)                                            \
                vec[2][i] -= q * vec[0][i];                                    \
    } while (0)

#define lll_swap()                                                             \
    do                                                                         \
    {                                                                          \
        /* sprintf (output, "Do SWAP (%ld).\n", k); out (); */              \
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

/* LLL -- double-version.
   The function expects the three linear forms to be given in l.
   It returns in v a set of l short vectors(?).

   The steps of the algorithm are (apart from an error in gram() ) taken from H. Cohen.

   We apply LLL, with l computed already in a particular basis.
   Wir wendn lll in der Weise an, dasz l schon bezueglich einer bestimmten
   Basis ausgerechnet ist. The returned matrix array is the transformation
   matrix from to the new basis. 

   WARNING: This function works only if l consists of not too large numbers.
   We step along the curve and hope that the vectors from the last tile that
   were short, are now not very long. */
inline void lll(double l[3][3], long vec[3][3])
{
    long i, k;
    double B[3], vec_gram[3][3];
    double mu[3][3];
    long q, tmp;
    double muc, Bc, t, tm;

    /* Standard basis. Initial guess for the transformation matrix. */
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
    gram();

    while (k < 3)
    {
        /* sprintf (output, "%ld\n", k); out (); */
        red_k_k1();

        /* Test LLL condition */
        if (B[k] < ((0.99 - mu[k][k - 1] * mu[k][k - 1]) * B[k - 1]))
        {
            /* if (B[k] < 0.5 * B[k-1]) { */
            /* sprintf (output, "Swap mit k = %ld \n", k); out (); */
            lll_swap();
            /* out (); */
        }
        else
        {
            if (k == 2)
                red_2_0();
            k++;
        }
    }
}

#define lf_new(lf, l, e)                                                       \
    do                                                                         \
    {                                                                          \
        /* long  t1, t2;                                                       \
        for (t1 = 0; t1 < 3; t1++)                                             \
         for (t2 = 0; t2 < 3; t2++)                                            \
          lf[t1][t2]                                                           \
          = l[t1][0]*e[t2][0] + l[t1][1]*e[t2][1] + l[t1][2]*e[t2][2]; */      \
                                                                               \
        lf[0][0] = l[0][0] * e[0][0] + l[0][1] * e[0][1] + l[0][2] * e[0][2];  \
        lf[0][1] = l[0][0] * e[1][0] + l[0][1] * e[1][1] + l[0][2] * e[1][2];  \
        lf[0][2] = l[0][0] * e[2][0] + l[0][1] * e[2][1] + l[0][2] * e[2][2];  \
                                                                               \
        lf[1][0] = l[1][0] * e[0][0] + l[1][1] * e[0][1] + l[1][2] * e[0][2];  \
        lf[1][1] = l[1][0] * e[1][0] + l[1][1] * e[1][1] + l[1][2] * e[1][2];  \
        lf[1][2] = l[1][0] * e[2][0] + l[1][1] * e[2][1] + l[1][2] * e[2][2];  \
                                                                               \
        lf[2][0] = l[2][0] * e[0][0] + l[2][1] * e[0][1] + l[2][2] * e[0][2];  \
        lf[2][1] = l[2][0] * e[1][0] + l[2][1] * e[1][1] + l[2][2] * e[1][2];  \
        lf[2][2] = l[2][0] * e[2][0] + l[2][1] * e[2][1] + l[2][2] * e[2][2];  \
    } while (0);

/******************************************************************************
 *
 * The code for finding grid points in the pyramid.
 *
 *****************************************************************************/

/* res = m * v, 3x3-Matrix times column-vector. */
#define MMUL(res, m, v)                                                        \
    do                                                                         \
    {                                                                          \
        /* long  k;                                                            \
                                                                               \
        for (k = 0; k < 3; k++)                                                \
         res[k] = m[k][0]*v[0] + m[k][1]*v[1] + m[k][2]*v[2]; */               \
                                                                               \
        res[0] = m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2];             \
        res[1] = m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2];             \
        res[2] = m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2];             \
    } while (0)

/* Compute the four corners of the pyramid. */
inline void pyr_corners(double pt1[3], double pt2[3], double pt3[3],
                      double pt4[3], double lf[3][3])
{
    long i, j;
    double det_inv, adj[3][3], inv[3][3];

    /* The corners are given by (l1, l2, l3) (pti) = vi.
       Important: v1, ..., v4 are listed in order of circulation. */
    double v1[3] = {1.0, 1.0, 1.0};
    double v2[3] = {1.0, -1.0, 1.0};
    double v3[3] = {1.0, -1.0, -1.0};
    /* double      v4[3] = {1.0,  1.0, -1.0}; */

    /* Inverse Matrix (unscaled). */
    adj[0][0] = lf[1][1] * lf[2][2] - lf[1][2] * lf[2][1];
    adj[0][1] = -(lf[0][1] * lf[2][2] - lf[0][2] * lf[2][1]);
    adj[0][2] = lf[0][1] * lf[1][2] - lf[0][2] * lf[1][1];

    adj[1][0] = -(lf[1][0] * lf[2][2] - lf[1][2] * lf[2][0]);
    adj[1][1] = lf[0][0] * lf[2][2] - lf[0][2] * lf[2][0];
    adj[1][2] = -(lf[0][0] * lf[1][2] - lf[0][2] * lf[1][0]);

    adj[2][0] = lf[1][0] * lf[2][1] - lf[1][1] * lf[2][0];
    adj[2][1] = -(lf[0][0] * lf[2][1] - lf[0][1] * lf[2][0]);
    adj[2][2] = lf[0][0] * lf[1][1] - lf[0][1] * lf[1][0];

    /* Determinant */
    det_inv = 1.0 / (lf[0][0] * adj[0][0] + lf[0][1] * adj[1][0] +
                     lf[0][2] * adj[2][0]);
    /* sprintf (output, "det_inv = %.15f\n", det_inv); out (); */
    /* Inverse */
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            inv[i][j] = adj[i][j] * det_inv;

    MMUL(pt1, inv, v1);
    MMUL(pt2, inv, v2);
    MMUL(pt3, inv, v3);
    /* MMUL (pt4, inv, v4); */
    pt4[0] = pt1[0] + pt3[0] - pt2[0];
    pt4[1] = pt1[1] + pt3[1] - pt2[1];
    pt4[2] = pt1[2] + pt3[2] - pt2[2];

    /* sprintf (output, "Corners:\n"); out ();
    sprintf (output, "pt1 = (%f, %f, %f)\n", pt1[0], pt1[1], pt1[2]); out ();
    sprintf (output, "pt2 = (%f, %f, %f)\n", pt2[0], pt2[1], pt2[2]); out ();
    sprintf (output, "pt3 = (%f, %f, %f)\n", pt3[0], pt3[1], pt3[2]); out ();
    sprintf (output, "pt4 = (%f, %f, %f)\n", pt4[0], pt4[1], pt4[2]); out (); */
}

/* The outermost loop runs over z.
   So need the maximum and minimum of the z-coordinates of the
   5 corners of the pyramid. */
inline void z_bounds(long *z_start, long *z_end, double *p1, double *p2,
                        double *p3, double *p4)
{
    double tmp1, tmp2;

    tmp1 = MIN(MIN(MIN(p1[2], p2[2]), MIN(p3[2], p4[2])), 0.0);
    tmp2 = MAX(MAX(MAX(p1[2], p2[2]), MAX(p3[2], p4[2])), 0.0);

    *z_start = lround(tmp1 + 0.5); /* Round up. */
    *z_end = lround(tmp2 - 0.5); /* Round down. */
    /* Whole numbers may be rounded incorrectly.
    This is not a bug, however, since we only need the points inside
    the pyramid. */

    /* sprintf (output, "z-bounds: [%ld, %ld]\n", *z_start, *z_end); out (); */
}

/* One edge is represntated as an array of 6 values.
   y = edge[2]*z + edge[3] and
   x = edge[4]*z + edge[5]
   for edge[0] <= z <= edge[1]. */
inline void edge(double *edge, double *start, double *end)
{
    double dz_inv;

    dz_inv = 1 / (start[2] - end[2]);

    if (start[2] < end[2])
    {
        edge[0] = start[2];
        edge[1] = end[2];
    }
    else
    {
        edge[0] = end[2];
        edge[1] = start[2];
    }
    /* If this causes a division by 0, then edge[0] = edge[1] anyway
    (and this is hopefully not an integer). */
    edge[2] = (start[1] - end[1]) * dz_inv; /* dy/dz */
    edge[3] = start[1] - edge[2] * start[2];

    edge[4] = (start[0] - end[0]) * dz_inv; /* dx/dz */
    edge[5] = start[0] - edge[4] * start[2];

    /* sprintf (output, "Edge: y = %f*z + %f, x = %f*z + %f for %f <= z <=
    %f\n",
                  edge[2], edge[3], edge[4], edge[5], edge[0], edge[1]);
    out (); */
}

/* Compute the eight edges of the pyramid */
inline void pyr_edges(double edges[8][6], double *p1, double *p2, double *p3,
                       double *p4)
{
    double s[3] = {0.0, 0.0, 0.0};

    edge(edges[0], s, p1);
    edge(edges[1], s, p2);
    edge(edges[2], s, p3);
    edge(edges[3], s, p4);
    edge(edges[4], p1, p2);
    edge(edges[5], p2, p3);
    edge(edges[6], p3, p4);
    edge(edges[7], p4, p1);
}

/* The two inner loops run over x und y.
   So we cut the pyramid with the z = ... plane and have a polygon.
   Polygon. We need the maximum and minimum of the x and y coordinates
   of this polygon.
   Corners are created by the intersection of the plane with
   the edges */
inline void x_and_y_bounds(long *x_start, long *x_end, long *y_start,
                              long *y_end, double edge[8][6], long z)
{
    int i;
    double *kan, k;
    double xstart, xend, ystart, yend;

    xstart = LONG_MAX;
    xend = LONG_MIN;
    ystart = LONG_MAX;
    yend = LONG_MIN;

    /* Loop over the eight edges */
    for (i = 0; i < 8; i++)
    {
        kan = edge[i]; /* this edge */
        /* Do we intersect this edge at all? */
        if ((kan[0] <= z) && (z <= kan[1]))
        {
            k = kan[2] * z + kan[3];
            ystart = MIN(ystart, k);
            yend = MAX(yend, k);
            k = kan[4] * z + kan[5];
            xstart = MIN(xstart, k);
            xend = MAX(xend, k);
        }
    }
    *y_start = lround(ystart + 0.5); /* Round up. */
    *y_end = lround(yend - 0.5); /* Round down. */
    *x_start = lround(xstart + 0.5); /* Round up. */
    *x_end = lround(xend - 0.5); /* Round down. */

    /* sprintf (output, "x_y-bounds for z = %ld: [%ld, %ld] x [%ld, %ld]\n",
                                            z, *x_start, *x_end, *y_start, *y_end);
    out (); */
}

/* Output the solutions. */
void post_proc(long vec_out0, long vec_out1, long vec_out2, long f0)
{
    /* double  quot, x0, x1; */

    /* quot = ((double) vec_out0) / ((double) vec_out2); */
    /* sprintf (output, "%.18f\n", quot); out (); */

    /* x0 = mpx_get_d (x_0);
    x0 -= half_step;
    x1 = x0 + 2 * half_step; */
    /* sprintf (output, "Interval [%.18f,%.18f]\n\n", x0, x1); out (); */

    /* if ((x0 < quot) && (x1 > quot)) */
    {
        sprintf(output, "(%ld, %ld, %ld) Solution for k = %ld.\n", vec_out0,
                vec_out1, vec_out2, f0);
        out();
        sprintf(output, "\n");
        out();
        /* sprintf (output, "Solution found.\n"); out (); */
        /* exit (0); */
    }
}

/* What is v2**3 - v0**3 - v1**3 for the current values of z, y and x?
   All computation is modulo 2**64.
   Output triples (v0, v1, v2), where v0**3 + v1**3 - v2**3 equals \pm 3(???). */
inline ulong evaluate_function(long x, long v00, long v01, long v02,
                               long vec_h0, long vec_h1, long vec_h2)
{
    long vec_out0, vec_out1, vec_out2;
    ulong f0;

    /* sprintf (output, "x = %ld.\n", x); out (); */

    vec_out0 = x * v00 + vec_h0; /* Use ulong arithmetic here? */
    vec_out1 = x * v01 + vec_h1; /* Don't have any overflows, anyway. */
    vec_out2 = x * v02 + vec_h2;

    /* Compute modulo 2**64. Therefore, cast to unsigned. */
    f0 = ((ulong)vec_out2) * ((ulong)vec_out2) * ((ulong)vec_out2) -
         ((ulong)vec_out0) * ((ulong)vec_out0) * ((ulong)vec_out0) -
         ((ulong)vec_out1) * ((ulong)vec_out1) * ((ulong)vec_out1);

    /* small_vect() computes only the case z >= 0.
    Which is why we need both plus and minus at this point.
    Output, of course, as signed long ints. */
    if ((f0 < MAX_K) || (-f0 < MAX_K))
        if (f0 != 0)
            post_proc(vec_out0, vec_out1, vec_out2, f0); /* Output */

    return (f0);
}

#define XY_LOOP                                                                \
    {                                                                          \
        for (y = y_start; y <= y_end; y++)                                     \
        {                                                                      \
            /* sprintf (output, "Compute y = %ld.\n", y); out (); */           \
            vec_h0 = y * v[1][0] + z * v[2][0];                                \
            vec_h1 = y * v[1][1] + z * v[2][1];                                \
            vec_h2 = y * v[1][2] + z * v[2][2];                                \
                                                                               \
            x = x_start;                                                       \
            /* These function calls causes output if f0 or -f0 is              \
               less than MAX_K. */                                             \
            f0 = evaluate_function(x, v00, v01, v02, vec_h0, vec_h1, vec_h2);  \
            ans++;                                                             \
                                                                               \
            x++;                                                               \
            f1 = evaluate_function(x, v00, v01, v02, vec_h0, vec_h1, vec_h2);  \
            ans++;                                                             \
                                                                               \
            x++;                                                               \
            f2 = evaluate_function(x, v00, v01, v02, vec_h0, vec_h1, vec_h2);  \
            ans++;                                                             \
                                                                               \
            x++;                                                               \
                                                                               \
            /* Difference scheme.*/                                            \
            d1 = f1 - f0;                                                      \
            d = f2 - f1;                                                       \
            dd = d - d1;                                                       \
            f = f2;                                                            \
            dd += ddd; /* increment dd one step. */                            \
            for (; x <= x_end; x++)                                            \
            {                                                                  \
                d += dd;                                                       \
                dd += ddd;                                                     \
                f += d;                                                        \
                /* if (f != evaluate_function (x, v00, v01, v02, vec_h0,       \
                vec_h1, vec_h2)) {                                             \
                 sprintf (output, "%lu %lu\n", f,                              \
                             evaluate_function (x, v00, v01, v02, vec_h0,      \
                vec_h1, vec_h2));                                              \
                 out ();                                                       \
                 sprintf (output, "Bug in Difference scheme!\n"); out ();      \
                 exit (0);                                                     \
                } */                                                           \
                ans++;                                                         \
                if ((f < MAX_K) || (-f < MAX_K))                               \
                    /* Only for output, compute the function value again.      \
                    Is faster than making an extra functioncall. */            \
                    evaluate_function(x, v00, v01, v02, vec_h0, vec_h1,        \
                                      vec_h2);                                 \
            }                                                                  \
        }                                                                      \
    }                                                                          \
    while (0)                                                                  \
        ;

/* Search for whole vectors
        0 < l1 < 1, |l2| < l1 und |l3| < l1.
   All calculations take place with respect to the reduced basis v.
   If (vec[0], vec[1], vec[2]) are divisible by 3, we do not call post_proc().
   (Because vec[0]**3 + vec[1]**3 + vec[2]**3 = 0 (mod 3) either all three
   components are divisible by 3 or none at all.) */
inline long small_vect(double lf[3][3], long v[3][3])
{
    long x, y, z;
    long x_start, x_end;
    long y_start, y_end;
    long z_start, z_end;
    long vec_h0, vec_h1, vec_h2;
    long ans;
    long v00, v01, v02;
    double p1[3], p2[3], p3[3], p4[3];
    double edge[8][6];
    ulong f, f0, f1, f2;
    ulong d, dd, ddd, d1;

    /* sprintf (output, "lll-Basis:\n[%ld %ld %ld],\n[%ld %ld %ld],\n[%ld %ld
    %ld]\n",
              v[0][0], v[0][1], v[0][2],
              v[1][0], v[1][1], v[1][2],
              v[2][0], v[2][1], v[2][2]);
    out ();
    sprintf (output, "Linear forms of the first vector: %f %f %f\n",
                                                  lf[0][0], lf[1][0], lf[2][0]);
    out (); */

    /* Compute the corners of the pyramid. The peak is the origin. */
    pyr_corners(p1, p2, p3, p4, lf);
    /* Compute the bound for z. */
    z_bounds(&z_start, &z_end, p1, p2, p3, p4);

    /* Compute the eight edges of the pyramid. */
    pyr_edges(edge, p1, p2, p3, p4);

    /* Triple-nested loop through the pyramid. */
    ans = 0;
    v00 = v[0][0];
    v01 = v[0][1];
    v02 = v[0][2];
    /* Difference scheme, first part.
    ddd is independnt of z, constant over the whole tile. */
    ddd = ((ulong)6) * (((ulong)v02) * ((ulong)v02) * ((ulong)v02) -
                        ((ulong)v00) * ((ulong)v00) * ((ulong)v00) -
                        ((ulong)v01) * ((ulong)v01) * ((ulong)v01));

    for (z = z_start; z <= z_end; z++)
    {
        /* Compute bounds for x and y. */
        x_and_y_bounds(&x_start, &x_end, &y_start, &y_end, edge, z);
        XY_LOOP;
    }
    return (ans);
}

/* Initialise the outer loop. */
void compute_interval(mpx_t x_0_start, mpx_t x_0_end)
{
    mpx_t step;
    mpx_t y_0;
    double lf[3][3], ln[3][3];
    long e[3][3], v[3][3];
    long counter, ans, err;

    x_0[0] = x_0_end[0];
    x_0[1] = x_0_end[1];
    /* Run backwards from x_0_end to x_0_start. */

    compute_tile_params(step);

    /* First loop-pass. Here we compute LLL with high precision.
       y_0 is also initialised. */
    init(v, y_0, x_0, half_step, tile_offset, half_tilewidth, UPPER_BOUND);
    y_inv_init(y_0);
    y_diff[0] = 0;
    y_diff[1] = 0;
    y_inv_diff[0] = 0;
    y_inv_diff[1] = 0;
    /* y_diff = y_inv_diff = 0. */
    sprintf(output, "Initialisation complete.\n");
    out();

    /* Loop. Here almost everything is done in double precision. */
    counter = NEW_TILES;
    tiles = -NEW_TILES;
    err = 0;
    /* while (x_0 >= x_0_start) ... . */
    while ((x_0[1] > x_0_start[1]) ||
           ((x_0[1] == x_0_start[1]) && (x_0[0] >= x_0_start[0])))
    {
        /* mpf_set_mpx (tmp1, x_0);
        gmp_sprintf (output, "\nx_0 = %.*Ff.\n", 40, tmp1); out (); */

        compute_y_value(y_0);
        /* The first time y_diff = 0 => y_0 is computed correctly. */
        /* mpf_set_mpx (tmp1, y_0);
        gmp_sprintf (output, "y_0 = %.*Ff.\n\n", 40, tmp1); out (); */
        compute_y_deriv(y_0);
        /* The first time y_inv_diff = 0 => y'(x_0) is computed correctly. */
        /* mpf_set_mpx (tmp1, Ax);
        gmp_sprintf (output, "A = %.*Ff.\n\n", 40, tmp1); out (); */

        if (counter == NEW_TILES)
        {
            counter = 0;
            tiles += NEW_TILES;
            compute_tile_params(step);
            y_inv_init(y_0);
            y_diff_init(y_0, step);
        }

        /* Last v was not correct due to overflow. */
        if (err > 0)
        {
            init(v, y_0, x_0, half_step, tile_offset, half_tilewidth,
                 UPPER_BOUND);
            mpf_set_mpx(tmp1, x_0);
            gmp_sprintf(output, "Neustart bei x_0 = %.*Ff.\n", 25, tmp1);
            out();
        }

        compute_three_linearf(lf, v, y_0, half_tilewidth, UPPER_BOUND);

        lll(lf, e);
        err = matrix_prod(v, e,
                          v); /* v is definitely modulo 2**64. */

        lf_new(ln, lf, e);
        ans = small_vect(ln, v);
        /* sprintf (output, "Find %ld Gridpoints.\n", ans); out (); */

        /* Tile goes left and right from x_0. */
        mpx_sub(x_0, x_0, step);
        counter++;
    }

    sprintf(output, "%ld Tiles computed total.\n", tiles + counter);
    out();
}

/* Run with something like:
      elkies_allg 0.4 0.000001

   Where 0.4 is the initial value.
   0.000001 is the interval length that should be processed. */
int main(int argc, char *argv[])
{
    mpx_t diff, x_0_start, x_0_end;

    mpf_set_default_prec(128);
    mpf_init(tmp1);
    mpf_init(tmp2);

    mpf_set_str(tmp1, argv[1], 10);
    mpx_set_mpf(x_0_start, tmp1);
    mpf_set_str(tmp2, argv[2], 10);
    mpx_set_mpf(diff, tmp2);

    row = 1;
    mpx_add(x_0_end, x_0_start, diff);

    mpf_set_mpx(tmp1, x_0_start);
    mpf_set_mpx(tmp2, x_0_end);
    /* Name of the output file */
    gmp_sprintf(file, "liste_allg_%.*Ff_%.*Ff.txt", 8, tmp1, 8, tmp2);

    /* First output */
    gmp_sprintf(output, "Start computing from %.*Ff to %.*Ff.\n", 15, tmp1, 15,
                tmp2);
    out();

    compute_interval(x_0_start, x_0_end);

    mpf_clear(tmp1);
    mpf_clear(tmp2);
    return 0;
}
