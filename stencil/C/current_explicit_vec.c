/*
  Ab-initio Real-Time Electron Dynamics Simulator, ARTED
  Copyright (C) 2016  ARTED developers

  This file is part of current_explicit_vec.c.

  current_explicit_vec.c is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  current_explicit_vec.c is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with current_explicit_vec.c.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <complex.h>
#include "./interop.h"

extern int NLx, NLy, NLz;

#ifndef ARTED_DOMAIN_POWER_OF_TWO
extern int *modx, *mody, *modz;
#endif

#ifdef ARTED_STENCIL_LOOP_BLOCKING
extern int BX, BY;
#endif

void current_stencil_( double         const* restrict C
                     , double complex const* restrict E
                     , double              * restrict F
                     , double              * restrict G
                     , double              * restrict H
) {
#if defined(__KNC__) || defined(__AVX512F__)
/* Hand-Code Vector processing for Xeon Phi */
  double complex const* e;

  int ix, iy, iz, n;

#ifdef ARTED_STENCIL_LOOP_BLOCKING
  int bx, by;
#endif

  __m512d CONJ = _mm512_set_pd(-1, 1, -1, 1, -1, 1, -1, 1);

  __m512d tt[3];
  for(n = 0 ; n < 3 ; ++n)
    tt[n] = _mm512_setzero_pd();

  __m512d wm[4];
  __m512d wp[4];
  __m512d v1, v2, v3;

  __m512i nly = _mm512_set1_epi32(NLy);
  __m512i nlz = _mm512_set1_epi32(NLz);
  __m512i nyx = _mm512_mask_blend_epi32(0xFF00, _mm512_set1_epi32(NLy    ), _mm512_set1_epi32(NLx    ));
#ifdef ARTED_DOMAIN_POWER_OF_TWO
  __m512i myx = _mm512_mask_blend_epi32(0xFF00, _mm512_set1_epi32(NLy - 1), _mm512_set1_epi32(NLx - 1));
#endif

  __declspec(align(64)) int yx_table_org[16] = { -4, -3, -2, -1, 1, 2, 3, 4, -4, -3, -2, -1, 1, 2, 3, 4};
  __declspec(align(64)) int yx_table[16];
  __m512i *yx_org = (__m512i*) yx_table_org;
  __m512i *yx     = (__m512i*) yx_table;

#pragma noprefetch
#pragma novector
#ifdef ARTED_STENCIL_LOOP_BLOCKING
  for(bx = 0 ; bx < NLx ; bx += BX)
  for(by = 0 ; by < NLy ; by += BY)
  for(ix = bx ; ix < MIN(bx+BX,NLx) ; ++ix)
#else
  for(ix = 0 ; ix < NLx ; ++ix)
#endif
  {
    __m512i tix = _mm512_set1_epi32(ix);
#ifndef ARTED_DOMAIN_POWER_OF_TWO
    __m512i mxm = _mm512_unaligned_load_epi32(modx + (ix - 4 + NLx));
    __m512i mxp = _mm512_alignr_epi32(mxm, mxm, 1);
    __m512i xmp = _mm512_mask_blend_epi32(0xF0F0, mxm, mxp);
            xmp = _mm512_permute4f128_epi32(xmp, _MM_PERM_BADC);
#endif
#ifdef ARTED_STENCIL_LOOP_BLOCKING
  for(iy = by ; iy < MIN(by+BY,NLy) ; ++iy)
#else
  for(iy = 0 ; iy < NLy ; ++iy)
#endif
  {
    __m512i tiy = _mm512_set1_epi32(iy);
#ifndef ARTED_DOMAIN_POWER_OF_TWO
    __m512i mym = _mm512_unaligned_load_epi32(mody + (iy - 4 + NLy));
    __m512i myp = _mm512_alignr_epi32(mym, mym, 1);
    __m512i ymp = _mm512_mask_blend_epi32(0xF0F0, mym, myp);
    __m512i uyx = _mm512_mask_blend_epi32(0xFF00, ymp, xmp);
#endif
    __m512i tyx = _mm512_mask_blend_epi32(0xFF00, tiy, tix);

    e = &E[ix*NLy*NLz + iy*NLz];

    for(iz = 0 ; iz < NLz ; iz += 4)
    {
      __m512i tiz = _mm512_set1_epi32(iz);
#ifdef ARTED_DOMAIN_POWER_OF_TWO
      __m512i mm  = _mm512_sub_epi32(tyx, _mm512_and_epi32(_mm512_add_epi32(_mm512_add_epi32(tyx, nyx), *yx_org), myx));
#else
      __m512i mm  = _mm512_sub_epi32(tyx, uyx);
#endif
      *yx = _mm512_sub_epi32(tiz, _mm512_mullo_epi32(_mm512_mask_mullo_epi32(mm, 0xFF00, mm, nly), nlz));

      // conj(e[iz])
      __m512d ez = _mm512_load_pd(e + iz);
      __m512d w  = _mm512_mul_pd(ez, CONJ);

      /* z-dimension (unit stride) */
      {
        __m512i z0, z2;
#ifdef ARTED_DOMAIN_POWER_OF_TWO
        z0 = _mm512_load_epi64(e + ((iz - 4 + NLz) & (NLz - 1)));
        z2 = _mm512_load_epi64(e + ((iz + 4 + NLz) & (NLz - 1)));
#else
        z0 = _mm512_load_epi64(e + modz[iz - 4 + NLz]);
        z2 = _mm512_load_epi64(e + modz[iz + 4 + NLz]);
#endif

        wm[3] = (__m512d) z0;
        wm[2] = (__m512d) _mm512_alignr_epi32((__m512i) ez, z0,  4);
        wm[1] = (__m512d) _mm512_alignr_epi32((__m512i) ez, z0,  8);
        wm[0] = (__m512d) _mm512_alignr_epi32((__m512i) ez, z0, 12);
        wp[0] = (__m512d) _mm512_alignr_epi32(z2, (__m512i) ez,  4);
        wp[1] = (__m512d) _mm512_alignr_epi32(z2, (__m512i) ez,  8);
        wp[2] = (__m512d) _mm512_alignr_epi32(z2, (__m512i) ez, 12);
        wp[3] = (__m512d) z2;

        v2 = _mm512_setzero_pd();
#pragma unroll(4)
        for(n = 0 ; n < 4 ; ++n) {
          v1 = _mm512_sub_pd(wp[n], wm[n]);
          v2 = _mm512_fmadd_pd(_mm512_set1_pd(C[n+8]), v1, v2);
        }
        v3    = (__m512d) _mm512_shuffle_epi32((__m512i) v2, _MM_PERM_BADC);
        tt[2] = _mm512_fmadd_pd(w, v3, tt[2]);
      }

      /* y-dimension (NLz stride) */
      {
        wm[3] = _mm512_load_pd(e + yx_table[0]);
        wm[2] = _mm512_load_pd(e + yx_table[1]);
        wm[1] = _mm512_load_pd(e + yx_table[2]);
        wm[0] = _mm512_load_pd(e + yx_table[3]);
        wp[0] = _mm512_load_pd(e + yx_table[4]);
        wp[1] = _mm512_load_pd(e + yx_table[5]);
        wp[2] = _mm512_load_pd(e + yx_table[6]);
        wp[3] = _mm512_load_pd(e + yx_table[7]);

        v2 = _mm512_setzero_pd();
#pragma unroll(4)
        for(n = 0 ; n < 4 ; ++n) {
          v1 = _mm512_sub_pd(wp[n], wm[n]);
          v2 = _mm512_fmadd_pd(_mm512_set1_pd(C[n+4]), v1, v2);
        }
        v3    = (__m512d) _mm512_shuffle_epi32((__m512i) v2, _MM_PERM_BADC);
        tt[1] = _mm512_fmadd_pd(w, v3, tt[1]);
      }

      /* x-dimension (NLy*NLz stride)  */
      {
        wm[3] = _mm512_load_pd(e + yx_table[ 8]);
        wm[2] = _mm512_load_pd(e + yx_table[ 9]);
        wm[1] = _mm512_load_pd(e + yx_table[10]);
        wm[0] = _mm512_load_pd(e + yx_table[11]);
        wp[0] = _mm512_load_pd(e + yx_table[12]);
        wp[1] = _mm512_load_pd(e + yx_table[13]);
        wp[2] = _mm512_load_pd(e + yx_table[14]);
        wp[3] = _mm512_load_pd(e + yx_table[15]);

        v2 = _mm512_setzero_pd();
#pragma unroll(4)
        for(n = 0 ; n < 4 ; ++n) {
          v1 = _mm512_sub_pd(wp[n], wm[n]);
          v2 = _mm512_fmadd_pd(_mm512_set1_pd(C[n]), v1, v2);
        }
        v3    = (__m512d) _mm512_shuffle_epi32((__m512i) v2, _MM_PERM_BADC);
        tt[0] = _mm512_fmadd_pd(w, v3, tt[0]);
      }
    }  /* NLz */
  } /* NLy */
  } /* NLx */

  *F = _mm512_reduce_add_pd(tt[0]);
  *G = _mm512_reduce_add_pd(tt[1]);
  *H = _mm512_reduce_add_pd(tt[2]);
#else
/* Hand-Code Vector processing for Intel Ivy-Bridge (AVX) */
  double complex const* e;

  int ix, iy, iz, n;

#ifdef ARTED_STENCIL_LOOP_BLOCKING
  int bx, by;
#endif

  const __m256i CONJ = _mm256_set_epi64x(1LL << 63, 0, 1LL << 63, 0);

  __m256d tx = _mm256_setzero_pd();
  __m256d ty = _mm256_setzero_pd();
  __m256d tz = _mm256_setzero_pd();

#pragma novector
#ifdef ARTED_STENCIL_LOOP_BLOCKING
  for(bx = 0 ; bx < NLx ; bx += BX)
  for(by = 0 ; by < NLy ; by += BY)
  for(ix = bx ; ix < MIN(bx+BX,NLx) ; ++ix) {
  for(iy = by ; iy < MIN(by+BY,NLy) ; ++iy)
#else
  for(ix = 0 ; ix < NLx ; ++ix) {
  for(iy = 0 ; iy < NLy ; ++iy)
#endif
  {
    e = &E[ix*NLy*NLz + iy*NLz];

    for(iz = 0 ; iz < NLz ; iz += 2)
    {
#define STENCIL_CALC(MM,PP,CC) \
      v2 = _mm256_mul_pd(_mm256_broadcast_sd(CC),  _mm256_sub_pd(PP, MM)); \
      v1 = _mm256_add_pd(v1, v2);

#define STENCIL(IDM,IDP,CC) \
      m = _mm256_load_pd((double *) (e + iz - IDM));     \
      p = _mm256_load_pd((double *) (e + iz - IDP));     \
      STENCIL_CALC(m,p,CC) \

      __m256d m, p;
      __m256d v1, v2, v3, v4;

      // conj(e[iz])
      __m256d ez = _mm256_load_pd((double *)(e + iz));
      __m256d w = (__m256d) _mm256_xor_si256((__m256i) ez, CONJ);

      /* x-dimension (NLy*NLz stride)  */
      {
        v1 = _mm256_setzero_pd();
        STENCIL(IDX(-1), IDX(1), C+0);
        STENCIL(IDX(-2), IDX(2), C+1);
        STENCIL(IDX(-3), IDX(3), C+2);
        STENCIL(IDX(-4), IDX(4), C+3);
        v3 = _mm256_shuffle_pd(v1, v1, 0x05);
        v4 = _mm256_mul_pd(w, v3);
        tx = _mm256_add_pd(v4, tx);
      }

      /* y-dimension (NLz stride) */
      {
        v1 = _mm256_setzero_pd();
        STENCIL(IDY(-1), IDY(1), C+4);
        STENCIL(IDY(-2), IDY(2), C+5);
        STENCIL(IDY(-3), IDY(3), C+6);
        STENCIL(IDY(-4), IDY(4), C+7);
        v3 = _mm256_shuffle_pd(v1, v1, 0x05);
        v4 = _mm256_mul_pd(w, v3);
        ty = _mm256_add_pd(v4, ty);
      }

      /* z-dimension (unit stride) */
      {
        __m256d z0,z1,z2,z3,z4,z5,z6,z7;
#ifdef ARTED_DOMAIN_POWER_OF_TWO
        z1 = _mm256_load_pd((double *)(e + ((iz - 2 + NLz) & (NLz - 1))));
        z2 = _mm256_load_pd((double *)(e + ((iz + 2 + NLz) & (NLz - 1))));
        z0 = _mm256_load_pd((double *)(e + ((iz - 4 + NLz) & (NLz - 1))));
        z3 = _mm256_load_pd((double *)(e + ((iz + 4 + NLz) & (NLz - 1))));
#else
        z1 = _mm256_load_pd((double *)(e + modz[iz - 2 + NLz]));
        z2 = _mm256_load_pd((double *)(e + modz[iz + 2 + NLz]));
        z0 = _mm256_load_pd((double *)(e + modz[iz - 4 + NLz]));
        z3 = _mm256_load_pd((double *)(e + modz[iz + 4 + NLz]));
#endif
        z6 = _mm256_permute2f128_pd(z0, z1, 0x21);
        z4 = _mm256_permute2f128_pd(z1, ez, 0x21);
        z5 = _mm256_permute2f128_pd(ez, z2, 0x21);
        z7 = _mm256_permute2f128_pd(z2, z3, 0x21);

        v1 = _mm256_setzero_pd();
        STENCIL_CALC(z1, z2, C+ 9);
        STENCIL_CALC(z0, z3, C+11);
        STENCIL_CALC(z4, z5, C+ 8);
        STENCIL_CALC(z6, z7, C+10);
        v3 = _mm256_shuffle_pd(v1, v1, 0x05);
        v4 = _mm256_mul_pd(w, v3);
        tz = _mm256_add_pd(v4, tz);
      }
    }  /* NLz */
  } /* NLy */
  } /* NLx */

  __declspec(align(64)) double r[4];

  _mm256_store_pd(r, tx);
  *F = r[0] + r[1] + r[2] + r[3];

  _mm256_store_pd(r, ty);
  *G = r[0] + r[1] + r[2] + r[3];

  _mm256_store_pd(r, tz);
  *H = r[0] + r[1] + r[2] + r[3];
#endif /* if defined(__KNC__) || defined(__AVX512F__) */
}
