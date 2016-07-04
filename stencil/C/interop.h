/*
  Ab-initio Real-Time Electron Dynamics Simulator, ARTED
  Copyright (C) 2016  ARTED developers

  This file is part of interop.h.

  interop.h is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  interop.h is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with interop.h.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef ARTED_INTEROP
#define ARTED_INTEROP

/* Stencil computation code with C supports Intel compiler only. */

#if defined(__KNC__) || defined(__AVX512F__)
# define MEM_ALIGNED 64
# define VECTOR_SIZE 4
#else
# define MEM_ALIGNED 32
# define VECTOR_SIZE 2
#endif

#ifdef ARTED_DOMAIN_POWER_OF_TWO
# ifdef ENABLE_STENCIL_CODE_WITH_PADDING
#   define IDX(npt) ((ix - ((ix + (npt) + NLx) & (NLx-1))) * PNLy * PNLz)
#   define IDY(npt) ((iy - ((iy + (npt) + NLy) & (NLy-1))) * PNLz)
#   define IDZ(npt) ((iz - ((iz + (npt) + NLz) & (NLz-1))))
# else
#   define IDX(npt) ((ix - ((ix + (npt) + NLx) & (NLx-1))) * NLy * NLz)
#   define IDY(npt) ((iy - ((iy + (npt) + NLy) & (NLy-1))) * NLz)
#   define IDZ(npt) ((iz - ((iz + (npt) + NLz) & (NLz-1))))
# endif
#else
# ifdef ENABLE_STENCIL_CODE_WITH_PADDING
#   define IDX(npt) ((ix - modx[ix + (npt) + NLx]) * PNLy * PNLz)
#   define IDY(npt) ((iy - mody[iy + (npt) + NLy]) * PNLz)
#   define IDZ(npt) ((iz - modz[iz + (npt) + NLz]))
# else
#   define IDX(npt) ((ix - modx[ix + (npt) + NLx]) * NLy * NLz)
#   define IDY(npt) ((iy - mody[iy + (npt) + NLy]) * NLz)
#   define IDZ(npt) ((iz - modz[iz + (npt) + NLz]))
# endif
#endif /* ARTED_DOMAIN_POWER_OF_TWO */

#define MIN(n,m) (n < m ? n : m)

#include <immintrin.h>

#ifdef ARTED_STENCIL_LOOP_BLOCKING
#define BX   opt_variables_mp_stencil_blocking_x_
#define BY   opt_variables_mp_stencil_blocking_y_
#endif

#define NL   global_variables_mp_nl_
#define NLx  global_variables_mp_nlx_
#define NLy  global_variables_mp_nly_
#define NLz  global_variables_mp_nlz_
#define PNLx opt_variables_mp_pnlx_
#define PNLy opt_variables_mp_pnly_
#define PNLz opt_variables_mp_pnlz_

#define modx opt_variables_mp_modx_
#define mody opt_variables_mp_mody_
#define modz opt_variables_mp_modz_

#define IDXF opt_variables_mp_zifdt_

#ifdef ARTED_MIGRATE_TO_KNL
# include "./knc2knl.h"
#else
# ifdef __KNC__
/* Knights Corner */
inline
__m512i _mm512_unaligned_load_epi32(int const* v) {
  __m512i w;
  w = _mm512_loadunpacklo_epi32(w, v + 0);
  w = _mm512_loadunpackhi_epi32(w, v + 16);
  return w;
}

inline
__m512d dcast_to_dcmplx(double const *v) {
  __m512d w = _mm512_loadunpacklo_pd(_mm512_setzero_pd(), v);
  __m512i a = _mm512_permute4f128_epi32((__m512i) w, _MM_PERM_BBAA);
  __m512d b = (__m512d) _mm512_shuffle_epi32(a, _MM_SHUFFLE(1,0,3,2));
  return _mm512_mask_blend_pd(0x66, (__m512d) a, b);
}
# endif
#endif

#ifdef __KNC__
/* Knights Corner */
inline
__m512d dcomplex_mul(__m512d a, __m512d b) {
  __m512d ze = _mm512_setzero_pd();
  __m512d s0 = _mm512_swizzle_pd(a, _MM_SWIZ_REG_CDAB);  /* s0 = [a.i a.r] */
  __m512d re = _mm512_mask_blend_pd(0xAA, a, s0);        /* re = [a.r a.r] */
  __m512d im = _mm512_mask_sub_pd(a, 0x55, ze, s0);      /* im = [-a.i a.i] */
  __m512d t0 = _mm512_mul_pd(re, b);                     /* t0 = [a.r*b.r a.r*b.i] */
  __m512d s1 = _mm512_swizzle_pd(b, _MM_SWIZ_REG_CDAB);  /* s1 = [b.i b.r] */
  return       _mm512_fmadd_pd(im, s1, t0);              /* [-a.i*b.i+a.r*b.r a.i*b.r+a.r*b.i] */
}

inline
__m512d dcomplex_gather(void const* v, __m512i idx) {
  /* swizzle index */
  __m512i one = _mm512_set_epi32(1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0);

  __m512i a   = _mm512_swizzle_epi64(idx, _MM_SWIZ_REG_AAAA);  // ABABABABABABABAB
  __m512i b   = _mm512_swizzle_epi64(idx, _MM_SWIZ_REG_BBBB);  // CDCDCDCDCDCDCDCD
  __m512i ab  = _mm512_shuffle_epi32(a, _MM_SHUFFLE(3,1,2,0)); // AABBAABBAABBAABB
  __m512i cd  = _mm512_shuffle_epi32(b, _MM_SHUFFLE(3,1,2,0)); // CCDDCCDDCCDDCCDD

  __m512i x   = _mm512_mask_blend_epi32(0xF0F0, ab, cd);       // AABBCCDDAABBCCDD
  __m512i y   = _mm512_slli_epi32(x, 1);
  __m512i z   = _mm512_xor_si512(y, one);

  /* gather */
  return _mm512_i32loextgather_pd(z, v, _MM_UPCONV_PD_NONE, 8, _MM_HINT_NONE);
}
#endif /* ifdef __KNC__ */

#if defined(__AVX__) && !defined(__AVX512F__)
/* Sandy-Bridge or higher processors */
inline
__m256d dcast_to_dcmplx(double const *v) {
  __m256d a = _mm256_loadu_pd(v);
  __m256d b = _mm256_permute2f128_pd(a, a, 0x0);
  return _mm256_shuffle_pd(b, b, 0xC);
}
#endif /* defined(__AVX__) && !defined(__AVX512F__) */

#endif /* ARTED_INTEROP */
