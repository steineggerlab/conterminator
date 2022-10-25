/* SPDX-License-Identifier: MIT
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * Copyright:
 *   2021      Evan Nemerson <evan@nemerson.com>
 */

#if !defined(SIMDE_WASM_SIMD128_H)
#define SIMDE_WASM_SIMD128_H

#include "../simde-common.h"

HEDLEY_DIAGNOSTIC_PUSH
SIMDE_DISABLE_UNWANTED_DIAGNOSTICS
SIMDE_BEGIN_DECLS_

typedef union {
  #if defined(SIMDE_VECTOR_SUBSCRIPT)
    SIMDE_ALIGN_TO_16 int8_t          i8 SIMDE_VECTOR(16) SIMDE_MAY_ALIAS;
    SIMDE_ALIGN_TO_16 int16_t        i16 SIMDE_VECTOR(16) SIMDE_MAY_ALIAS;
    SIMDE_ALIGN_TO_16 int32_t        i32 SIMDE_VECTOR(16) SIMDE_MAY_ALIAS;
    SIMDE_ALIGN_TO_16 int64_t        i64 SIMDE_VECTOR(16) SIMDE_MAY_ALIAS;
    SIMDE_ALIGN_TO_16 uint8_t         u8 SIMDE_VECTOR(16) SIMDE_MAY_ALIAS;
    SIMDE_ALIGN_TO_16 uint16_t       u16 SIMDE_VECTOR(16) SIMDE_MAY_ALIAS;
    SIMDE_ALIGN_TO_16 uint32_t       u32 SIMDE_VECTOR(16) SIMDE_MAY_ALIAS;
    SIMDE_ALIGN_TO_16 uint64_t       u64 SIMDE_VECTOR(16) SIMDE_MAY_ALIAS;
    #if defined(SIMDE_HAVE_INT128_)
    SIMDE_ALIGN_TO_16 simde_int128  i128 SIMDE_VECTOR(16) SIMDE_MAY_ALIAS;
    SIMDE_ALIGN_TO_16 simde_uint128 u128 SIMDE_VECTOR(16) SIMDE_MAY_ALIAS;
    #endif
    SIMDE_ALIGN_TO_16 simde_float32  f32 SIMDE_VECTOR(16) SIMDE_MAY_ALIAS;
    SIMDE_ALIGN_TO_16 simde_float64  f64 SIMDE_VECTOR(16) SIMDE_MAY_ALIAS;
    SIMDE_ALIGN_TO_16 int_fast32_t  i32f SIMDE_VECTOR(16) SIMDE_MAY_ALIAS;
    SIMDE_ALIGN_TO_16 uint_fast32_t u32f SIMDE_VECTOR(16) SIMDE_MAY_ALIAS;
  #else
    SIMDE_ALIGN_TO_16 int8_t         i8[16];
    SIMDE_ALIGN_TO_16 int16_t        i16[8];
    SIMDE_ALIGN_TO_16 int32_t        i32[4];
    SIMDE_ALIGN_TO_16 int64_t        i64[2];
    SIMDE_ALIGN_TO_16 uint8_t        u8[16];
    SIMDE_ALIGN_TO_16 uint16_t       u16[8];
    SIMDE_ALIGN_TO_16 uint32_t       u32[4];
    SIMDE_ALIGN_TO_16 uint64_t       u64[2];
    #if defined(SIMDE_HAVE_INT128_)
    SIMDE_ALIGN_TO_16 simde_int128  i128[1];
    SIMDE_ALIGN_TO_16 simde_uint128 u128[1];
    #endif
    SIMDE_ALIGN_TO_16 simde_float32  f32[4];
    SIMDE_ALIGN_TO_16 simde_float64  f64[2];
    SIMDE_ALIGN_TO_16 int_fast32_t  i32f[16 / sizeof(int_fast32_t)];
    SIMDE_ALIGN_TO_16 uint_fast32_t u32f[16 / sizeof(uint_fast32_t)];
  #endif

  #if defined(SIMDE_X86_SSE_NATIVE)
    SIMDE_ALIGN_TO_16 __m128         sse_m128;
    #if defined(SIMDE_X86_SSE2_NATIVE)
      SIMDE_ALIGN_TO_16 __m128i      sse_m128i;
      SIMDE_ALIGN_TO_16 __m128d      sse_m128d;
    #endif
  #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
    SIMDE_ALIGN_TO_16 int8x16_t      neon_i8;
    SIMDE_ALIGN_TO_16 int16x8_t      neon_i16;
    SIMDE_ALIGN_TO_16 int32x4_t      neon_i32;
    SIMDE_ALIGN_TO_16 int64x2_t      neon_i64;
    SIMDE_ALIGN_TO_16 uint8x16_t     neon_u8;
    SIMDE_ALIGN_TO_16 uint16x8_t     neon_u16;
    SIMDE_ALIGN_TO_16 uint32x4_t     neon_u32;
    SIMDE_ALIGN_TO_16 uint64x2_t     neon_u64;
    SIMDE_ALIGN_TO_16 float32x4_t    neon_f32;
    #if defined(SIMDE_ARM_NEON_A64V8_NATIVE)
      SIMDE_ALIGN_TO_16 float64x2_t    neon_f64;
    #endif
  #elif defined(SIMDE_WASM_SIMD128_NATIVE)
    SIMDE_ALIGN_TO_16 v128_t         wasm_v128;
  #elif defined(SIMDE_POWER_ALTIVEC_P6_NATIVE) || defined(SIMDE_ZARCH_ZVECTOR_13_NATIVE)
    SIMDE_ALIGN_TO_16 SIMDE_POWER_ALTIVEC_VECTOR(unsigned char)      altivec_u8;
    SIMDE_ALIGN_TO_16 SIMDE_POWER_ALTIVEC_VECTOR(unsigned short)     altivec_u16;
    SIMDE_ALIGN_TO_16 SIMDE_POWER_ALTIVEC_VECTOR(unsigned int)       altivec_u32;
    SIMDE_ALIGN_TO_16 SIMDE_POWER_ALTIVEC_VECTOR(signed char)        altivec_i8;
    SIMDE_ALIGN_TO_16 SIMDE_POWER_ALTIVEC_VECTOR(signed short)       altivec_i16;
    SIMDE_ALIGN_TO_16 SIMDE_POWER_ALTIVEC_VECTOR(signed int)         altivec_i32;
    SIMDE_ALIGN_TO_16 SIMDE_POWER_ALTIVEC_VECTOR(float)              altivec_f32;
    #if defined(SIMDE_POWER_ALTIVEC_P7_NATIVE) || defined(SIMDE_ZARCH_ZVECTOR_13_NATIVE)
      SIMDE_ALIGN_TO_16 SIMDE_POWER_ALTIVEC_VECTOR(unsigned long long) altivec_u64;
      SIMDE_ALIGN_TO_16 SIMDE_POWER_ALTIVEC_VECTOR(signed long long)   altivec_i64;
      SIMDE_ALIGN_TO_16 SIMDE_POWER_ALTIVEC_VECTOR(double)             altivec_f64;
    #endif
  #endif
} simde_v128_private;

#if defined(SIMDE_WASM_SIMD128_NATIVE)
  typedef v128_t simde_v128_t;
#elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
   typedef int32x4_t simde_v128_t;
#elif defined(SIMDE_X86_SSE2_NATIVE)
   typedef __m128i simde_v128_t;
#elif defined(SIMDE_X86_SSE_NATIVE)
   typedef __m128 simde_v128_t;
#elif defined(SIMDE_POWER_ALTIVEC_P6_NATIVE) || defined(SIMDE_ZARCH_ZVECTOR_13_NATIVE)
   typedef SIMDE_POWER_ALTIVEC_VECTOR(signed int) simde_v128_t;
#elif defined(SIMDE_VECTOR_SUBSCRIPT)
  typedef int32_t simde_v128_t SIMDE_ALIGN_TO_16 SIMDE_VECTOR(16) SIMDE_MAY_ALIAS;
#else
  typedef simde_v128_private simde_v128_t;
#endif

#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  typedef simde_v128_t v128_t;
#endif

HEDLEY_STATIC_ASSERT(16 == sizeof(simde_v128_t), "simde_v128_t size incorrect");
HEDLEY_STATIC_ASSERT(16 == sizeof(simde_v128_private), "simde_v128_private size incorrect");
#if defined(SIMDE_CHECK_ALIGNMENT) && defined(SIMDE_ALIGN_OF)
HEDLEY_STATIC_ASSERT(SIMDE_ALIGN_OF(simde_v128_t) == 16, "simde_v128_t is not 16-byte aligned");
HEDLEY_STATIC_ASSERT(SIMDE_ALIGN_OF(simde_v128_private) == 16, "simde_v128_private is not 16-byte aligned");
#endif

#define SIMDE_WASM_SIMD128_GENERATE_CONVERSION_FUNCTIONS(Other_Type, SIMDe_Type, To_Name, From_Name) \
  SIMDE_FUNCTION_ATTRIBUTES \
  Other_Type To_Name(SIMDe_Type v) { \
    Other_Type r; \
    simde_memcpy(&r, &v, sizeof(r)); \
    return r; \
  } \
  \
  SIMDE_FUNCTION_ATTRIBUTES \
  SIMDe_Type From_Name(Other_Type v) { \
    SIMDe_Type r; \
    simde_memcpy(&r, &v, sizeof(r)); \
    return r; \
  }

SIMDE_WASM_SIMD128_GENERATE_CONVERSION_FUNCTIONS(simde_v128_private, simde_v128_t, simde_v128_to_private, simde_v128_from_private)

#if defined(SIMDE_X86_SSE2_NATIVE)
  SIMDE_WASM_SIMD128_GENERATE_CONVERSION_FUNCTIONS(__m128 , simde_v128_t, simde_v128_to_m128 , simde_v128_from_m128 )
  SIMDE_WASM_SIMD128_GENERATE_CONVERSION_FUNCTIONS(__m128i, simde_v128_t, simde_v128_to_m128i, simde_v128_from_m128i)
  SIMDE_WASM_SIMD128_GENERATE_CONVERSION_FUNCTIONS(__m128d, simde_v128_t, simde_v128_to_m128d, simde_v128_from_m128d)
#endif

#if defined(SIMDE_ARM_NEON_A32V7_NATIVE)
  SIMDE_WASM_SIMD128_GENERATE_CONVERSION_FUNCTIONS(  int8x16_t, simde_v128_t, simde_v128_to_neon_i8 , simde_v128_from_neon_i8 )
  SIMDE_WASM_SIMD128_GENERATE_CONVERSION_FUNCTIONS(  int16x8_t, simde_v128_t, simde_v128_to_neon_i16, simde_v128_from_neon_i16)
  SIMDE_WASM_SIMD128_GENERATE_CONVERSION_FUNCTIONS(  int32x4_t, simde_v128_t, simde_v128_to_neon_i32, simde_v128_from_neon_i32)
  SIMDE_WASM_SIMD128_GENERATE_CONVERSION_FUNCTIONS(  int64x2_t, simde_v128_t, simde_v128_to_neon_i64, simde_v128_from_neon_i64)
  SIMDE_WASM_SIMD128_GENERATE_CONVERSION_FUNCTIONS( uint8x16_t, simde_v128_t, simde_v128_to_neon_u8 , simde_v128_from_neon_u8 )
  SIMDE_WASM_SIMD128_GENERATE_CONVERSION_FUNCTIONS( uint16x8_t, simde_v128_t, simde_v128_to_neon_u16, simde_v128_from_neon_u16)
  SIMDE_WASM_SIMD128_GENERATE_CONVERSION_FUNCTIONS( uint32x4_t, simde_v128_t, simde_v128_to_neon_u32, simde_v128_from_neon_u32)
  SIMDE_WASM_SIMD128_GENERATE_CONVERSION_FUNCTIONS( uint64x2_t, simde_v128_t, simde_v128_to_neon_u64, simde_v128_from_neon_u64)
  SIMDE_WASM_SIMD128_GENERATE_CONVERSION_FUNCTIONS(float32x4_t, simde_v128_t, simde_v128_to_neon_f32, simde_v128_from_neon_f32)
  #if defined(SIMDE_ARM_NEON_A64V8_NATIVE)
    SIMDE_WASM_SIMD128_GENERATE_CONVERSION_FUNCTIONS(float64x2_t, simde_v128_t, simde_v128_to_neon_f64, simde_v128_from_neon_f64)
  #endif
#endif /* defined(SIMDE_ARM_NEON_A32V7_NATIVE) */

#if defined(SIMDE_POWER_ALTIVEC_P6_NATIVE) || defined(SIMDE_ZARCH_ZVECTOR_13_NATIVE)
  SIMDE_WASM_SIMD128_GENERATE_CONVERSION_FUNCTIONS(SIMDE_POWER_ALTIVEC_VECTOR(  signed  char), simde_v128_t, simde_v128_to_altivec_i8 , simde_v128_from_altivec_i8 )
  SIMDE_WASM_SIMD128_GENERATE_CONVERSION_FUNCTIONS(SIMDE_POWER_ALTIVEC_VECTOR(  signed short), simde_v128_t, simde_v128_to_altivec_i16, simde_v128_from_altivec_i16)
  SIMDE_WASM_SIMD128_GENERATE_CONVERSION_FUNCTIONS(SIMDE_POWER_ALTIVEC_VECTOR(  signed   int), simde_v128_t, simde_v128_to_altivec_i32, simde_v128_from_altivec_i32)
  SIMDE_WASM_SIMD128_GENERATE_CONVERSION_FUNCTIONS(SIMDE_POWER_ALTIVEC_VECTOR(unsigned  char), simde_v128_t, simde_v128_to_altivec_u8 , simde_v128_from_altivec_u8 )
  SIMDE_WASM_SIMD128_GENERATE_CONVERSION_FUNCTIONS(SIMDE_POWER_ALTIVEC_VECTOR(unsigned short), simde_v128_t, simde_v128_to_altivec_u16, simde_v128_from_altivec_u16)
  SIMDE_WASM_SIMD128_GENERATE_CONVERSION_FUNCTIONS(SIMDE_POWER_ALTIVEC_VECTOR(unsigned   int), simde_v128_t, simde_v128_to_altivec_u32, simde_v128_from_altivec_u32)
  #if defined(SIMDE_POWER_ALTIVEC_P7_NATIVE) || defined(SIMDE_ZARCH_ZVECTOR_13_NATIVE)
    SIMDE_WASM_SIMD128_GENERATE_CONVERSION_FUNCTIONS(SIMDE_POWER_ALTIVEC_VECTOR(  signed  long), simde_v128_t, simde_v128_to_altivec_i64, simde_v128_from_altivec_i64)
    SIMDE_WASM_SIMD128_GENERATE_CONVERSION_FUNCTIONS(SIMDE_POWER_ALTIVEC_VECTOR(unsigned  long), simde_v128_t, simde_v128_to_altivec_u64, simde_v128_from_altivec_u64)
  #endif

  #if defined(SIMDE_BUG_GCC_95782)
    SIMDE_FUNCTION_ATTRIBUTES
    SIMDE_POWER_ALTIVEC_VECTOR(float)
    simde_v128_to_altivec_f32(simde_v128_t value) {
      simde_v128_private r_ = simde_v128_to_private(value);
      return r_.altivec_f32;
    }

    SIMDE_FUNCTION_ATTRIBUTES
    simde_v128_t
    simde_v128_from_altivec_f32(SIMDE_POWER_ALTIVEC_VECTOR(float) value) {
      simde_v128_private r_;
      r_.altivec_f32 = value;
      return simde_v128_from_private(r_);
    }
  #else
    SIMDE_WASM_SIMD128_GENERATE_CONVERSION_FUNCTIONS(SIMDE_POWER_ALTIVEC_VECTOR(float), simde_v128_t, simde_v128_to_altivec_f32, simde_v128_from_altivec_f32)
  #endif
#endif /* defined(SIMDE_POWER_ALTIVEC_P6_NATIVE) */

/*
 * Begin function implementations
 */

/* load */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_v128_load (const void * mem) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_v128_load(mem);
  #elif defined(SIMDE_X86_SSE2_NATIVE)
    return _mm_loadu_si128(HEDLEY_REINTERPRET_CAST(const __m128i*, mem));
  #else
    simde_v128_t r;
    simde_memcpy(&r, mem, sizeof(r));
    return r;
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_v128_load(mem) simde_wasm_v128_load((mem))
#endif

/* store */

SIMDE_FUNCTION_ATTRIBUTES
void
simde_wasm_v128_store (void * mem, simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    wasm_v128_store(mem, a);
  #elif defined(SIMDE_X86_SSE2_NATIVE)
    _mm_storeu_si128(HEDLEY_REINTERPRET_CAST(__m128i*, mem), a);
  #else
    simde_memcpy(mem, &a, sizeof(a));
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_v128_store(mem, a) simde_wasm_v128_store((mem), (a))
#endif

/* make */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i8x16_make (
    int8_t c0, int8_t c1, int8_t  c2, int8_t  c3, int8_t  c4, int8_t  c5, int8_t  c6, int8_t  c7,
    int8_t c8, int8_t c9, int8_t c10, int8_t c11, int8_t c12, int8_t c13, int8_t c14, int8_t c15) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return
      wasm_i8x16_make(
        c0, c1,  c2,  c3,  c4,  c5,  c6,  c7,
        c8, c9, c10, c11, c12, c13, c14, c15);
  #elif defined(SIMDE_X86_SSE2_NATIVE)
    return
      _mm_setr_epi8(
        c0, c1,  c2,  c3,  c4,  c5,  c6,  c7,
        c8, c9, c10, c11, c12, c13, c14, c15);
  #else
    simde_v128_private r_;

    r_.i8[ 0] =  c0;
    r_.i8[ 1] =  c1;
    r_.i8[ 2] =  c2;
    r_.i8[ 3] =  c3;
    r_.i8[ 4] =  c4;
    r_.i8[ 5] =  c5;
    r_.i8[ 6] =  c6;
    r_.i8[ 7] =  c7;
    r_.i8[ 8] =  c8;
    r_.i8[ 9] =  c9;
    r_.i8[10] = c10;
    r_.i8[11] = c11;
    r_.i8[12] = c12;
    r_.i8[13] = c13;
    r_.i8[14] = c14;
    r_.i8[15] = c15;

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define \
    wasm_i8x16_make( \
        c0, c1,  c2,  c3,  c4,  c5,  c6,  c7, \
        c8, c9, c10, c11, c12, c13, c14, c15) \
      simde_wasm_i8x16_make( \
          (c0), (c1),  (c2),  (c3),  (c4),  (c5),  (c6),  (c7), \
          (c8), (c9), (c10), (c11), (c12), (c13), (c14), (c15))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_make (
    int16_t c0, int16_t c1, int16_t  c2, int16_t  c3, int16_t  c4, int16_t  c5, int16_t  c6, int16_t  c7) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_make(c0, c1, c2, c3, c4, c5, c6, c7);
  #elif defined(SIMDE_X86_SSE2_NATIVE)
    return _mm_setr_epi16(c0, c1, c2, c3, c4, c5, c6, c7);
  #else
    simde_v128_private r_;

    r_.i16[0] = c0;
    r_.i16[1] = c1;
    r_.i16[2] = c2;
    r_.i16[3] = c3;
    r_.i16[4] = c4;
    r_.i16[5] = c5;
    r_.i16[6] = c6;
    r_.i16[7] = c7;

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define \
    wasm_i16x8_make(c0, c1,  c2,  c3,  c4,  c5,  c6,  c7) \
      simde_wasm_i16x8_make((c0), (c1),  (c2),  (c3),  (c4),  (c5),  (c6),  (c7))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i32x4_make (int32_t c0, int32_t c1, int32_t  c2, int32_t  c3) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i32x4_make(c0, c1, c2, c3);
  #elif defined(SIMDE_X86_SSE2_NATIVE)
    return _mm_setr_epi32(c0, c1, c2, c3);
  #else
    simde_v128_private r_;

    r_.i32[0] = c0;
    r_.i32[1] = c1;
    r_.i32[2] = c2;
    r_.i32[3] = c3;

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_make(c0, c1,  c2,  c3) simde_wasm_i32x4_make((c0), (c1),  (c2),  (c3))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i64x2_make (int64_t c0, int64_t c1) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i64x2_make(c0, c1);
  #elif defined(SIMDE_X86_SSE2_NATIVE)
    return _mm_set_epi64x(c1, c0);
  #else
    simde_v128_private r_;

    r_.i64[ 0] = c0;
    r_.i64[ 1] = c1;

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i64x2_make(c0, c1) simde_wasm_i64x2_make((c0), (c1))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f32x4_make (simde_float32 c0, simde_float32 c1, simde_float32  c2, simde_float32  c3) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f32x4_make(c0, c1, c2, c3);
  #else
    simde_v128_private r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128 = _mm_setr_ps(c0, c1, c2, c3);
    #else
      r_.f32[0] = c0;
      r_.f32[1] = c1;
      r_.f32[2] = c2;
      r_.f32[3] = c3;
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f32x4_make(c0, c1,  c2,  c3) simde_wasm_f32x4_make((c0), (c1),  (c2),  (c3))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f64x2_make (simde_float64 c0, simde_float64 c1) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f64x2_make(c0, c1);
  #else
    simde_v128_private r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128d = _mm_set_pd(c1, c0);
    #else
      r_.f64[ 0] = c0;
      r_.f64[ 1] = c1;
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f64x2_make(c0, c1) simde_wasm_f64x2_make((c0), (c1))
#endif

/* const */

#if defined(SIMDE_WASM_SIMD128_NATIVE)
  #define \
    simde_wasm_i8x16_const( \
        c0, c1,  c2,  c3,  c4,  c5,  c6,  c7, \
        c8, c9, c10, c11, c12, c13, c14, c15) \
      wasm_i8x16_const( \
          (c0), (c1),  (c2),  (c3),  (c4),  (c5),  (c6),  (c7), \
          (c8), (c9), (c10), (c11), (c12), (c13), (c14), (c15))
#elif defined(SIMDE_STATEMENT_EXPR_) && defined(SIMDE_ASSERT_CONSTANT_) && defined(SIMDE_STATIC_ASSERT)
  #define \
    simde_wasm_i8x16_const( \
        c0, c1,  c2,  c3,  c4,  c5,  c6,  c7, \
        c8, c9, c10, c11, c12, c13, c14, c15) \
    SIMDE_STATEMENT_EXPR_(({ \
      SIMDE_ASSERT_CONSTANT_(c0); \
      SIMDE_ASSERT_CONSTANT_(c1); \
      SIMDE_ASSERT_CONSTANT_(c2); \
      SIMDE_ASSERT_CONSTANT_(c3); \
      SIMDE_ASSERT_CONSTANT_(c4); \
      SIMDE_ASSERT_CONSTANT_(c5); \
      SIMDE_ASSERT_CONSTANT_(c6); \
      SIMDE_ASSERT_CONSTANT_(c7); \
      SIMDE_ASSERT_CONSTANT_(c8); \
      SIMDE_ASSERT_CONSTANT_(c9); \
      SIMDE_ASSERT_CONSTANT_(c10); \
      SIMDE_ASSERT_CONSTANT_(c11); \
      SIMDE_ASSERT_CONSTANT_(c12); \
      SIMDE_ASSERT_CONSTANT_(c13); \
      SIMDE_ASSERT_CONSTANT_(c13); \
      SIMDE_ASSERT_CONSTANT_(c15); \
      \
      simde_wasm_i8x16_make( \
        c0, c1,  c2,  c3,  c4,  c5,  c6,  c7, \
        c8, c9, c10, c11, c12, c13, c14, c15); \
    }))
#else
  SIMDE_FUNCTION_ATTRIBUTES
  simde_v128_t
  simde_wasm_i8x16_const (
      int8_t c0, int8_t c1, int8_t  c2, int8_t  c3, int8_t  c4, int8_t  c5, int8_t  c6, int8_t  c7,
      int8_t c8, int8_t c9, int8_t c10, int8_t c11, int8_t c12, int8_t c13, int8_t c14, int8_t c15) {
    return simde_wasm_i8x16_make(
        c0, c1,  c2,  c3,  c4,  c5,  c6,  c7,
        c8, c9, c10, c11, c12, c13, c14, c15);
  }
#endif
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define \
    wasm_i8x16_const( \
        c0, c1,  c2,  c3,  c4,  c5,  c6,  c7, \
        c8, c9, c10, c11, c12, c13, c14, c15) \
      simde_wasm_i8x16_const( \
          (c0), (c1),  (c2),  (c3),  (c4),  (c5),  (c6),  (c7), \
          (c8), (c9), (c10), (c11), (c12), (c13), (c14), (c15))
#endif

#if defined(SIMDE_WASM_SIMD128_NATIVE)
  #define \
    simde_wasm_i16x8_const( \
        c0, c1,  c2,  c3,  c4,  c5,  c6,  c7) \
      wasm_i16x8_const( \
          (c0), (c1),  (c2),  (c3),  (c4),  (c5),  (c6),  (c7))
#elif defined(SIMDE_STATEMENT_EXPR_) && defined(SIMDE_ASSERT_CONSTANT_) && defined(SIMDE_STATIC_ASSERT)
  #define \
    simde_wasm_i16x8_const( \
        c0, c1,  c2,  c3,  c4,  c5,  c6,  c7) \
    SIMDE_STATEMENT_EXPR_(({ \
      SIMDE_ASSERT_CONSTANT_(c0); \
      SIMDE_ASSERT_CONSTANT_(c1); \
      SIMDE_ASSERT_CONSTANT_(c2); \
      SIMDE_ASSERT_CONSTANT_(c3); \
      SIMDE_ASSERT_CONSTANT_(c4); \
      SIMDE_ASSERT_CONSTANT_(c5); \
      SIMDE_ASSERT_CONSTANT_(c6); \
      SIMDE_ASSERT_CONSTANT_(c7); \
      \
      simde_wasm_i16x8_make( \
        c0, c1,  c2,  c3,  c4,  c5,  c6,  c7); \
    }))
#else
  SIMDE_FUNCTION_ATTRIBUTES
  simde_v128_t
  simde_wasm_i16x8_const (
      int16_t c0, int16_t c1, int16_t  c2, int16_t  c3, int16_t  c4, int16_t  c5, int16_t  c6, int16_t  c7) {
    return simde_wasm_i16x8_make(
        c0, c1,  c2,  c3,  c4,  c5,  c6,  c7);
  }
#endif
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define \
    wasm_i16x8_const( \
        c0, c1, c2, c3, c4, c5, c6, c7) \
      simde_wasm_i16x8_const( \
          (c0), (c1), (c2), (c3), (c4), (c5), (c6), (c7))
#endif

#if defined(SIMDE_WASM_SIMD128_NATIVE)
  #define \
    simde_wasm_i32x4_const( \
        c0, c1,  c2,  c3) \
      wasm_i32x4_const( \
          (c0), (c1),  (c2),  (c3))
#elif defined(SIMDE_STATEMENT_EXPR_) && defined(SIMDE_ASSERT_CONSTANT_) && defined(SIMDE_STATIC_ASSERT)
  #define \
    simde_wasm_i32x4_const( \
        c0, c1,  c2,  c3) \
    SIMDE_STATEMENT_EXPR_(({ \
      SIMDE_ASSERT_CONSTANT_(c0); \
      SIMDE_ASSERT_CONSTANT_(c1); \
      SIMDE_ASSERT_CONSTANT_(c2); \
      SIMDE_ASSERT_CONSTANT_(c3); \
      \
      simde_wasm_i32x4_make( \
        c0, c1,  c2,  c3); \
    }))
#else
  SIMDE_FUNCTION_ATTRIBUTES
  simde_v128_t
  simde_wasm_i32x4_const (
      int32_t c0, int32_t c1, int32_t  c2, int32_t  c3) {
    return simde_wasm_i32x4_make(
        c0, c1,  c2,  c3);
  }
#endif
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define \
    wasm_i32x4_const( \
        c0, c1,  c2,  c3) \
      simde_wasm_i32x4_const( \
          (c0), (c1),  (c2),  (c3))
#endif

#if defined(SIMDE_WASM_SIMD128_NATIVE)
  #define \
    simde_wasm_i64x2_const( \
        c0, c1) \
      wasm_i64x2_const( \
          (c0), (c1))
#elif defined(SIMDE_STATEMENT_EXPR_) && defined(SIMDE_ASSERT_CONSTANT_) && defined(SIMDE_STATIC_ASSERT)
  #define \
    simde_wasm_i64x2_const( \
        c0, c1) \
    SIMDE_STATEMENT_EXPR_(({ \
      SIMDE_ASSERT_CONSTANT_(c0); \
      SIMDE_ASSERT_CONSTANT_(c1); \
      \
      simde_wasm_i64x2_make( \
        c0, c1); \
    }))
#else
  SIMDE_FUNCTION_ATTRIBUTES
  simde_v128_t
  simde_wasm_i64x2_const (
      int64_t c0, int64_t c1) {
    return simde_wasm_i64x2_make(
        c0, c1);
  }
#endif
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define \
    wasm_i64x2_const( \
        c0, c1) \
      simde_wasm_i64x2_const( \
          (c0), (c1))
#endif

/* splat */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i8x16_splat (int8_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i8x16_splat(a);
  #else
    simde_v128_private r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_set1_epi8(a);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_i8 = vdupq_n_s8(a);
    #elif defined(SIMDE_POWER_ALTIVEC_P6_NATIVE) || defined(SIMDE_ZARCH_ZVECTOR_13_NATIVE)
      r_.altivec_i8 = vec_splats(a);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i8) / sizeof(r_.i8[0])) ; i++) {
        r_.i8[i] = a;
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i8x16_splat(a) simde_wasm_i8x16_splat((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_splat (int16_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_splat(a);
  #else
    simde_v128_private r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_set1_epi16(a);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_i16 = vdupq_n_s16(a);
    #elif defined(SIMDE_POWER_ALTIVEC_P6_NATIVE) || defined(SIMDE_ZARCH_ZVECTOR_13_NATIVE)
      r_.altivec_i16 = vec_splats(a);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i16) / sizeof(r_.i16[0])) ; i++) {
        r_.i16[i] = a;
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_splat(a) simde_wasm_i16x8_splat((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i32x4_splat (int32_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i32x4_splat(a);
  #else
    simde_v128_private r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_set1_epi32(a);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_i32 = vdupq_n_s32(a);
    #elif defined(SIMDE_POWER_ALTIVEC_P6_NATIVE) || defined(SIMDE_ZARCH_ZVECTOR_13_NATIVE)
      r_.altivec_i32 = vec_splats(a);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32) / sizeof(r_.i32[0])) ; i++) {
        r_.i32[i] = a;
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_splat(a) simde_wasm_i32x4_splat((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i64x2_splat (int64_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i64x2_splat(a);
  #else
    simde_v128_private r_;

    #if defined(SIMDE_X86_SSE2_NATIVE) && (!defined(HEDLEY_MSVC_VERSION) || HEDLEY_MSVC_VERSION_CHECK(19,0,0))
      r_.sse_m128i = _mm_set1_epi64x(a);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_i64 = vdupq_n_s64(a);
    #elif defined(SIMDE_POWER_ALTIVEC_P7_NATIVE) || defined(SIMDE_ZARCH_ZVECTOR_13_NATIVE)
      r_.altivec_i64 = vec_splats(HEDLEY_STATIC_CAST(signed long long, a));
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i64) / sizeof(r_.i64[0])) ; i++) {
        r_.i64[i] = a;
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i64x2_splat(a) simde_wasm_i64x2_splat((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f32x4_splat (simde_float32 a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f32x4_splat(a);
  #else
    simde_v128_private r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128 = _mm_set1_ps(a);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_f32 = vdupq_n_f32(a);
    #elif defined(SIMDE_POWER_ALTIVEC_P7_NATIVE) || defined(SIMDE_ZARCH_ZVECTOR_13_NATIVE)
      r_.altivec_f32 = vec_splats(a);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f32) / sizeof(r_.f32[0])) ; i++) {
        r_.f32[i] = a;
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f32x4_splat(a) simde_wasm_f32x4_splat((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f64x2_splat (simde_float64 a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f64x2_splat(a);
  #else
    simde_v128_private r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128d = _mm_set1_pd(a);
    #elif defined(SIMDE_ARM_NEON_A64V8_NATIVE)
      r_.neon_f64 = vdupq_n_f64(a);
    #elif defined(SIMDE_POWER_ALTIVEC_P7_NATIVE) || defined(SIMDE_ZARCH_ZVECTOR_13_NATIVE)
      r_.altivec_f64 = vec_splats(a);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f64) / sizeof(r_.f64[0])) ; i++) {
        r_.f64[i] = a;
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f64x2_splat(a) simde_wasm_f64x2_splat((a))
#endif

/* load_splat */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_v8x16_load_splat (const void * mem) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_v8x16_load_splat(mem);
  #else
    int8_t v;
    simde_memcpy(&v, mem, sizeof(v));
    return simde_wasm_i8x16_splat(v);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_v8x16_load_splat(mem) simde_wasm_v8x16_load_splat((mem))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_v16x8_load_splat (const void * mem) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_v16x8_load_splat(mem);
  #else
    int16_t v;
    simde_memcpy(&v, mem, sizeof(v));
    return simde_wasm_i16x8_splat(v);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_v16x8_load_splat(mem) simde_wasm_v16x8_load_splat((mem))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_v32x4_load_splat (const void * mem) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_v32x4_load_splat(mem);
  #else
    int32_t v;
    simde_memcpy(&v, mem, sizeof(v));
    return simde_wasm_i32x4_splat(v);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_v32x4_load_splat(mem) simde_wasm_v32x4_load_splat((mem))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_v64x2_load_splat (const void * mem) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_v64x2_load_splat(mem);
  #else
    int64_t v;
    simde_memcpy(&v, mem, sizeof(v));
    return simde_wasm_i64x2_splat(v);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_v64x2_load_splat(mem) simde_wasm_v64x2_load_splat((mem))
#endif

/* extract_lane
 *
 * Note that, unlike normal WASM SIMD128 we return intN_t instead of
 * int for sizeof(X) <= sizeof(int).  This is done for portability;
 * the regular API doesn't have to worry about things like int being
 * 16 bits (like on AVR).
 *
 * This does mean that code which works in SIMDe may not work without
 * changes on WASM, but luckily the necessary changes (i.e., casting
 * the return values to smaller type when assigning to the smaller
 * type) mean the code will work in *both* SIMDe and a native
 * implementation.  If you use the simde_* prefixed functions it will
 * always work. */

SIMDE_FUNCTION_ATTRIBUTES
int8_t
simde_wasm_i8x16_extract_lane (simde_v128_t a, const int lane) {
  simde_v128_private a_ = simde_v128_to_private(a);
  return a_.i8[lane & 15];
}
#if defined(SIMDE_WASM_SIMD128_NATIVE)
  #define simde_wasm_i8x16_extract_lane(a, lane) HEDLEY_STATIC_CAST(int8_t, wasm_i8x16_extract_lane((a), (lane)))
#elif defined(SIMDE_X86_SSE4_1_NATIVE)
  #define simde_wasm_i8x16_extract_lane(a, lane) HEDLEY_STATIC_CAST(int8_t, _mm_extract_epi8(simde_v128_to_m128i(a), (lane) & 15))
#elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
  #define simde_wasm_i8x16_extract_lane(a, lane) vgetq_lane_s8(simde_v128_to_neon_i8(a), (lane) & 15)
#endif
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i8x16_extract_lane(a, lane) simde_wasm_i8x16_extract_lane((a), (lane))
#endif

SIMDE_FUNCTION_ATTRIBUTES
int16_t
simde_wasm_i16x8_extract_lane (simde_v128_t a, const int lane) {
  simde_v128_private a_ = simde_v128_to_private(a);
  return a_.i16[lane & 7];
}
#if defined(SIMDE_WASM_SIMD128_NATIVE)
  #define simde_wasm_i16x8_extract_lane(a, lane) HEDLEY_STATIC_CAST(int16_t, wasm_i16x8_extract_lane((a), (lane)))
#elif defined(SIMDE_X86_SSE2_NATIVE)
  #define simde_wasm_i16x8_extract_lane(a, lane) HEDLEY_STATIC_CAST(int16_t, _mm_extract_epi16((a), (lane) & 7))
#elif defined(SIMDE_ARM_NEON_A32V7_NATIVE) && !defined(SIMDE_BUG_CLANG_BAD_VGET_SET_LANE_TYPES)
  #define simde_wasm_i16x8_extract_lane(a, lane) vgetq_lane_s16(simde_v128_to_neon_i16(a), (lane) & 7)
#endif
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_extract_lane(a, lane) simde_wasm_i16x8_extract_lane((a), (lane))
#endif

SIMDE_FUNCTION_ATTRIBUTES
int32_t
simde_wasm_i32x4_extract_lane (simde_v128_t a, const int lane) {
  simde_v128_private a_ = simde_v128_to_private(a);
  return a_.i32[lane & 3];
}
#if defined(SIMDE_WASM_SIMD128_NATIVE)
  #define simde_wasm_i32x4_extract_lane(a, lane) HEDLEY_STATIC_CAST(int32_t, wasm_i32x4_extract_lane((a), (lane)))
#elif defined(SIMDE_X86_SSE4_1_NATIVE)
  #define simde_wasm_i32x4_extract_lane(a, lane) HEDLEY_STATIC_CAST(int32_t, _mm_extract_epi32((a), (lane) & 3))
#elif defined(SIMDE_ARM_NEON_A32V7_NATIVE) && !defined(SIMDE_BUG_CLANG_BAD_VGET_SET_LANE_TYPES)
  #define simde_wasm_i32x4_extract_lane(a, lane) vgetq_lane_s32(simde_v128_to_neon_i32(a), (lane) & 3)
#endif
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_extract_lane(a, lane) simde_wasm_i32x4_extract_lane((a), (lane))
#endif

SIMDE_FUNCTION_ATTRIBUTES
int64_t
simde_wasm_i64x2_extract_lane (simde_v128_t a, const int lane) {
  simde_v128_private a_ = simde_v128_to_private(a);
  return a_.i64[lane & 1];
}
#if defined(SIMDE_WASM_SIMD128_NATIVE)
  #define simde_wasm_i64x2_extract_lane(a, lane) HEDLEY_STATIC_CAST(int64_t, wasm_i64x2_extract_lane((a), (lane)))
#elif defined(SIMDE_X86_SSE4_1_NATIVE) && defined(SIMDE_ARCH_AMD64)
  #define simde_wasm_i64x2_extract_lane(a, lane) HEDLEY_STATIC_CAST(int64_t, _mm_extract_epi64((a), (lane) & 1))
#elif defined(SIMDE_ARM_NEON_A32V7_NATIVE) && !defined(SIMDE_BUG_CLANG_BAD_VGET_SET_LANE_TYPES)
  #define simde_wasm_i64x2_extract_lane(a, lane) vgetq_lane_s64(simde_v128_to_neon_i64(a), (lane) & 1)
#endif
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i64x2_extract_lane(a, lane) simde_wasm_i64x2_extract_lane((a), (lane))
#endif

#if defined(SIMDE_WASM_SIMD128_NATIVE)
  #define simde_wasm_u8x16_extract_lane(a, lane) HEDLEY_STATIC_CAST(uint8_t, wasm_u8x16_extract_lane((a), (lane)))
#elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
  #define simde_wasm_u8x16_extract_lane(a, lane) vgetq_lane_u8(simde_v128_to_neon_u8(a), (lane) & 15)
#else
  #define simde_wasm_u8x16_extract_lane(a, lane) HEDLEY_STATIC_CAST(uint8_t, simde_wasm_i8x16_extract_lane((a), (lane)))
#endif
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u8x16_extract_lane(a, lane) simde_wasm_u8x16_extract_lane((a), (lane))
#endif

#if defined(SIMDE_WASM_SIMD128_NATIVE)
  #define simde_wasm_u16x8_extract_lane(a, lane) HEDLEY_STATIC_CAST(uint16_t, wasm_u16x8_extract_lane((a), (lane)))
#elif defined(SIMDE_ARM_NEON_A32V7_NATIVE) && !defined(SIMDE_BUG_CLANG_BAD_VGET_SET_LANE_TYPES)
  #define simde_wasm_u16x8_extract_lane(a, lane) vgetq_lane_u16(simde_v128_to_neon_u16(a), (lane) & 7)
#else
  #define simde_wasm_u16x8_extract_lane(a, lane) HEDLEY_STATIC_CAST(uint16_t, simde_wasm_i16x8_extract_lane((a), (lane)))
#endif
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u16x8_extract_lane(a, lane) simde_wasm_u16x8_extract_lane((a), (lane))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_float32
simde_wasm_f32x4_extract_lane (simde_v128_t a, const int lane) {
  simde_v128_private a_ = simde_v128_to_private(a);
  return a_.f32[lane & 3];
}
#if defined(SIMDE_WASM_SIMD128_NATIVE)
  #define simde_wasm_f32x4_extract_lane(a, lane) wasm_f32x4_extract_lane((a), (lane))
#elif defined(SIMDE_X86_SSE4_1_NATIVE)
  #define simde_wasm_f32x4(a, lane) _mm_extract_ps(simde_v128_to_m128(a), (lane) & 3)
#elif defined(SIMDE_ARM_NEON_A32V7_NATIVE) && !defined(SIMDE_BUG_CLANG_BAD_VGET_SET_LANE_TYPES)
  #define simde_wasm_f32x4_extract_lane(a, lane) vgetq_lane_f32(simde_v128_to_neon_f32(a), (lane) & 3)
#endif
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f32x4_extract_lane(a, lane) simde_wasm_f32x4_extract_lane((a), (lane))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_float64
simde_wasm_f64x2_extract_lane (simde_v128_t a, const int lane) {
  simde_v128_private a_ = simde_v128_to_private(a);
  return a_.f64[lane & 1];
}
#if defined(SIMDE_WASM_SIMD128_NATIVE)
  #define simde_wasm_f64x2_extract_lane(a, lane) wasm_f64x2_extract_lane((a), (lane))
#elif defined(SIMDE_ARM_NEON_A64V8_NATIVE) && !defined(SIMDE_BUG_CLANG_BAD_VGET_SET_LANE_TYPES)
  #define simde_wasm_f64x2_extract_lane(a, lane) vgetq_lane_f64(simde_v128_to_neon_f64(a), (lane) & 1)
#endif
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f64x2_extract_lane(a, lane) simde_wasm_f64x2_extract_lane((a), (lane))
#endif

/* replace_lane */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i8x16_replace_lane (simde_v128_t a, const int lane, int8_t value) {
  simde_v128_private a_ = simde_v128_to_private(a);
  a_.i8[lane & 15] = value;
  return simde_v128_from_private(a_);
}
#if defined(SIMDE_WASM_SIMD128_NATIVE)
  #define simde_wasm_i8x16_replace_lane(a, lane, value) wasm_i8x16_replace_lane((a), (lane), (value))
#elif defined(SIMDE_X86_SSE4_1_NATIVE)
  #define simde_wasm_i8x16_replace_lane(a, lane, value) _mm_insert_epi8((a), (value), (lane) & 15)
#elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
  #define simde_wasm_i8x16_replace_lane(a, lane, value) simde_v128_from_neon_i8(vsetq_lane_s8((value), simde_v128_to_neon_i8(a), (lane) & 15))
#endif
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i8x16_replace_lane(a, lane, value) simde_wasm_i8x16_replace_lane((a), (lane), (value))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_replace_lane (simde_v128_t a, const int lane, int16_t value) {
  simde_v128_private a_ = simde_v128_to_private(a);
  a_.i16[lane & 7] = value;
  return simde_v128_from_private(a_);
}
#if defined(SIMDE_WASM_SIMD128_NATIVE)
  #define simde_wasm_i16x8_replace_lane(a, lane, value) wasm_i16x8_replace_lane((a), (lane), (value))
#elif defined(SIMDE_X86_SSE2_NATIVE)
  #define simde_wasm_i16x8_replace_lane(a, lane, value) _mm_insert_epi16((a), (value), (lane) & 7)
#elif defined(SIMDE_ARM_NEON_A32V7_NATIVE) && !defined(SIMDE_BUG_CLANG_BAD_VGET_SET_LANE_TYPES)
  #define simde_wasm_i16x8_replace_lane(a, lane, value) simde_v128_from_neon_i16(vsetq_lane_s16((value), simde_v128_to_neon_i16(a), (lane) & 7))
#endif
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_replace_lane(a, lane, value) simde_wasm_i16x8_replace_lane((a), (lane), (value))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i32x4_replace_lane (simde_v128_t a, const int lane, int32_t value) {
  simde_v128_private a_ = simde_v128_to_private(a);
  a_.i32[lane & 3] = value;
  return simde_v128_from_private(a_);
}
#if defined(SIMDE_WASM_SIMD128_NATIVE)
  #define simde_wasm_i32x4_replace_lane(a, lane, value) wasm_i32x4_replace_lane((a), (lane), (value))
#elif defined(SIMDE_X86_SSE4_1_NATIVE)
  #define simde_wasm_i32x4_replace_lane(a, lane, value) _mm_insert_epi32((a), (value), (lane) & 3)
#elif defined(SIMDE_ARM_NEON_A32V7_NATIVE) && !defined(SIMDE_BUG_CLANG_BAD_VGET_SET_LANE_TYPES)
  #define simde_wasm_i32x4_replace_lane(a, lane, value) simde_v128_from_neon_i32(vsetq_lane_s32((value), simde_v128_to_neon_i32(a), (lane) & 3))
#endif
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_replace_lane(a, lane, value) simde_wasm_i32x4_replace_lane((a), (lane), (value))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i64x2_replace_lane (simde_v128_t a, const int lane, int64_t value) {
  simde_v128_private a_ = simde_v128_to_private(a);
  a_.i64[lane & 1] = value;
  return simde_v128_from_private(a_);
}
#if defined(SIMDE_WASM_SIMD128_NATIVE)
  #define simde_wasm_i64x2_replace_lane(a, lane, value) wasm_i64x2_replace_lane((a), (lane), (value))
#elif defined(SIMDE_X86_SSE4_1_NATIVE) && defined(SIMDE_ARCH_AMD64)
  #define simde_wasm_i64x2_replace_lane(a, lane, value) _mm_insert_epi64((a), (value), (lane) & 1)
#elif defined(SIMDE_ARM_NEON_A32V7_NATIVE) && !defined(SIMDE_BUG_CLANG_BAD_VGET_SET_LANE_TYPES)
  #define simde_wasm_i64x2_replace_lane(a, lane, value) simde_v128_from_neon_i64(vsetq_lane_s64((value), simde_v128_to_neon_i64(a), (lane) & 1))
#endif
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i64x2_replace_lane(a, lane, value) simde_wasm_i64x2_replace_lane((a), (lane), (value))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f32x4_replace_lane (simde_v128_t a, const int lane, simde_float32 value) {
  simde_v128_private a_ = simde_v128_to_private(a);
  a_.f32[lane & 3] = value;
  return simde_v128_from_private(a_);
}
#if defined(SIMDE_WASM_SIMD128_NATIVE)
  #define simde_wasm_f32x4_replace_lane(a, lane, value) wasm_f32x4_replace_lane((a), (lane), (value))
#elif defined(SIMDE_ARM_NEON_A32V7_NATIVE) && !defined(SIMDE_BUG_CLANG_BAD_VGET_SET_LANE_TYPES)
  #define simde_wasm_f32x4_replace_lane(a, lane, value) simde_v128_from_neon_f32(vsetq_lane_f32((value), simde_v128_to_neon_f32(a), (lane) & 3))
#endif
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f32x4_replace_lane(a, lane, value) simde_wasm_f32x4_replace_lane((a), (lane), (value))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f64x2_replace_lane (simde_v128_t a, const int lane, simde_float64 value) {
  simde_v128_private a_ = simde_v128_to_private(a);
  a_.f64[lane & 1] = value;
  return simde_v128_from_private(a_);
}
#if defined(SIMDE_WASM_SIMD128_NATIVE)
  #define simde_wasm_f64x2_replace_lane(a, lane, value) wasm_f64x2_replace_lane((a), (lane), (value))
#elif defined(SIMDE_ARM_NEON_A64V8_NATIVE) && !defined(SIMDE_BUG_CLANG_BAD_VGET_SET_LANE_TYPES)
  #define simde_wasm_f64x2_replace_lane(a, lane, value) simde_v128_from_neon_f64(vsetq_lane_f64((value), simde_v128_to_neon_f64(a), (lane) & 1))
#endif
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f64x2_replace_lane(a, lane, value) simde_wasm_f64x2_replace_lane((a), (lane), (value))
#endif

/* eq */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i8x16_eq (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i8x16_eq(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_cmpeq_epi8(a_.sse_m128i, b_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u8 = vceqq_s8(a_.neon_i8, b_.neon_i8);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i8 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i8), a_.i8 == b_.i8);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i8) / sizeof(r_.i8[0])) ; i++) {
        r_.i8[i] = (a_.i8[i] == b_.i8[i]) ? ~INT8_C(0) : INT8_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i8x16_eq(a, b) simde_wasm_i8x16_eq((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_eq (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_eq(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_cmpeq_epi16(a_.sse_m128i, b_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u16 = vceqq_s16(a_.neon_i16, b_.neon_i16);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i16 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i16), a_.i16 == b_.i16);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i16) / sizeof(r_.i16[0])) ; i++) {
        r_.i16[i] = (a_.i16[i] == b_.i16[i]) ? ~INT16_C(0) : INT16_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_eq(a, b) simde_wasm_i16x8_eq((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i32x4_eq (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i32x4_eq(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_cmpeq_epi32(a_.sse_m128i, b_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u32 = vceqq_s32(a_.neon_i32, b_.neon_i32);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i32), a_.i32 == b_.i32);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32) / sizeof(r_.i32[0])) ; i++) {
        r_.i32[i] = (a_.i32[i] == b_.i32[i]) ? ~INT32_C(0) : INT32_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_eq(a, b) simde_wasm_i32x4_eq((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f32x4_eq (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f32x4_eq(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128 = _mm_cmpeq_ps(a_.sse_m128, b_.sse_m128);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u32 = vceqq_f32(a_.neon_f32, b_.neon_f32);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i32), a_.f32 == b_.f32);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f32) / sizeof(r_.f32[0])) ; i++) {
        r_.i32[i] = (a_.f32[i] == b_.f32[i]) ? ~INT32_C(0) : INT32_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f32x4_eq(a, b) simde_wasm_f32x4_eq((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f64x2_eq (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f64x2_eq(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128d = _mm_cmpeq_pd(a_.sse_m128d, b_.sse_m128d);
    #elif defined(SIMDE_ARM_NEON_A64V8_NATIVE)
      r_.neon_u64 = vceqq_f64(a_.neon_f64, b_.neon_f64);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i32), a_.f64 == b_.f64);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f64) / sizeof(r_.f64[0])) ; i++) {
        r_.i64[i] = (a_.f64[i] == b_.f64[i]) ? ~INT64_C(0) : INT64_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f64x2_eq(a, b) simde_wasm_f64x2_eq((a), (b))
#endif

/* ne */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i8x16_ne (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i8x16_ne(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u8 = vmvnq_u8(vceqq_s8(a_.neon_i8, b_.neon_i8));
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i8 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i8), a_.i8 != b_.i8);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i8) / sizeof(r_.i8[0])) ; i++) {
        r_.i8[i] = (a_.i8[i] != b_.i8[i]) ? ~INT8_C(0) : INT8_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i8x16_ne(a, b) simde_wasm_i8x16_ne((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_ne (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_ne(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u16 = vmvnq_u16(vceqq_s16(a_.neon_i16, b_.neon_i16));
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i16 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i16), a_.i16 != b_.i16);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i16) / sizeof(r_.i16[0])) ; i++) {
        r_.i16[i] = (a_.i16[i] != b_.i16[i]) ? ~INT16_C(0) : INT16_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_ne(a, b) simde_wasm_i16x8_ne((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i32x4_ne (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i32x4_ne(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u32 = vmvnq_u32(vceqq_s32(a_.neon_i32, b_.neon_i32));
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i32), a_.i32 != b_.i32);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32) / sizeof(r_.i32[0])) ; i++) {
        r_.i32[i] = (a_.i32[i] != b_.i32[i]) ? ~INT32_C(0) : INT32_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_ne(a, b) simde_wasm_i32x4_ne((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f32x4_ne (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f32x4_ne(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128 = _mm_cmpneq_ps(a_.sse_m128, b_.sse_m128);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u32 = vmvnq_u32(vceqq_f32(a_.neon_f32, b_.neon_f32));
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i32), a_.f32 != b_.f32);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f32) / sizeof(r_.f32[0])) ; i++) {
        r_.i32[i] = (a_.f32[i] != b_.f32[i]) ? ~INT32_C(0) : INT32_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f32x4_ne(a, b) simde_wasm_f32x4_ne((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f64x2_ne (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f64x2_ne(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128d = _mm_cmpneq_pd(a_.sse_m128d, b_.sse_m128d);
    #elif defined(SIMDE_ARM_NEON_A64V8_NATIVE)
      r_.neon_u32 = vmvnq_u32(vreinterpretq_u32_u64(vceqq_f64(a_.neon_f64, b_.neon_f64)));
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i32), a_.f64 != b_.f64);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f64) / sizeof(r_.f64[0])) ; i++) {
        r_.i64[i] = (a_.f64[i] != b_.f64[i]) ? ~INT64_C(0) : INT64_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f64x2_ne(a, b) simde_wasm_f64x2_ne((a), (b))
#endif

/* lt */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i8x16_lt (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i8x16_lt(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_cmplt_epi8(a_.sse_m128i, b_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u8 = vcltq_s8(a_.neon_i8, b_.neon_i8);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i8 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i8), a_.i8 < b_.i8);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i8) / sizeof(r_.i8[0])) ; i++) {
        r_.i8[i] = (a_.i8[i] < b_.i8[i]) ? ~INT8_C(0) : INT8_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i8x16_lt(a, b) simde_wasm_i8x16_lt((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_lt (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_lt(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_cmplt_epi16(a_.sse_m128i, b_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u16 = vcltq_s16(a_.neon_i16, b_.neon_i16);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i16 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i16), a_.i16 < b_.i16);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i16) / sizeof(r_.i16[0])) ; i++) {
        r_.i16[i] = (a_.i16[i] < b_.i16[i]) ? ~INT16_C(0) : INT16_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_lt(a, b) simde_wasm_i16x8_lt((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i32x4_lt (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i32x4_lt(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_cmplt_epi32(a_.sse_m128i, b_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u32 = vcltq_s32(a_.neon_i32, b_.neon_i32);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i32), a_.i32 < b_.i32);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32) / sizeof(r_.i32[0])) ; i++) {
        r_.i32[i] = (a_.i32[i] < b_.i32[i]) ? ~INT32_C(0) : INT32_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_lt(a, b) simde_wasm_i32x4_lt((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u8x16_lt (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u8x16_lt(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u8 = vcltq_u8(a_.neon_u8, b_.neon_u8);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.u8 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.u8), a_.u8 < b_.u8);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u8) / sizeof(r_.u8[0])) ; i++) {
        r_.u8[i] = (a_.u8[i] < b_.u8[i]) ? ~UINT8_C(0) : UINT8_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u8x16_lt(a, b) simde_wasm_u8x16_lt((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u16x8_lt (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u16x8_lt(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u16 = vcltq_u16(a_.neon_u16, b_.neon_u16);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.u16 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.u16), a_.u16 < b_.u16);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u16) / sizeof(r_.u16[0])) ; i++) {
        r_.u16[i] = (a_.u16[i] < b_.u16[i]) ? ~UINT16_C(0) : UINT16_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u16x8_lt(a, b) simde_wasm_u16x8_lt((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u32x4_lt (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u32x4_lt(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u32 = vcltq_u32(a_.neon_u32, b_.neon_u32);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.u32 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.u32), a_.u32 < b_.u32);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u32) / sizeof(r_.u32[0])) ; i++) {
        r_.u32[i] = (a_.u32[i] < b_.u32[i]) ? ~UINT32_C(0) : UINT32_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u32x4_lt(a, b) simde_wasm_u32x4_lt((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f32x4_lt (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f32x4_lt(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128 = _mm_cmplt_ps(a_.sse_m128, b_.sse_m128);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u32 = vcltq_f32(a_.neon_f32, b_.neon_f32);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i32), a_.f32 < b_.f32);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f32) / sizeof(r_.f32[0])) ; i++) {
        r_.i32[i] = (a_.f32[i] < b_.f32[i]) ? ~INT32_C(0) : INT32_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f32x4_lt(a, b) simde_wasm_f32x4_lt((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f64x2_lt (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f64x2_lt(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128d = _mm_cmplt_pd(a_.sse_m128d, b_.sse_m128d);
    #elif defined(SIMDE_ARM_NEON_A64V8_NATIVE)
      r_.neon_u64 = vcltq_f64(a_.neon_f64, b_.neon_f64);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i32), a_.f64 < b_.f64);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f64) / sizeof(r_.f64[0])) ; i++) {
        r_.i64[i] = (a_.f64[i] < b_.f64[i]) ? ~INT64_C(0) : INT64_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f64x2_lt(a, b) simde_wasm_f64x2_lt((a), (b))
#endif

/* gt */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i8x16_gt (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i8x16_gt(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_cmpgt_epi8(a_.sse_m128i, b_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u8 = vcgtq_s8(a_.neon_i8, b_.neon_i8);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i8 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i8), a_.i8 > b_.i8);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i8) / sizeof(r_.i8[0])) ; i++) {
        r_.i8[i] = (a_.i8[i] > b_.i8[i]) ? ~INT8_C(0) : INT8_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i8x16_gt(a, b) simde_wasm_i8x16_gt((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_gt (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_gt(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_cmpgt_epi16(a_.sse_m128i, b_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u16 = vcgtq_s16(a_.neon_i16, b_.neon_i16);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i16 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i16), a_.i16 > b_.i16);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i16) / sizeof(r_.i16[0])) ; i++) {
        r_.i16[i] = (a_.i16[i] > b_.i16[i]) ? ~INT16_C(0) : INT16_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_gt(a, b) simde_wasm_i16x8_gt((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i32x4_gt (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i32x4_gt(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_cmpgt_epi32(a_.sse_m128i, b_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u32 = vcgtq_s32(a_.neon_i32, b_.neon_i32);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i32), a_.i32 > b_.i32);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32) / sizeof(r_.i32[0])) ; i++) {
        r_.i32[i] = (a_.i32[i] > b_.i32[i]) ? ~INT32_C(0) : INT32_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_gt(a, b) simde_wasm_i32x4_gt((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u8x16_gt (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u8x16_gt(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u8 = vcgtq_u8(a_.neon_u8, b_.neon_u8);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.u8 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.u8), a_.u8 > b_.u8);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u8) / sizeof(r_.u8[0])) ; i++) {
        r_.u8[i] = (a_.u8[i] > b_.u8[i]) ? ~UINT8_C(0) : UINT8_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u8x16_gt(a, b) simde_wasm_u8x16_gt((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u16x8_gt (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u16x8_gt(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u16 = vcgtq_u16(a_.neon_u16, b_.neon_u16);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.u16 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.u16), a_.u16 > b_.u16);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u16) / sizeof(r_.u16[0])) ; i++) {
        r_.u16[i] = (a_.u16[i] > b_.u16[i]) ? ~UINT16_C(0) : UINT16_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u16x8_gt(a, b) simde_wasm_u16x8_gt((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u32x4_gt (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u32x4_gt(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u32 = vcgtq_u32(a_.neon_u32, b_.neon_u32);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.u32 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.u32), a_.u32 > b_.u32);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u32) / sizeof(r_.u32[0])) ; i++) {
        r_.u32[i] = (a_.u32[i] > b_.u32[i]) ? ~UINT32_C(0) : UINT32_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u32x4_gt(a, b) simde_wasm_u32x4_gt((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f32x4_gt (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f32x4_gt(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128 = _mm_cmpgt_ps(a_.sse_m128, b_.sse_m128);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u32 = vcgtq_f32(a_.neon_f32, b_.neon_f32);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i32), a_.f32 > b_.f32);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f32) / sizeof(r_.f32[0])) ; i++) {
        r_.i32[i] = (a_.f32[i] > b_.f32[i]) ? ~INT32_C(0) : INT32_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f32x4_gt(a, b) simde_wasm_f32x4_gt((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f64x2_gt (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f64x2_gt(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128d = _mm_cmpgt_pd(a_.sse_m128d, b_.sse_m128d);
    #elif defined(SIMDE_ARM_NEON_A64V8_NATIVE)
      r_.neon_u64 = vcgtq_f64(a_.neon_f64, b_.neon_f64);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i64 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i64), a_.f64 > b_.f64);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f64) / sizeof(r_.f64[0])) ; i++) {
        r_.i64[i] = (a_.f64[i] > b_.f64[i]) ? ~INT64_C(0) : INT64_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f64x2_gt(a, b) simde_wasm_f64x2_gt((a), (b))
#endif

/* le */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i8x16_le (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i8x16_le(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r_.sse_m128i = _mm_cmpeq_epi8(a_.sse_m128i, _mm_min_epi8(a_.sse_m128i, b_.sse_m128i));
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u8 = vcleq_s8(a_.neon_i8, b_.neon_i8);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i8 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i8), a_.i8 <= b_.i8);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i8) / sizeof(r_.i8[0])) ; i++) {
        r_.i8[i] = (a_.i8[i] <= b_.i8[i]) ? ~INT8_C(0) : INT8_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i8x16_le(a, b) simde_wasm_i8x16_le((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_le (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_le(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r_.sse_m128i = _mm_cmpeq_epi16(a_.sse_m128i, _mm_min_epi16(a_.sse_m128i, b_.sse_m128i));
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u16 = vcleq_s16(a_.neon_i16, b_.neon_i16);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i16 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i16), a_.i16 <= b_.i16);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i16) / sizeof(r_.i16[0])) ; i++) {
        r_.i16[i] = (a_.i16[i] <= b_.i16[i]) ? ~INT16_C(0) : INT16_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_le(a, b) simde_wasm_i16x8_le((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i32x4_le (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i32x4_le(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r_.sse_m128i = _mm_cmpeq_epi32(a_.sse_m128i, _mm_min_epi32(a_.sse_m128i, b_.sse_m128i));
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u32 = vcleq_s32(a_.neon_i32, b_.neon_i32);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i32), a_.i32 <= b_.i32);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32) / sizeof(r_.i32[0])) ; i++) {
        r_.i32[i] = (a_.i32[i] <= b_.i32[i]) ? ~INT32_C(0) : INT32_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_le(a, b) simde_wasm_i32x4_le((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u8x16_le (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u8x16_le(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u8 = vcleq_u8(a_.neon_u8, b_.neon_u8);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.u8 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.u8), a_.u8 <= b_.u8);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u8) / sizeof(r_.u8[0])) ; i++) {
        r_.u8[i] = (a_.u8[i] <= b_.u8[i]) ? ~UINT8_C(0) : UINT8_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u8x16_le(a, b) simde_wasm_u8x16_le((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u16x8_le (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u16x8_le(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u16 = vcleq_u16(a_.neon_u16, b_.neon_u16);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.u16 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.u16), a_.u16 <= b_.u16);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u16) / sizeof(r_.u16[0])) ; i++) {
        r_.u16[i] = (a_.u16[i] <= b_.u16[i]) ? ~UINT16_C(0) : UINT16_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u16x8_le(a, b) simde_wasm_u16x8_le((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u32x4_le (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u32x4_le(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u32 = vcleq_u32(a_.neon_u32, b_.neon_u32);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.u32 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.u32), a_.u32 <= b_.u32);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u32) / sizeof(r_.u32[0])) ; i++) {
        r_.u32[i] = (a_.u32[i] <= b_.u32[i]) ? ~UINT32_C(0) : UINT32_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u32x4_le(a, b) simde_wasm_u32x4_le((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f32x4_le (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f32x4_le(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128 = _mm_cmple_ps(a_.sse_m128, b_.sse_m128);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u32 = vcleq_f32(a_.neon_f32, b_.neon_f32);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i32), a_.f32 <= b_.f32);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f32) / sizeof(r_.f32[0])) ; i++) {
        r_.i32[i] = (a_.f32[i] <= b_.f32[i]) ? ~INT32_C(0) : INT32_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f32x4_le(a, b) simde_wasm_f32x4_le((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f64x2_le (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f64x2_le(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128d = _mm_cmple_pd(a_.sse_m128d, b_.sse_m128d);
    #elif defined(SIMDE_ARM_NEON_A64V8_NATIVE)
      r_.neon_u64 = vcleq_f64(a_.neon_f64, b_.neon_f64);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i32), a_.f64 <= b_.f64);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f64) / sizeof(r_.f64[0])) ; i++) {
        r_.i64[i] = (a_.f64[i] <= b_.f64[i]) ? ~INT64_C(0) : INT64_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f64x2_le(a, b) simde_wasm_f64x2_le((a), (b))
#endif

/* ge */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i8x16_ge (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i8x16_ge(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r_.sse_m128i = _mm_cmpeq_epi8(_mm_min_epi8(a_.sse_m128i, b_.sse_m128i), b_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u8 = vcgeq_s8(a_.neon_i8, b_.neon_i8);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i8 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i8), a_.i8 >= b_.i8);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i8) / sizeof(r_.i8[0])) ; i++) {
        r_.i8[i] = (a_.i8[i] >= b_.i8[i]) ? ~INT8_C(0) : INT8_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i8x16_ge(a, b) simde_wasm_i8x16_ge((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_ge (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_ge(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r_.sse_m128i = _mm_cmpeq_epi16(_mm_min_epi16(a_.sse_m128i, b_.sse_m128i), b_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u16 = vcgeq_s16(a_.neon_i16, b_.neon_i16);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i16 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i16), a_.i16 >= b_.i16);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i16) / sizeof(r_.i16[0])) ; i++) {
        r_.i16[i] = (a_.i16[i] >= b_.i16[i]) ? ~INT16_C(0) : INT16_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_ge(a, b) simde_wasm_i16x8_ge((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i32x4_ge (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i32x4_ge(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r_.sse_m128i = _mm_cmpeq_epi32(_mm_min_epi32(a_.sse_m128i, b_.sse_m128i), b_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u32 = vcgeq_s32(a_.neon_i32, b_.neon_i32);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i32), a_.i32 >= b_.i32);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32) / sizeof(r_.i32[0])) ; i++) {
        r_.i32[i] = (a_.i32[i] >= b_.i32[i]) ? ~INT32_C(0) : INT32_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_ge(a, b) simde_wasm_i32x4_ge((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u8x16_ge (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u8x16_ge(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r_.sse_m128i = _mm_cmpeq_epi8(_mm_min_epu8(a_.sse_m128i, b_.sse_m128i), b_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u8 = vcgeq_u8(a_.neon_u8, b_.neon_u8);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.u8 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.u8), a_.u8 >= b_.u8);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u8) / sizeof(r_.u8[0])) ; i++) {
        r_.u8[i] = (a_.u8[i] >= b_.u8[i]) ? ~UINT8_C(0) : UINT8_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u8x16_ge(a, b) simde_wasm_u8x16_ge((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u16x8_ge (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u16x8_ge(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r_.sse_m128i = _mm_cmpeq_epi16(_mm_min_epu16(a_.sse_m128i, b_.sse_m128i), b_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u16 = vcgeq_u16(a_.neon_u16, b_.neon_u16);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.u16 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.u16), a_.u16 >= b_.u16);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u16) / sizeof(r_.u16[0])) ; i++) {
        r_.u16[i] = (a_.u16[i] >= b_.u16[i]) ? ~UINT16_C(0) : UINT16_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u16x8_ge(a, b) simde_wasm_u16x8_ge((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u32x4_ge (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u32x4_ge(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r_.sse_m128i = _mm_cmpeq_epi32(_mm_min_epu32(a_.sse_m128i, b_.sse_m128i), b_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u32 = vcgeq_u32(a_.neon_u32, b_.neon_u32);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.u32 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.u32), a_.u32 >= b_.u32);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u32) / sizeof(r_.u32[0])) ; i++) {
        r_.u32[i] = (a_.u32[i] >= b_.u32[i]) ? ~UINT32_C(0) : UINT32_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u32x4_ge(a, b) simde_wasm_u32x4_ge((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f32x4_ge (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f32x4_ge(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128 = _mm_cmpge_ps(a_.sse_m128, b_.sse_m128);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_u32 = vcgeq_f32(a_.neon_f32, b_.neon_f32);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i32), a_.f32 >= b_.f32);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f32) / sizeof(r_.f32[0])) ; i++) {
        r_.i32[i] = (a_.f32[i] >= b_.f32[i]) ? ~INT32_C(0) : INT32_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f32x4_ge(a, b) simde_wasm_f32x4_ge((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f64x2_ge (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f64x2_ge(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128d = _mm_cmpge_pd(a_.sse_m128d, b_.sse_m128d);
    #elif defined(SIMDE_ARM_NEON_A64V8_NATIVE)
      r_.neon_u64 = vcgeq_f64(a_.neon_f64, b_.neon_f64);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32 = HEDLEY_REINTERPRET_CAST(__typeof__(r_.i32), a_.f64 >= b_.f64);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f64) / sizeof(r_.f64[0])) ; i++) {
        r_.i64[i] = (a_.f64[i] >= b_.f64[i]) ? ~INT64_C(0) : INT64_C(0);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f64x2_ge(a, b) simde_wasm_f64x2_ge((a), (b))
#endif

/* not */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_v128_not (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_v128_not(a);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_xor_si128(a_.sse_m128i, _mm_set1_epi32(~INT32_C(0)));
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_i32 = vmvnq_s32(a_.neon_i32);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32f = ~a_.i32f;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32f) / sizeof(r_.i32f[0])) ; i++) {
        r_.i32f[i] = ~(a_.i32f[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_v128_not(a) simde_wasm_v128_not((a))
#endif

/* and */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_v128_and (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_v128_and(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_and_si128(a_.sse_m128i, b_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_i32 = vandq_s32(a_.neon_i32, b_.neon_i32);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32f = a_.i32f & b_.i32f;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32f) / sizeof(r_.i32f[0])) ; i++) {
        r_.i32f[i] = a_.i32f[i] & b_.i32f[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_v128_and(a, b) simde_wasm_v128_and((a), (b))
#endif

/* or */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_v128_or (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_v128_or(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_or_si128(a_.sse_m128i, b_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_i32 = vorrq_s32(a_.neon_i32, b_.neon_i32);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32f = a_.i32f | b_.i32f;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32f) / sizeof(r_.i32f[0])) ; i++) {
        r_.i32f[i] = a_.i32f[i] | b_.i32f[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_v128_or(a, b) simde_wasm_v128_or((a), (b))
#endif

/* xor */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_v128_xor (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_v128_xor(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_xor_si128(a_.sse_m128i, b_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_i32 = veorq_s32(a_.neon_i32, b_.neon_i32);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32f = a_.i32f ^ b_.i32f;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32f) / sizeof(r_.i32f[0])) ; i++) {
        r_.i32f[i] = a_.i32f[i] ^ b_.i32f[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_v128_xor(a, b) simde_wasm_v128_xor((a), (b))
#endif

/* andnot */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_v128_andnot (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_v128_andnot(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_andnot_si128(b_.sse_m128i, a_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_i32 = vbicq_s32(a_.neon_i32, b_.neon_i32);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32f = a_.i32f & ~b_.i32f;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32f) / sizeof(r_.i32f[0])) ; i++) {
        r_.i32f[i] = a_.i32f[i] & ~b_.i32f[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_v128_andnot(a, b) simde_wasm_v128_andnot((a), (b))
#endif

/* bitselect */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_v128_bitselect(simde_v128_t a, simde_v128_t b, simde_v128_t mask) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_v128_bitselect(a, b, mask);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      mask_ = simde_v128_to_private(mask),
      r_;

    #if defined(SIMDE_X86_AVX512VL_NATIVE)
      r_.sse_m128i = _mm_ternarylogic_epi32(mask_.sse_m128i, a_.sse_m128i, b_.sse_m128i, 0xca);
    #elif defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i =
        _mm_or_si128(
          _mm_and_si128   (mask_.sse_m128i, a_.sse_m128i),
          _mm_andnot_si128(mask_.sse_m128i, b_.sse_m128i));
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_i32 = vbslq_s32(mask_.neon_u32, a_.neon_i32, b_.neon_i32);
    #elif defined(SIMDE_POWER_ALTIVEC_P6_NATIVE) || defined(SIMDE_ZARCH_ZVECTOR_13_NATIVE)
      r_.altivec_i32 = vec_sel(b_.altivec_i32, a_.altivec_i32, mask_.altivec_u32);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32f = (a_.i32f & mask_.i32f) | (b_.i32f & ~mask_.i32f);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32f) / sizeof(r_.i32f[0])) ; i++) {
        r_.i32f[i] = (a_.i32f[i] & mask_.i32f[i]) | (b_.i32f[i] & ~mask_.i32f[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_v128_bitselect(a, b, c) simde_wasm_v128_bitselect((a), (b), (c))
#endif

/* abs */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i8x16_abs (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i8x16_abs(a);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r_.sse_m128i = _mm_abs_epi8(a_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_i8 = vabsq_s8(a_.neon_i8);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i8) / sizeof(r_.i8[0])) ; i++) {
        r_.i8[i] = (a_.i8[i] < INT8_C(0)) ? -a_.i8[i] : a_.i8[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i8x16_abs(a) simde_wasm_i8x16_abs((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_abs (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_abs(a);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r_.sse_m128i = _mm_abs_epi16(a_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_i16 = vabsq_s16(a_.neon_i16);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i16) / sizeof(r_.i16[0])) ; i++) {
        r_.i16[i] = (a_.i16[i] < INT8_C(0)) ? -a_.i16[i] : a_.i16[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_abs(a) simde_wasm_i16x8_abs((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i32x4_abs (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i32x4_abs(a);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r_.sse_m128i = _mm_abs_epi32(a_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A64V8_NATIVE)
      r_.neon_i32 = vabsq_s32(a_.neon_i32);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32) / sizeof(r_.i32[0])) ; i++) {
        r_.i32[i] = (a_.i32[i] < INT8_C(0)) ? -a_.i32[i] : a_.i32[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_abs(a) simde_wasm_i32x4_abs((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f32x4_abs (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f32x4_abs(a);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_andnot_si128(_mm_set1_epi32(HEDLEY_STATIC_CAST(int32_t, UINT32_C(1) << 31)), a_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_f32 = vabsq_f32(a_.neon_f32);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f32) / sizeof(r_.f32[0])) ; i++) {
        r_.f32[i] = (a_.f32[i] < SIMDE_FLOAT32_C(0.0)) ? -a_.f32[i] : a_.f32[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f32x4_abs(a) simde_wasm_f32x4_abs((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f64x2_abs (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f64x2_abs(a);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_andnot_si128(_mm_set1_epi64x(HEDLEY_STATIC_CAST(int64_t, UINT64_C(1) << 63)), a_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A64V8_NATIVE)
      r_.neon_f64 = vabsq_f64(a_.neon_f64);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f64) / sizeof(r_.f64[0])) ; i++) {
        r_.f64[i] = (a_.f64[i] < SIMDE_FLOAT64_C(0.0)) ? -a_.f64[i] : a_.f64[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f64x2_abs(a) simde_wasm_f64x2_abs((a))
#endif

/* neg */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i8x16_neg (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i8x16_neg(a);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_sub_epi8(_mm_setzero_si128(), a_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_i8 = vnegq_s8(a_.neon_i8);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i8 = -a_.i8;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i8) / sizeof(r_.i8[0])) ; i++) {
        r_.i8[i] = -a_.i8[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i8x16_neg(a) simde_wasm_i8x16_neg((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_neg (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_neg(a);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_sub_epi16(_mm_setzero_si128(), a_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_i16 = vnegq_s16(a_.neon_i16);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i16 = -a_.i16;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i16) / sizeof(r_.i16[0])) ; i++) {
        r_.i16[i] = -a_.i16[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_neg(a) simde_wasm_i16x8_neg((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i32x4_neg (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i32x4_neg(a);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_sub_epi32(_mm_setzero_si128(), a_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_i32 = vnegq_s32(a_.neon_i32);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32 = -a_.i32;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32) / sizeof(r_.i32[0])) ; i++) {
        r_.i32[i] = -a_.i32[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_neg(a) simde_wasm_i32x4_neg((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i64x2_neg (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i64x2_neg(a);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_sub_epi64(_mm_setzero_si128(), a_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A64V8_NATIVE)
      r_.neon_i64 = vnegq_s64(a_.neon_i64);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i64 = -a_.i64;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i64) / sizeof(r_.i64[0])) ; i++) {
        r_.i64[i] = -a_.i64[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i64x2_neg(a) simde_wasm_i64x2_neg((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f32x4_neg (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f32x4_neg(a);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_xor_si128(_mm_set1_epi32(HEDLEY_STATIC_CAST(int32_t, UINT32_C(1) << 31)), a_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A32V7_NATIVE)
      r_.neon_f32 = vnegq_f32(a_.neon_f32);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.f32 = -a_.f32;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f32) / sizeof(r_.f32[0])) ; i++) {
        r_.f32[i] = -a_.f32[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f32x4_neg(a) simde_wasm_f32x4_neg((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f64x2_neg (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f64x2_neg(a);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_xor_si128(_mm_set1_epi64x(HEDLEY_STATIC_CAST(int64_t, UINT64_C(1) << 63)), a_.sse_m128i);
    #elif defined(SIMDE_ARM_NEON_A64V8_NATIVE)
      r_.neon_f64 = vnegq_f64(a_.neon_f64);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.f64 = -a_.f64;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f64) / sizeof(r_.f64[0])) ; i++) {
        r_.f64[i] = -a_.f64[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f64x2_neg(a) simde_wasm_f64x2_neg((a))
#endif

/* any_true */

SIMDE_FUNCTION_ATTRIBUTES
simde_bool
simde_wasm_i8x16_any_true (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i8x16_any_true(a);
  #else
    simde_v128_private a_ = simde_v128_to_private(a);
    int_fast32_t r = 0;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r = !_mm_test_all_zeros(a_.sse_m128i, _mm_set1_epi32(~INT32_C(0)));
    #else
      SIMDE_VECTORIZE_REDUCTION(|:r)
      for (size_t i = 0 ; i < (sizeof(a_.i32f) / sizeof(a_.i32f[0])) ; i++) {
        r |= (a_.i32f[i]);
      }
    #endif

    return HEDLEY_STATIC_CAST(simde_bool, r);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i8x16_any_true(a) simde_wasm_i8x16_any_true((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_bool
simde_wasm_i16x8_any_true (simde_v128_t a) {
  return simde_wasm_i8x16_any_true(a);
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_any_true(a) simde_wasm_i16x8_any_true((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_bool
simde_wasm_i32x4_any_true (simde_v128_t a) {
  return simde_wasm_i8x16_any_true(a);
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_any_true(a) simde_wasm_i32x4_any_true((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_bool
simde_wasm_i64x2_any_true (simde_v128_t a) {
  return simde_wasm_i8x16_any_true(a);
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES) || (defined(SIMDE_ENABLE_NATIVE_ALIASES) && !defined(__wasm_unimplemented_simd128__))
  #define wasm_i64x2_any_true(a) simde_wasm_i64x2_any_true((a))
#endif

/* all_true */

SIMDE_FUNCTION_ATTRIBUTES
simde_bool
simde_wasm_i8x16_all_true (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i8x16_all_true(a);
  #else
    simde_v128_private a_ = simde_v128_to_private(a);

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      return _mm_test_all_zeros(_mm_cmpeq_epi8(a_.sse_m128i, _mm_set1_epi8(INT8_C(0))), _mm_set1_epi8(~INT8_C(0)));
    #else
      int8_t r = !INT8_C(0);

      SIMDE_VECTORIZE_REDUCTION(&:r)
      for (size_t i = 0 ; i < (sizeof(a_.i8) / sizeof(a_.i8[0])) ; i++) {
        r &= !!(a_.i8[i]);
      }

      return r;
    #endif
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i8x16_all_true(a) simde_wasm_i8x16_all_true((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_bool
simde_wasm_i16x8_all_true (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_all_true(a);
  #else
    simde_v128_private a_ = simde_v128_to_private(a);

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      return _mm_test_all_zeros(_mm_cmpeq_epi16(a_.sse_m128i, _mm_setzero_si128()), _mm_set1_epi16(~INT16_C(0)));
    #else
      int16_t r = !INT16_C(0);

      SIMDE_VECTORIZE_REDUCTION(&:r)
      for (size_t i = 0 ; i < (sizeof(a_.i16) / sizeof(a_.i16[0])) ; i++) {
        r &= !!(a_.i16[i]);
      }

      return r;
    #endif
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_all_true(a) simde_wasm_i16x8_all_true((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_bool
simde_wasm_i32x4_all_true (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i32x4_all_true(a);
  #else
    simde_v128_private a_ = simde_v128_to_private(a);

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      return _mm_test_all_zeros(_mm_cmpeq_epi32(a_.sse_m128i, _mm_setzero_si128()), _mm_set1_epi32(~INT32_C(0)));
    #else
      int32_t r = !INT32_C(0);

      SIMDE_VECTORIZE_REDUCTION(&:r)
      for (size_t i = 0 ; i < (sizeof(a_.i32) / sizeof(a_.i32[0])) ; i++) {
        r &= !!(a_.i32[i]);
      }

      return r;
    #endif
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_all_true(a) simde_wasm_i32x4_all_true((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_bool
simde_wasm_i64x2_all_true (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE) && defined(__wasm_unimplemented_simd128__)
    return wasm_i64x2_all_true(a);
  #else
    simde_v128_private a_ = simde_v128_to_private(a);

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      return _mm_test_all_zeros(_mm_cmpeq_epi64(a_.sse_m128i, _mm_setzero_si128()), _mm_set1_epi32(~INT32_C(0)));
    #else
      int64_t r = !INT32_C(0);

      SIMDE_VECTORIZE_REDUCTION(&:r)
      for (size_t i = 0 ; i < (sizeof(a_.i64) / sizeof(a_.i64[0])) ; i++) {
        r &= !!(a_.i64[i]);
      }

      return r;
    #endif
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES) || (defined(SIMDE_ENABLE_NATIVE_ALIASES) && !defined(__wasm_unimplemented_simd128__))
  #define wasm_i64x2_all_true(a) simde_wasm_i64x2_all_true((a))
#endif

/* shl
 *
 * Note: LLVM's implementation currently doesn't operate modulo
 * lane width, but the spec now says it should. */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i8x16_shl (simde_v128_t a, int32_t count) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i8x16_shl(a, count);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_VECTOR_SUBSCRIPT) && defined(SIMDE_VECTOR_SCALAR)
      r_.i8 = a_.i8 << (count & 7);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i8) / sizeof(r_.i8[0])) ; i++) {
        r_.i8[i] = HEDLEY_STATIC_CAST(int8_t, a_.i8[i] << (count & 7));
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i8x16_shl(a, count) simde_wasm_i8x16_shl((a), (count))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_shl (simde_v128_t a, int32_t count) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_shl(a, count);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      return _mm_sll_epi16(a_.sse_m128i, _mm_cvtsi32_si128(count & 15));
    #elif defined(SIMDE_VECTOR_SUBSCRIPT) && defined(SIMDE_VECTOR_SCALAR)
      r_.i16 = a_.i16 << (count & 15);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i16) / sizeof(r_.i16[0])) ; i++) {
        r_.i16[i] = HEDLEY_STATIC_CAST(int16_t, a_.i16[i] << (count & 15));
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_shl(a, count) simde_wasm_i16x8_shl((a), (count))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i32x4_shl (simde_v128_t a, int32_t count) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i32x4_shl(a, count);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      return _mm_sll_epi32(a_.sse_m128i, _mm_cvtsi32_si128(count & 31));
    #elif defined(SIMDE_VECTOR_SUBSCRIPT) && defined(SIMDE_VECTOR_SCALAR)
      r_.i32 = a_.i32 << (count & 31);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32) / sizeof(r_.i32[0])) ; i++) {
        r_.i32[i] = HEDLEY_STATIC_CAST(int32_t, a_.i32[i] << (count & 31));
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_shl(a, count) simde_wasm_i32x4_shl((a), (count))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i64x2_shl (simde_v128_t a, int32_t count) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i64x2_shl(a, count);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      return _mm_sll_epi64(a_.sse_m128i, _mm_cvtsi32_si128(count & 63));
    #elif defined(SIMDE_VECTOR_SUBSCRIPT) && defined(SIMDE_VECTOR_SCALAR)
      r_.i64 = a_.i64 << (count & 63);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i64) / sizeof(r_.i64[0])) ; i++) {
        r_.i64[i] = HEDLEY_STATIC_CAST(int64_t, a_.i64[i] << (count & 63));
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i64x2_shl(a, count) simde_wasm_i64x2_shl((a), (count))
#endif

/* shr */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i8x16_shr (simde_v128_t a, int32_t count) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i8x16_shr(a, count);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_VECTOR_SUBSCRIPT) && defined(SIMDE_VECTOR_SCALAR)
      r_.i8 = a_.i8 >> (count & 7);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i8) / sizeof(r_.i8[0])) ; i++) {
        r_.i8[i] = HEDLEY_STATIC_CAST(int8_t, a_.i8[i] >> (count & 7));
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i8x16_shr(a, count) simde_wasm_i8x16_shr((a), (count))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_shr (simde_v128_t a, int32_t count) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_shr(a, count);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      return _mm_sra_epi16(a_.sse_m128i, _mm_cvtsi32_si128(count & 15));
    #elif defined(SIMDE_VECTOR_SUBSCRIPT) && defined(SIMDE_VECTOR_SCALAR)
      r_.i16 = a_.i16 >> (count & 15);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i16) / sizeof(r_.i16[0])) ; i++) {
        r_.i16[i] = HEDLEY_STATIC_CAST(int16_t, a_.i16[i] >> (count & 15));
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_shr(a, count) simde_wasm_i16x8_shr((a), (count))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i32x4_shr (simde_v128_t a, int32_t count) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i32x4_shr(a, count);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      return _mm_sra_epi32(a_.sse_m128i, _mm_cvtsi32_si128(count & 31));
    #elif defined(SIMDE_VECTOR_SUBSCRIPT) && defined(SIMDE_VECTOR_SCALAR)
      r_.i32 = a_.i32 >> (count & 31);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32) / sizeof(r_.i32[0])) ; i++) {
        r_.i32[i] = HEDLEY_STATIC_CAST(int32_t, a_.i32[i] >> (count & 31));
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_shr(a, count) simde_wasm_i32x4_shr((a), (count))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i64x2_shr (simde_v128_t a, int32_t count) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i64x2_shr(a, count);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_X86_AVX512VL_NATIVE)
      return _mm_sra_epi64(a_.sse_m128i, _mm_cvtsi32_si128(count & 63));
    #elif defined(SIMDE_VECTOR_SUBSCRIPT) && defined(SIMDE_VECTOR_SCALAR)
      r_.i64 = a_.i64 >> (count & 63);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i64) / sizeof(r_.i64[0])) ; i++) {
        r_.i64[i] = HEDLEY_STATIC_CAST(int64_t, a_.i64[i] >> (count & 63));
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i64x2_shr(a, count) simde_wasm_i64x2_shr((a), (count))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u8x16_shr (simde_v128_t a, int32_t count) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u8x16_shr(a, count);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_VECTOR_SUBSCRIPT) && defined(SIMDE_VECTOR_SCALAR)
      r_.u8 = a_.u8 >> (count & 7);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u8) / sizeof(r_.u8[0])) ; i++) {
        r_.u8[i] = HEDLEY_STATIC_CAST(uint8_t, a_.u8[i] >> (count & 7));
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u8x16_shr(a, count) simde_wasm_u8x16_shr((a), (count))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u16x8_shr (simde_v128_t a, int32_t count) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u16x8_shr(a, count);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      return _mm_srl_epi16(a_.sse_m128i, _mm_cvtsi32_si128(count & 15));
    #elif defined(SIMDE_VECTOR_SUBSCRIPT) && defined(SIMDE_VECTOR_SCALAR)
      r_.u16 = a_.u16 >> (count & 15);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u16) / sizeof(r_.u16[0])) ; i++) {
        r_.u16[i] = HEDLEY_STATIC_CAST(uint16_t, a_.u16[i] >> (count & 15));
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u16x8_shr(a, count) simde_wasm_u16x8_shr((a), (count))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u32x4_shr (simde_v128_t a, int32_t count) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u32x4_shr(a, count);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      return _mm_srl_epi32(a_.sse_m128i, _mm_cvtsi32_si128(count & 31));
    #elif defined(SIMDE_VECTOR_SUBSCRIPT) && defined(SIMDE_VECTOR_SCALAR)
      r_.u32 = a_.u32 >> (count & 31);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u32) / sizeof(r_.u32[0])) ; i++) {
        r_.u32[i] = HEDLEY_STATIC_CAST(uint32_t, a_.u32[i] >> (count & 31));
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u32x4_shr(a, count) simde_wasm_u32x4_shr((a), (count))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u64x2_shr (simde_v128_t a, int32_t count) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u64x2_shr(a, count);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      return _mm_srl_epi64(a_.sse_m128i, _mm_cvtsi32_si128(count & 63));
    #elif defined(SIMDE_VECTOR_SUBSCRIPT) && defined(SIMDE_VECTOR_SCALAR)
      r_.u64 = a_.u64 >> (count & 63);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u64) / sizeof(r_.u64[0])) ; i++) {
        r_.u64[i] = HEDLEY_STATIC_CAST(uint64_t, a_.u64[i] >> (count & 63));
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u64x2_shr(a, count) simde_wasm_u64x2_shr((a), (count))
#endif

/* add */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i8x16_add (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i8x16_add(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_add_epi8(a_.sse_m128i, b_.sse_m128i);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i8 = a_.i8 + b_.i8;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i8) / sizeof(r_.i8[0])) ; i++) {
        r_.i8[i] = a_.i8[i] + b_.i8[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i8x16_add(a, b) simde_wasm_i8x16_add((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_add (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_add(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_add_epi16(a_.sse_m128i, b_.sse_m128i);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i16 = a_.i16 + b_.i16;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i16) / sizeof(r_.i16[0])) ; i++) {
        r_.i16[i] = a_.i16[i] + b_.i16[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_add(a, b) simde_wasm_i16x8_add((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i32x4_add (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i32x4_add(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_add_epi32(a_.sse_m128i, b_.sse_m128i);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32 = a_.i32 + b_.i32;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32) / sizeof(r_.i32[0])) ; i++) {
        r_.i32[i] = a_.i32[i] + b_.i32[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_add(a, b) simde_wasm_i32x4_add((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i64x2_add (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i64x2_add(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_add_epi64(a_.sse_m128i, b_.sse_m128i);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i64 = a_.i64 + b_.i64;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i64) / sizeof(r_.i64[0])) ; i++) {
        r_.i64[i] = a_.i64[i] + b_.i64[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i64x2_add(a, b) simde_wasm_i64x2_add((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f32x4_add (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f32x4_add(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128 = _mm_add_ps(a_.sse_m128, b_.sse_m128);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.f32 = a_.f32 + b_.f32;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f32) / sizeof(r_.f32[0])) ; i++) {
        r_.f32[i] = a_.f32[i] + b_.f32[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f32x4_add(a, b) simde_wasm_f32x4_add((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f64x2_add (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f64x2_add(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128d = _mm_add_pd(a_.sse_m128d, b_.sse_m128d);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.f64 = a_.f64 + b_.f64;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f64) / sizeof(r_.f64[0])) ; i++) {
        r_.f64[i] = a_.f64[i] + b_.f64[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f64x2_add(a, b) simde_wasm_f64x2_add((a), (b))
#endif

/* sub */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i8x16_sub (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i8x16_sub(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_sub_epi8(a_.sse_m128i, b_.sse_m128i);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i8 = a_.i8 - b_.i8;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i8) / sizeof(r_.i8[0])) ; i++) {
        r_.i8[i] = a_.i8[i] - b_.i8[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i8x16_sub(a, b) simde_wasm_i8x16_sub((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_sub (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_sub(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_sub_epi16(a_.sse_m128i, b_.sse_m128i);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i16 = a_.i16 - b_.i16;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i16) / sizeof(r_.i16[0])) ; i++) {
        r_.i16[i] = a_.i16[i] - b_.i16[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_sub(a, b) simde_wasm_i16x8_sub((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i32x4_sub (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i32x4_sub(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_sub_epi32(a_.sse_m128i, b_.sse_m128i);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32 = a_.i32 - b_.i32;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32) / sizeof(r_.i32[0])) ; i++) {
        r_.i32[i] = a_.i32[i] - b_.i32[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_sub(a, b) simde_wasm_i32x4_sub((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i64x2_sub (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i64x2_sub(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_sub_epi64(a_.sse_m128i, b_.sse_m128i);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i64 = a_.i64 - b_.i64;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i64) / sizeof(r_.i64[0])) ; i++) {
        r_.i64[i] = a_.i64[i] - b_.i64[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i64x2_sub(a, b) simde_wasm_i64x2_sub((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f32x4_sub (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f32x4_sub(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128 = _mm_sub_ps(a_.sse_m128, b_.sse_m128);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.f32 = a_.f32 - b_.f32;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f32) / sizeof(r_.f32[0])) ; i++) {
        r_.f32[i] = a_.f32[i] - b_.f32[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f32x4_sub(a, b) simde_wasm_f32x4_sub((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f64x2_sub (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f64x2_sub(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128d = _mm_sub_pd(a_.sse_m128d, b_.sse_m128d);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.f64 = a_.f64 - b_.f64;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f64) / sizeof(r_.f64[0])) ; i++) {
        r_.f64[i] = a_.f64[i] - b_.f64[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f64x2_sub(a, b) simde_wasm_f64x2_sub((a), (b))
#endif

/* mul */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_mul (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_mul(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_mullo_epi16(a_.sse_m128i, b_.sse_m128i);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i16 = a_.i16 * b_.i16;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i16) / sizeof(r_.i16[0])) ; i++) {
        r_.i16[i] = a_.i16[i] * b_.i16[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_mul(a, b) simde_wasm_i16x8_mul((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i32x4_mul (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i32x4_mul(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r_.sse_m128i = _mm_mullo_epi32(a_.sse_m128i, b_.sse_m128i);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i32 = a_.i32 * b_.i32;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32) / sizeof(r_.i32[0])) ; i++) {
        r_.i32[i] = a_.i32[i] * b_.i32[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_mul(a, b) simde_wasm_i32x4_mul((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i64x2_mul (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i64x2_mul(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_AVX512VL_NATIVE) && defined(SIMDE_X86_AVX512DQ_NATIVE)
      r_.sse_m128i = _mm_mullo_epi64(a_.sse_m128i, b_.sse_m128i);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.i64 = a_.i64 * b_.i64;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i64) / sizeof(r_.i64[0])) ; i++) {
        r_.i64[i] = a_.i64[i] * b_.i64[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i64x2_mul(a, b) simde_wasm_i64x2_mul((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f32x4_mul (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f32x4_mul(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128 = _mm_mul_ps(a_.sse_m128, b_.sse_m128);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.f32 = a_.f32 * b_.f32;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f32) / sizeof(r_.f32[0])) ; i++) {
        r_.f32[i] = a_.f32[i] * b_.f32[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f32x4_mul(a, b) simde_wasm_f32x4_mul((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f64x2_mul (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f64x2_mul(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128d = _mm_mul_pd(a_.sse_m128d, b_.sse_m128d);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.f64 = a_.f64 * b_.f64;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f64) / sizeof(r_.f64[0])) ; i++) {
        r_.f64[i] = a_.f64[i] * b_.f64[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f64x2_mul(a, b) simde_wasm_f64x2_mul((a), (b))
#endif

/* min */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i8x16_min (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i8x16_min(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r_.sse_m128i = _mm_min_epi8(a_.sse_m128i, b_.sse_m128i);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i8) / sizeof(r_.i8[0])) ; i++) {
        r_.i8[i] = (a_.i8[i] < b_.i8[i]) ? a_.i8[i] : b_.i8[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i8x16_min(a, b) simde_wasm_i8x16_min((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_min (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_min(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_min_epi16(a_.sse_m128i, b_.sse_m128i);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i16) / sizeof(r_.i16[0])) ; i++) {
        r_.i16[i] = (a_.i16[i] < b_.i16[i]) ? a_.i16[i] : b_.i16[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_min(a, b) simde_wasm_i16x8_min((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i32x4_min (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i32x4_min(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r_.sse_m128i = _mm_min_epi32(a_.sse_m128i, b_.sse_m128i);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32) / sizeof(r_.i32[0])) ; i++) {
        r_.i32[i] = (a_.i32[i] < b_.i32[i]) ? a_.i32[i] : b_.i32[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_min(a, b) simde_wasm_i32x4_min((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u8x16_min (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u8x16_min(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_min_epu8(a_.sse_m128i, b_.sse_m128i);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u8) / sizeof(r_.u8[0])) ; i++) {
        r_.u8[i] = (a_.u8[i] < b_.u8[i]) ? a_.u8[i] : b_.u8[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u8x16_min(a, b) simde_wasm_u8x16_min((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u16x8_min (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u16x8_min(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r_.sse_m128i = _mm_min_epu16(a_.sse_m128i, b_.sse_m128i);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u16) / sizeof(r_.u16[0])) ; i++) {
        r_.u16[i] = (a_.u16[i] < b_.u16[i]) ? a_.u16[i] : b_.u16[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u16x8_min(a, b) simde_wasm_u16x8_min((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u32x4_min (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u32x4_min(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r_.sse_m128i = _mm_min_epu32(a_.sse_m128i, b_.sse_m128i);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32) / sizeof(r_.i32[0])) ; i++) {
        r_.u32[i] = (a_.u32[i] < b_.u32[i]) ? a_.u32[i] : b_.u32[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u32x4_min(a, b) simde_wasm_u32x4_min((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f32x4_min (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f32x4_min(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r_.sse_m128 = _mm_blendv_ps(
          _mm_set1_ps(SIMDE_MATH_NANF),
          _mm_min_ps(a_.sse_m128, b_.sse_m128),
          _mm_cmpord_ps(a_.sse_m128, b_.sse_m128));
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f32) / sizeof(r_.f32[0])) ; i++) {
        r_.f32[i] = simde_math_isnan(a_.f32[i]) ? a_.f32[i] : ((a_.f32[i] < b_.f32[i]) ? a_.f32[i] : b_.f32[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f32x4_min(a, b) simde_wasm_f32x4_min((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f64x2_min (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f64x2_min(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r_.sse_m128d = _mm_blendv_pd(
          _mm_set1_pd(SIMDE_MATH_NAN),
          _mm_min_pd(a_.sse_m128d, b_.sse_m128d),
          _mm_cmpord_pd(a_.sse_m128d, b_.sse_m128d));
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f64) / sizeof(r_.f64[0])) ; i++) {
        r_.f64[i] = simde_math_isnan(a_.f64[i]) ? a_.f64[i] : ((a_.f64[i] < b_.f64[i]) ? a_.f64[i] : b_.f64[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f64x2_min(a, b) simde_wasm_f64x2_min((a), (b))
#endif

/* max */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i8x16_max (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i8x16_max(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r_.sse_m128i = _mm_max_epi8(a_.sse_m128i, b_.sse_m128i);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i8) / sizeof(r_.i8[0])) ; i++) {
        r_.i8[i] = (a_.i8[i] > b_.i8[i]) ? a_.i8[i] : b_.i8[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i8x16_max(a, b) simde_wasm_i8x16_max((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_max (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_max(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_max_epi16(a_.sse_m128i, b_.sse_m128i);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i16) / sizeof(r_.i16[0])) ; i++) {
        r_.i16[i] = (a_.i16[i] > b_.i16[i]) ? a_.i16[i] : b_.i16[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_max(a, b) simde_wasm_i16x8_max((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i32x4_max (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i32x4_max(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r_.sse_m128i = _mm_max_epi32(a_.sse_m128i, b_.sse_m128i);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32) / sizeof(r_.i32[0])) ; i++) {
        r_.i32[i] = (a_.i32[i] > b_.i32[i]) ? a_.i32[i] : b_.i32[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_max(a, b) simde_wasm_i32x4_max((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u8x16_max (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u8x16_max(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_max_epu8(a_.sse_m128i, b_.sse_m128i);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u8) / sizeof(r_.u8[0])) ; i++) {
        r_.u8[i] = (a_.u8[i] > b_.u8[i]) ? a_.u8[i] : b_.u8[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u8x16_max(a, b) simde_wasm_u8x16_max((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u16x8_max (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u16x8_max(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r_.sse_m128i = _mm_max_epu16(a_.sse_m128i, b_.sse_m128i);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u16) / sizeof(r_.u16[0])) ; i++) {
        r_.u16[i] = (a_.u16[i] > b_.u16[i]) ? a_.u16[i] : b_.u16[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u16x8_max(a, b) simde_wasm_u16x8_max((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u32x4_max (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u32x4_max(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r_.sse_m128i = _mm_max_epu32(a_.sse_m128i, b_.sse_m128i);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32) / sizeof(r_.i32[0])) ; i++) {
        r_.u32[i] = (a_.u32[i] > b_.u32[i]) ? a_.u32[i] : b_.u32[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u32x4_max(a, b) simde_wasm_u32x4_max((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f32x4_max (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f32x4_max(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r_.sse_m128 = _mm_blendv_ps(
          _mm_set1_ps(SIMDE_MATH_NANF),
          _mm_max_ps(a_.sse_m128, b_.sse_m128),
          _mm_cmpord_ps(a_.sse_m128, b_.sse_m128));
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f32) / sizeof(r_.f32[0])) ; i++) {
        r_.f32[i] = simde_math_isnan(a_.f32[i]) ? a_.f32[i] : ((a_.f32[i] > b_.f32[i]) ? a_.f32[i] : b_.f32[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f32x4_max(a, b) simde_wasm_f32x4_max((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f64x2_max (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f64x2_max(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE4_1_NATIVE)
      r_.sse_m128d = _mm_blendv_pd(
          _mm_set1_pd(SIMDE_MATH_NAN),
          _mm_max_pd(a_.sse_m128d, b_.sse_m128d),
          _mm_cmpord_pd(a_.sse_m128d, b_.sse_m128d));
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f64) / sizeof(r_.f64[0])) ; i++) {
        r_.f64[i] = simde_math_isnan(a_.f64[i]) ? a_.f64[i] : ((a_.f64[i] > b_.f64[i]) ? a_.f64[i] : b_.f64[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f64x2_max(a, b) simde_wasm_f64x2_max((a), (b))
#endif

/* add_saturate */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i8x16_add_saturate (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i8x16_add_saturate(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_adds_epi8(a_.sse_m128i, b_.sse_m128i);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i8) / sizeof(r_.i8[0])) ; i++) {
        r_.i8[i] = simde_math_adds_i8(a_.i8[i], b_.i8[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i8x16_add_saturate(a, b) simde_wasm_i8x16_add_saturate((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_add_saturate (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_add_saturate(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_adds_epi16(a_.sse_m128i, b_.sse_m128i);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i16) / sizeof(r_.i16[0])) ; i++) {
        r_.i16[i] = simde_math_adds_i16(a_.i16[i], b_.i16[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_add_saturate(a, b) simde_wasm_i16x8_add_saturate((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u8x16_add_saturate (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u8x16_add_saturate(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_adds_epu8(a_.sse_m128i, b_.sse_m128i);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u8) / sizeof(r_.u8[0])) ; i++) {
        r_.u8[i] = simde_math_adds_u8(a_.u8[i], b_.u8[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u8x16_add_saturate(a, b) simde_wasm_u8x16_add_saturate((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u16x8_add_saturate (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u16x8_add_saturate(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_adds_epu16(a_.sse_m128i, b_.sse_m128i);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u16) / sizeof(r_.u16[0])) ; i++) {
        r_.u16[i] = simde_math_adds_u16(a_.u16[i], b_.u16[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u16x8_add_saturate(a, b) simde_wasm_u16x8_add_saturate((a), (b))
#endif

/* avgr */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u8x16_avgr (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u8x16_avgr(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_avg_epu8(a_.sse_m128i, b_.sse_m128i);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u8) / sizeof(r_.u8[0])) ; i++) {
        r_.u8[i] = (a_.u8[i] + b_.u8[i] + 1) >> 1;
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u8x16_avgr(a, b) simde_wasm_u8x16_avgr((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u16x8_avgr (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u16x8_avgr(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_avg_epu16(a_.sse_m128i, b_.sse_m128i);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u16) / sizeof(r_.u16[0])) ; i++) {
        r_.u16[i] = (a_.u16[i] + b_.u16[i] + 1) >> 1;
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u16x8_avgr(a, b) simde_wasm_u16x8_avgr((a), (b))
#endif

/* sub_saturate */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i8x16_sub_saturate (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i8x16_sub_saturate(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_subs_epi8(a_.sse_m128i, b_.sse_m128i);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i8) / sizeof(r_.i8[0])) ; i++) {
        r_.i8[i] = simde_math_subs_i8(a_.i8[i], b_.i8[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i8x16_sub_saturate(a, b) simde_wasm_i8x16_sub_saturate((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_sub_saturate (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_sub_saturate(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_subs_epi16(a_.sse_m128i, b_.sse_m128i);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i16) / sizeof(r_.i16[0])) ; i++) {
        r_.i16[i] = simde_math_subs_i16(a_.i16[i], b_.i16[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_sub_saturate(a, b) simde_wasm_i16x8_sub_saturate((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u8x16_sub_saturate (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u8x16_sub_saturate(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_subs_epu8(a_.sse_m128i, b_.sse_m128i);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u8) / sizeof(r_.u8[0])) ; i++) {
        r_.u8[i] = simde_math_subs_u8(a_.u8[i], b_.u8[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u8x16_sub_saturate(a, b) simde_wasm_u8x16_sub_saturate((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u16x8_sub_saturate (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u16x8_sub_saturate(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128i = _mm_subs_epu16(a_.sse_m128i, b_.sse_m128i);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u16) / sizeof(r_.u16[0])) ; i++) {
        r_.u16[i] = simde_math_subs_u16(a_.u16[i], b_.u16[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u16x8_sub_saturate(a, b) simde_wasm_u16x8_sub_saturate((a), (b))
#endif

/* pmin */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f32x4_pmin (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f32x4_pmin(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128 = _mm_min_ps(b_.sse_m128, a_.sse_m128);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f32) / sizeof(r_.f32[0])) ; i++) {
        r_.f32[i] = (b_.f32[i] < a_.f32[i]) ? b_.f32[i] : a_.f32[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f32x4_pmin(a, b) simde_wasm_f32x4_pmin((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f64x2_pmin (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f64x2_pmin(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128d = _mm_min_pd(b_.sse_m128d, a_.sse_m128d);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f64) / sizeof(r_.f64[0])) ; i++) {
        r_.f64[i] = (b_.f64[i] < a_.f64[i]) ? b_.f64[i] : a_.f64[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f64x2_pmin(a, b) simde_wasm_f64x2_pmin((a), (b))
#endif

/* pmax */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f32x4_pmax (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f32x4_pmax(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128 = _mm_max_ps(b_.sse_m128, a_.sse_m128);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f32) / sizeof(r_.f32[0])) ; i++) {
        r_.f32[i] = (a_.f32[i] < b_.f32[i]) ? b_.f32[i] : a_.f32[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f32x4_pmax(a, b) simde_wasm_f32x4_pmax((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f64x2_pmax (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f64x2_pmax(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128d = _mm_max_pd(b_.sse_m128d, a_.sse_m128d);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f64) / sizeof(r_.f64[0])) ; i++) {
        r_.f64[i] = (a_.f64[i] < b_.f64[i]) ? b_.f64[i] : a_.f64[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f64x2_pmax(a, b) simde_wasm_f64x2_pmax((a), (b))
#endif

/* div */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f32x4_div (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f32x4_div(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128 = _mm_div_ps(a_.sse_m128, b_.sse_m128);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.f32 = a_.f32 / b_.f32;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f32) / sizeof(r_.f32[0])) ; i++) {
        r_.f32[i] = a_.f32[i] / b_.f32[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f32x4_div(a, b) simde_wasm_f32x4_div((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f64x2_div (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f64x2_div(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_X86_SSE2_NATIVE)
      r_.sse_m128d = _mm_div_pd(a_.sse_m128d, b_.sse_m128d);
    #elif defined(SIMDE_VECTOR_SUBSCRIPT)
      r_.f64 = a_.f64 / b_.f64;
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f64) / sizeof(r_.f64[0])) ; i++) {
        r_.f64[i] = a_.f64[i] / b_.f64[i];
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f64x2_div(a, b) simde_wasm_f64x2_div((a), (b))
#endif

/* shuffle */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_v8x16_shuffle (
    simde_v128_t a, simde_v128_t b,
    const int c0, const int c1, const int  c2, const int  c3, const int  c4, const int  c5, const int  c6, const int  c7,
    const int c8, const int c9, const int c10, const int c11, const int c12, const int c13, const int c14, const int c15) {
  simde_v128_private
    a_ = simde_v128_to_private(a),
    b_ = simde_v128_to_private(b),
    r_;

  r_.i8[ 0] = ( c0 < 16) ? a_.i8[ c0] : b_.i8[ c0 & 15];
  r_.i8[ 1] = ( c1 < 16) ? a_.i8[ c1] : b_.i8[ c1 & 15];
  r_.i8[ 2] = ( c2 < 16) ? a_.i8[ c2] : b_.i8[ c2 & 15];
  r_.i8[ 3] = ( c3 < 16) ? a_.i8[ c3] : b_.i8[ c3 & 15];
  r_.i8[ 4] = ( c4 < 16) ? a_.i8[ c4] : b_.i8[ c4 & 15];
  r_.i8[ 5] = ( c5 < 16) ? a_.i8[ c5] : b_.i8[ c5 & 15];
  r_.i8[ 6] = ( c6 < 16) ? a_.i8[ c6] : b_.i8[ c6 & 15];
  r_.i8[ 7] = ( c7 < 16) ? a_.i8[ c7] : b_.i8[ c7 & 15];
  r_.i8[ 8] = ( c8 < 16) ? a_.i8[ c8] : b_.i8[ c8 & 15];
  r_.i8[ 9] = ( c9 < 16) ? a_.i8[ c9] : b_.i8[ c9 & 15];
  r_.i8[10] = (c10 < 16) ? a_.i8[c10] : b_.i8[c10 & 15];
  r_.i8[11] = (c11 < 16) ? a_.i8[c11] : b_.i8[c11 & 15];
  r_.i8[12] = (c12 < 16) ? a_.i8[c12] : b_.i8[c12 & 15];
  r_.i8[13] = (c13 < 16) ? a_.i8[c13] : b_.i8[c13 & 15];
  r_.i8[14] = (c14 < 16) ? a_.i8[c14] : b_.i8[c14 & 15];
  r_.i8[15] = (c15 < 16) ? a_.i8[c15] : b_.i8[c15 & 15];

  return simde_v128_from_private(r_);
}
#if defined(SIMDE_WASM_SIMD128_NATIVE)
  #define \
    simde_wasm_v8x16_shuffle( \
        a, b, \
        c0, c1,  c2,  c3,  c4,  c5,  c6,  c7, \
        c8, c9, c10, c11, c12, c13, c14, c15) \
    wasm_v8x16_shuffle( \
        a, b, \
        c0, c1,  c2,  c3,  c4,  c5,  c6,  c7, \
        c8, c9, c10, c11, c12, c13, c14, c15)
#elif defined(SIMDE_SHUFFLE_VECTOR_)
  #define \
    simde_wasm_v8x16_shuffle( \
        a, b, \
        c0, c1,  c2,  c3,  c4,  c5,  c6,  c7, \
        c8, c9, c10, c11, c12, c13, c14, c15) \
    (__extension__ ({ \
      HEDLEY_REINTERPRET_CAST(simde_v128_t, SIMDE_SHUFFLE_VECTOR_(8, 16, \
          HEDLEY_REINTERPRET_CAST(int8_t SIMDE_VECTOR(16), a), \
          HEDLEY_REINTERPRET_CAST(int8_t SIMDE_VECTOR(16), b), \
          c0, c1,  c2,  c3,  c4,  c5,  c6,  c7, \
          c8, c9, c10, c11, c12, c13, c14, c15)); \
    }))
#endif
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define \
    wasm_v8x16_shuffle(a, b, \
        c0, c1,  c2,  c3,  c4,  c5,  c6,  c7, \
        c8, c9, c10, c11, c12, c13, c14, c15) \
      simde_wasm_v8x16_shuffle((a), (b), \
          (c0), (c1),  (c2),  (c3),  (c4),  (c5),  (c6),  (c7), \
          (c8), (c9), (c10), (c11), (c12), (c13), (c14), (c15))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_v16x8_shuffle (
    simde_v128_t a, simde_v128_t b,
    const int c0, const int c1, const int  c2, const int  c3, const int  c4, const int  c5, const int  c6, const int  c7) {
  simde_v128_private
    a_ = simde_v128_to_private(a),
    b_ = simde_v128_to_private(b),
    r_;

  r_.i16[ 0] = (c0 < 8) ? a_.i16[ c0] : b_.i16[ c0 & 7];
  r_.i16[ 1] = (c1 < 8) ? a_.i16[ c1] : b_.i16[ c1 & 7];
  r_.i16[ 2] = (c2 < 8) ? a_.i16[ c2] : b_.i16[ c2 & 7];
  r_.i16[ 3] = (c3 < 8) ? a_.i16[ c3] : b_.i16[ c3 & 7];
  r_.i16[ 4] = (c4 < 8) ? a_.i16[ c4] : b_.i16[ c4 & 7];
  r_.i16[ 5] = (c5 < 8) ? a_.i16[ c5] : b_.i16[ c5 & 7];
  r_.i16[ 6] = (c6 < 8) ? a_.i16[ c6] : b_.i16[ c6 & 7];
  r_.i16[ 7] = (c7 < 8) ? a_.i16[ c7] : b_.i16[ c7 & 7];

  return simde_v128_from_private(r_);
}
#if defined(SIMDE_WASM_SIMD128_NATIVE)
  #define \
    simde_wasm_v16x8_shuffle( \
        a, b, \
        c0, c1,  c2,  c3,  c4,  c5,  c6,  c7) \
    wasm_v16x8_shuffle( \
        a, b, \
        c0, c1,  c2,  c3,  c4,  c5,  c6,  c7)
#elif defined(SIMDE_SHUFFLE_VECTOR_)
  #define \
    simde_wasm_v16x8_shuffle( \
        a, b, \
        c0, c1,  c2,  c3,  c4,  c5,  c6,  c7) \
    (__extension__ ({ \
      HEDLEY_REINTERPRET_CAST(simde_v128_t, SIMDE_SHUFFLE_VECTOR_(16, 16, \
          HEDLEY_REINTERPRET_CAST(int16_t SIMDE_VECTOR(16), a), \
          HEDLEY_REINTERPRET_CAST(int16_t SIMDE_VECTOR(16), b), \
          c0, c1,  c2,  c3,  c4,  c5,  c6,  c7)); \
    }))
#endif
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define \
    wasm_v16x8_shuffle(a, b, \
        c0, c1,  c2,  c3,  c4,  c5,  c6,  c7) \
      simde_wasm_v16x8_shuffle((a), (b), \
          (c0), (c1),  (c2),  (c3),  (c4),  (c5),  (c6),  (c7))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_v32x4_shuffle (
    simde_v128_t a, simde_v128_t b,
    const int c0, const int c1, const int  c2, const int  c3) {
  simde_v128_private
    a_ = simde_v128_to_private(a),
    b_ = simde_v128_to_private(b),
    r_;

  r_.i32[ 0] = (c0 < 4) ? a_.i32[ c0] : b_.i32[ c0 & 3];
  r_.i32[ 1] = (c1 < 4) ? a_.i32[ c1] : b_.i32[ c1 & 3];
  r_.i32[ 2] = (c2 < 4) ? a_.i32[ c2] : b_.i32[ c2 & 3];
  r_.i32[ 3] = (c3 < 4) ? a_.i32[ c3] : b_.i32[ c3 & 3];

  return simde_v128_from_private(r_);
}
#if defined(SIMDE_WASM_SIMD128_NATIVE)
  #define \
    simde_wasm_v32x4_shuffle( \
        a, b, \
        c0, c1,  c2,  c3) \
    wasm_v32x4_shuffle( \
        a, b, \
        c0, c1,  c2,  c3)
#elif defined(SIMDE_SHUFFLE_VECTOR_)
  #define \
    simde_wasm_v32x4_shuffle( \
        a, b, \
        c0, c1,  c2,  c3) \
    (__extension__ ({ \
      HEDLEY_REINTERPRET_CAST(simde_v128_t, SIMDE_SHUFFLE_VECTOR_(32, 16, \
          HEDLEY_REINTERPRET_CAST(int32_t SIMDE_VECTOR(16), a), \
          HEDLEY_REINTERPRET_CAST(int32_t SIMDE_VECTOR(16), b), \
          c0, c1,  c2,  c3)); \
    }))
#endif
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define \
    wasm_v32x4_shuffle(a, b, \
        c0, c1,  c2,  c3) \
      simde_wasm_v32x4_shuffle((a), (b), \
          (c0), (c1),  (c2),  (c3))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_v64x2_shuffle (
    simde_v128_t a, simde_v128_t b,
    const int c0, const int c1) {
  simde_v128_private
    a_ = simde_v128_to_private(a),
    b_ = simde_v128_to_private(b),
    r_;

  r_.i64[ 0] = (c0 < 2) ? a_.i64[ c0] : b_.i64[ c0 & 1];
  r_.i64[ 1] = (c1 < 2) ? a_.i64[ c1] : b_.i64[ c1 & 1];

  return simde_v128_from_private(r_);
}
#if defined(SIMDE_WASM_SIMD128_NATIVE)
  #define \
    simde_wasm_v64x2_shuffle( \
        a, b, \
        c0, c1) \
    wasm_v64x2_shuffle( \
        a, b, \
        c0, c1)
#elif defined(SIMDE_SHUFFLE_VECTOR_)
  #define \
    simde_wasm_v64x2_shuffle( \
        a, b, \
        c0, c1) \
    (__extension__ ({ \
      HEDLEY_REINTERPRET_CAST(simde_v128_t, SIMDE_SHUFFLE_VECTOR_(64, 16, \
          HEDLEY_REINTERPRET_CAST(int64_t SIMDE_VECTOR(16), a), \
          HEDLEY_REINTERPRET_CAST(int64_t SIMDE_VECTOR(16), b), \
          c0, c1)); \
    }))
#endif
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define \
    wasm_v64x2_shuffle(a, b, \
        c0, c1) \
      simde_wasm_v64x2_shuffle((a), (b), \
          (c0), (c1))
#endif

/* swizzle */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_v8x16_swizzle (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_v8x16_swizzle(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    SIMDE_VECTORIZE
    for (size_t i = 0 ; i < (sizeof(r_.i8) / sizeof(r_.i8[0])) ; i++) {
      r_.i8[i] = ((b_.i8[i] & 15) == b_.i8[i]) ? a_.i8[b_.i8[i]] : INT8_C(0);
    }

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_v8x16_swizzle(a, b) simde_wasm_v8x16_swizzle((a), (b))
#endif

/* narrow */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i8x16_narrow_i16x8 (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i8x16_narrow_i16x8(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_CONVERT_VECTOR_) && HEDLEY_HAS_BUILTIN(__builtin_shufflevector)
      int16_t v SIMDE_VECTOR(32) = SIMDE_SHUFFLE_VECTOR_(16, 32, a_.i16, b_.i16, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
      const int16_t min SIMDE_VECTOR(32) = { INT8_MIN, INT8_MIN, INT8_MIN, INT8_MIN, INT8_MIN, INT8_MIN, INT8_MIN, INT8_MIN, INT8_MIN, INT8_MIN, INT8_MIN, INT8_MIN, INT8_MIN, INT8_MIN, INT8_MIN, INT8_MIN };
      const int16_t max SIMDE_VECTOR(32) = { INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX };

      int16_t m SIMDE_VECTOR(32);
      m = HEDLEY_REINTERPRET_CAST(__typeof__(m), v < min);
      v = (v & ~m) | (min & m);

      m = v > max;
      v = (v & ~m) | (max & m);

      SIMDE_CONVERT_VECTOR_(r_.i8, v);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i8) / sizeof(r_.i8[0])) ; i++) {
        int16_t v = (i < (sizeof(a_.i16) / sizeof(a_.i16[0]))) ? a_.i16[i] : b_.i16[i & 7];
        r_.i8[i] = (v < INT8_MIN) ? INT8_MIN : ((v > INT8_MAX) ? INT8_MAX : HEDLEY_STATIC_CAST(int8_t, v));
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i8x16_narrow_i16x8(a, b) simde_wasm_i8x16_narrow_i16x8((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_narrow_i32x4 (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_narrow_i32x4(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_CONVERT_VECTOR_) && HEDLEY_HAS_BUILTIN(__builtin_shufflevector)
      int32_t v SIMDE_VECTOR(32) = SIMDE_SHUFFLE_VECTOR_(32, 32, a_.i32, b_.i32, 0, 1, 2, 3, 4, 5, 6, 7);
      const int32_t min SIMDE_VECTOR(32) = { INT16_MIN, INT16_MIN, INT16_MIN, INT16_MIN, INT16_MIN, INT16_MIN, INT16_MIN, INT16_MIN };
      const int32_t max SIMDE_VECTOR(32) = { INT16_MAX, INT16_MAX, INT16_MAX, INT16_MAX, INT16_MAX, INT16_MAX, INT16_MAX, INT16_MAX };

      int32_t m SIMDE_VECTOR(32);
      m = HEDLEY_REINTERPRET_CAST(__typeof__(m), v < min);
      v = (v & ~m) | (min & m);

      m = HEDLEY_REINTERPRET_CAST(__typeof__(m), v > max);
      v = (v & ~m) | (max & m);

      SIMDE_CONVERT_VECTOR_(r_.i16, v);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i16) / sizeof(r_.i16[0])) ; i++) {
        int32_t v = (i < (sizeof(a_.i32) / sizeof(a_.i32[0]))) ? a_.i32[i] : b_.i32[i & 3];
        r_.i16[i] = (v < INT16_MIN) ? INT16_MIN : ((v > INT16_MAX) ? INT16_MAX : HEDLEY_STATIC_CAST(int16_t, v));
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_narrow_i32x4(a, b) simde_wasm_i16x8_narrow_i32x4((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u8x16_narrow_i16x8 (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u8x16_narrow_i16x8(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_CONVERT_VECTOR_) && HEDLEY_HAS_BUILTIN(__builtin_shufflevector)
      int16_t v SIMDE_VECTOR(32) = SIMDE_SHUFFLE_VECTOR_(16, 32, a_.i16, b_.i16, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
      const int16_t min SIMDE_VECTOR(32) = { 0, };
      const int16_t max SIMDE_VECTOR(32) = { UINT8_MAX, UINT8_MAX, UINT8_MAX, UINT8_MAX, UINT8_MAX, UINT8_MAX, UINT8_MAX, UINT8_MAX, UINT8_MAX, UINT8_MAX, UINT8_MAX, UINT8_MAX, UINT8_MAX, UINT8_MAX, UINT8_MAX, UINT8_MAX };

      int16_t m SIMDE_VECTOR(32);
      m = HEDLEY_REINTERPRET_CAST(__typeof__(m), v < min);
      v = (v & ~m) | (min & m);

      m = HEDLEY_REINTERPRET_CAST(__typeof__(m), v > max);
      v = (v & ~m) | (max & m);

      SIMDE_CONVERT_VECTOR_(r_.i8, v);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i8) / sizeof(r_.i8[0])) ; i++) {
        int16_t v = (i < (sizeof(a_.i16) / sizeof(a_.i16[0]))) ? a_.i16[i] : b_.i16[i & 7];
        r_.u8[i] = (v < 0) ? UINT8_C(0) : ((v > UINT8_MAX) ? UINT8_MAX : HEDLEY_STATIC_CAST(uint8_t, v));
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u8x16_narrow_i16x8(a, b) simde_wasm_u8x16_narrow_i16x8((a), (b))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u16x8_narrow_i32x4 (simde_v128_t a, simde_v128_t b) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u16x8_narrow_i32x4(a, b);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      b_ = simde_v128_to_private(b),
      r_;

    #if defined(SIMDE_CONVERT_VECTOR_) && HEDLEY_HAS_BUILTIN(__builtin_shufflevector)
      int32_t v SIMDE_VECTOR(32) = SIMDE_SHUFFLE_VECTOR_(32, 32, a_.i32, b_.i32, 0, 1, 2, 3, 4, 5, 6, 7);
      const int32_t min SIMDE_VECTOR(32) = { 0, };
      const int32_t max SIMDE_VECTOR(32) = { UINT16_MAX, UINT16_MAX, UINT16_MAX, UINT16_MAX, UINT16_MAX, UINT16_MAX, UINT16_MAX, UINT16_MAX };

      int32_t m SIMDE_VECTOR(32);
      m = HEDLEY_REINTERPRET_CAST(__typeof__(m), v < min);
      v = (v & ~m) | (min & m);

      m = HEDLEY_REINTERPRET_CAST(__typeof__(m), v > max);
      v = (v & ~m) | (max & m);

      SIMDE_CONVERT_VECTOR_(r_.i16, v);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i16) / sizeof(r_.i16[0])) ; i++) {
        int32_t v = (i < (sizeof(a_.i32) / sizeof(a_.i32[0]))) ? a_.i32[i] : b_.i32[i & 3];
        r_.u16[i] = (v < 0) ? UINT16_C(0) : ((v > UINT16_MAX) ? UINT16_MAX : HEDLEY_STATIC_CAST(uint16_t, v));
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u16x8_narrow_i32x4(a, b) simde_wasm_u16x8_narrow_i32x4((a), (b))
#endif

/* widen_low */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_widen_low_i8x16 (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_widen_low_i8x16(a);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_CONVERT_VECTOR_)
      const int8_t v SIMDE_VECTOR(8) = {
        a_.i8[0], a_.i8[1], a_.i8[2], a_.i8[3],
        a_.i8[4], a_.i8[5], a_.i8[6], a_.i8[7]
      };

      SIMDE_CONVERT_VECTOR_(r_.i16, v);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i16) / sizeof(r_.i16[0])) ; i++) {
        r_.i16[i] = HEDLEY_STATIC_CAST(int16_t, a_.i8[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_widen_low_i8x16(a) simde_wasm_i16x8_widen_low_i8x16((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i32x4_widen_low_i16x8 (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i32x4_widen_low_i16x8(a);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_CONVERT_VECTOR_)
      const int16_t v SIMDE_VECTOR(8) = { a_.i16[0], a_.i16[1], a_.i16[2], a_.i16[3] };

      SIMDE_CONVERT_VECTOR_(r_.i32, v);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32) / sizeof(r_.i32[0])) ; i++) {
        r_.i32[i] = HEDLEY_STATIC_CAST(int32_t, a_.i16[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_widen_low_i16x8(a) simde_wasm_i32x4_widen_low_i16x8((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_widen_low_u8x16 (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_widen_low_u8x16(a);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_CONVERT_VECTOR_)
      const uint8_t v SIMDE_VECTOR(8) = {
        a_.u8[0], a_.u8[1], a_.u8[2], a_.u8[3],
        a_.u8[4], a_.u8[5], a_.u8[6], a_.u8[7]
      };

      SIMDE_CONVERT_VECTOR_(r_.i16, v);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i16) / sizeof(r_.i16[0])) ; i++) {
        r_.i16[i] = HEDLEY_STATIC_CAST(int16_t, a_.u8[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_widen_low_u8x16(a) simde_wasm_i16x8_widen_low_u8x16((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i32x4_widen_low_u16x8 (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i32x4_widen_low_u16x8(a);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_CONVERT_VECTOR_)
      const uint16_t v SIMDE_VECTOR(8) = { a_.u16[0], a_.u16[1], a_.u16[2], a_.u16[3] };

      SIMDE_CONVERT_VECTOR_(r_.i32, v);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32) / sizeof(r_.i32[0])) ; i++) {
        r_.i32[i] = HEDLEY_STATIC_CAST(int32_t, a_.u16[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_widen_low_u16x8(a) simde_wasm_i32x4_widen_low_u16x8((a))
#endif

/* widen_high */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_widen_high_i8x16 (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_widen_high_i8x16(a);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_CONVERT_VECTOR_)
      const int8_t v SIMDE_VECTOR(8) = {
        a_.i8[ 8], a_.i8[ 9], a_.i8[10], a_.i8[11],
        a_.i8[12], a_.i8[13], a_.i8[14], a_.i8[15]
      };

      SIMDE_CONVERT_VECTOR_(r_.i16, v);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i16) / sizeof(r_.i16[0])) ; i++) {
        r_.i16[i] = HEDLEY_STATIC_CAST(int16_t, a_.i8[i + 8]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_widen_high_i8x16(a) simde_wasm_i16x8_widen_high_i8x16((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i32x4_widen_high_i16x8 (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i32x4_widen_high_i16x8(a);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_CONVERT_VECTOR_)
      const int16_t v SIMDE_VECTOR(8) = { a_.i16[4], a_.i16[5], a_.i16[6], a_.i16[7] };

      SIMDE_CONVERT_VECTOR_(r_.i32, v);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32) / sizeof(r_.i32[0])) ; i++) {
        r_.i32[i] = HEDLEY_STATIC_CAST(int32_t, a_.i16[i + 4]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_widen_high_i16x8(a) simde_wasm_i32x4_widen_high_i16x8((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_widen_high_u8x16 (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_widen_high_u8x16(a);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_CONVERT_VECTOR_)
      const uint8_t v SIMDE_VECTOR(8) = {
        a_.u8[ 8], a_.u8[ 9], a_.u8[10], a_.u8[11],
        a_.u8[12], a_.u8[13], a_.u8[14], a_.u8[15]
      };

      SIMDE_CONVERT_VECTOR_(r_.i16, v);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i16) / sizeof(r_.i16[0])) ; i++) {
        r_.i16[i] = HEDLEY_STATIC_CAST(int16_t, a_.u8[i + 8]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_widen_high_u8x16(a) simde_wasm_i16x8_widen_high_u8x16((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i32x4_widen_high_u16x8 (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i32x4_widen_high_u16x8(a);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_CONVERT_VECTOR_)
      const uint16_t v SIMDE_VECTOR(8) = { a_.u16[4], a_.u16[5], a_.u16[6], a_.u16[7] };

      SIMDE_CONVERT_VECTOR_(r_.i32, v);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32) / sizeof(r_.i32[0])) ; i++) {
        r_.i32[i] = HEDLEY_STATIC_CAST(int32_t, a_.u16[i + 4]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_widen_high_u16x8(a) simde_wasm_i32x4_widen_high_u16x8((a))
#endif

/* X_load_Y */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i16x8_load_8x8 (const void * mem) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i16x8_load_8x8(mem);
  #else
    simde_v128_private r_;

    #if defined(SIMDE_CONVERT_VECTOR_)
      int8_t v SIMDE_VECTOR(8);
      simde_memcpy(&v, mem, sizeof(v));
      SIMDE_CONVERT_VECTOR_(r_.i16, v);
    #else
      SIMDE_ALIGN_TO_16 int8_t v[8];
      simde_memcpy(v, mem, sizeof(v));

      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i16) / sizeof(r_.i16[0])) ; i++) {
        r_.i16[i] = HEDLEY_STATIC_CAST(int16_t, v[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i16x8_load_8x8(mem) simde_wasm_i16x8_load_8x8((mem))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i32x4_load_16x4 (const void * mem) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i32x4_load_16x4(mem);
  #else
    simde_v128_private r_;

    #if defined(SIMDE_CONVERT_VECTOR_)
      int16_t v SIMDE_VECTOR(8);
      simde_memcpy(&v, mem, sizeof(v));
      SIMDE_CONVERT_VECTOR_(r_.i32, v);
    #else
      SIMDE_ALIGN_TO_16 int16_t v[4];
      simde_memcpy(v, mem, sizeof(v));

      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i32) / sizeof(r_.i32[0])) ; i++) {
        r_.i32[i] = HEDLEY_STATIC_CAST(int32_t, v[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_load_16x4(mem) simde_wasm_i32x4_load_16x4((mem))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i64x2_load_32x2 (const void * mem) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i64x2_load_32x2(mem);
  #else
    simde_v128_private r_;

    #if defined(SIMDE_CONVERT_VECTOR_)
      int32_t v SIMDE_VECTOR(8);
      simde_memcpy(&v, mem, sizeof(v));
      SIMDE_CONVERT_VECTOR_(r_.i64, v);
    #else
      SIMDE_ALIGN_TO_16 int32_t v[2];
      simde_memcpy(v, mem, sizeof(v));

      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.i64) / sizeof(r_.i64[0])) ; i++) {
        r_.i64[i] = HEDLEY_STATIC_CAST(int64_t, v[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i64x2_load_32x2(mem) simde_wasm_i64x2_load_32x2((mem))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u16x8_load_8x8 (const void * mem) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u16x8_load_8x8(mem);
  #else
    simde_v128_private r_;

    #if defined(SIMDE_CONVERT_VECTOR_)
      uint8_t v SIMDE_VECTOR(8);
      simde_memcpy(&v, mem, sizeof(v));
      SIMDE_CONVERT_VECTOR_(r_.u16, v);
    #else
      SIMDE_ALIGN_TO_16 uint8_t v[8];
      simde_memcpy(v, mem, sizeof(v));

      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u16) / sizeof(r_.u16[0])) ; i++) {
        r_.u16[i] = HEDLEY_STATIC_CAST(uint16_t, v[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u16x8_load_8x8(mem) simde_wasm_u16x8_load_8x8((mem))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u32x4_load_16x4 (const void * mem) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u32x4_load_16x4(mem);
  #else
    simde_v128_private r_;

    #if defined(SIMDE_CONVERT_VECTOR_)
      uint16_t v SIMDE_VECTOR(8);
      simde_memcpy(&v, mem, sizeof(v));
      SIMDE_CONVERT_VECTOR_(r_.u32, v);
    #else
      SIMDE_ALIGN_TO_16 uint16_t v[4];
      simde_memcpy(v, mem, sizeof(v));

      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u32) / sizeof(r_.u32[0])) ; i++) {
        r_.u32[i] = HEDLEY_STATIC_CAST(uint32_t, v[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u32x4_load_16x4(mem) simde_wasm_u32x4_load_16x4((mem))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u64x2_load_32x2 (const void * mem) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u64x2_load_32x2(mem);
  #else
    simde_v128_private r_;

    #if defined(SIMDE_CONVERT_VECTOR_)
      uint32_t v SIMDE_VECTOR(8);
      simde_memcpy(&v, mem, sizeof(v));
      SIMDE_CONVERT_VECTOR_(r_.u64, v);
    #else
      SIMDE_ALIGN_TO_16 uint32_t v[2];
      simde_memcpy(v, mem, sizeof(v));

      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.u64) / sizeof(r_.u64[0])) ; i++) {
        r_.u64[i] = HEDLEY_STATIC_CAST(uint64_t, v[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u64x2_load_32x2(mem) simde_wasm_u64x2_load_32x2((mem))
#endif

/* convert */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f32x4_convert_i32x4 (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f32x4_convert_i32x4(a);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_CONVERT_VECTOR_)
      SIMDE_CONVERT_VECTOR_(r_.f32, a_.i32);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f32) / sizeof(r_.f32[0])) ; i++) {
        r_.f32[i] = HEDLEY_STATIC_CAST(simde_float32, a_.i32[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f32x4_convert_i32x4(a) simde_wasm_f32x4_convert_i32x4((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_f32x4_convert_u32x4 (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_f32x4_convert_u32x4(a);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    #if defined(SIMDE_CONVERT_VECTOR_)
      SIMDE_CONVERT_VECTOR_(r_.f32, a_.u32);
    #else
      SIMDE_VECTORIZE
      for (size_t i = 0 ; i < (sizeof(r_.f32) / sizeof(r_.f32[0])) ; i++) {
        r_.f32[i] = HEDLEY_STATIC_CAST(simde_float32, a_.u32[i]);
      }
    #endif

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_f32x4_convert_u32x4(a) simde_wasm_f32x4_convert_u32x4((a))
#endif

/* trunc_saturate */

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_i32x4_trunc_saturate_f32x4 (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_i32x4_trunc_saturate_f32x4(a);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    SIMDE_VECTORIZE
    for (size_t i = 0 ; i < (sizeof(r_.i32) / sizeof(r_.i32[0])) ; i++) {
      if (simde_math_isnanf(a_.f32[i])) {
        r_.i32[i] = INT32_C(0);
      } else if (a_.f32[i] < HEDLEY_STATIC_CAST(simde_float32, INT32_MIN)) {
        r_.i32[i] = INT32_MIN;
      } else if (a_.f32[i] > HEDLEY_STATIC_CAST(simde_float32, INT32_MAX)) {
        r_.i32[i] = INT32_MAX;
      } else {
        r_.i32[i] = HEDLEY_STATIC_CAST(int32_t, a_.f32[i]);
      }
    }

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_i32x4_trunc_saturate_f32x4(a) simde_wasm_i32x4_trunc_saturate_f32x4((a))
#endif

SIMDE_FUNCTION_ATTRIBUTES
simde_v128_t
simde_wasm_u32x4_trunc_saturate_f32x4 (simde_v128_t a) {
  #if defined(SIMDE_WASM_SIMD128_NATIVE)
    return wasm_u32x4_trunc_saturate_f32x4(a);
  #else
    simde_v128_private
      a_ = simde_v128_to_private(a),
      r_;

    SIMDE_VECTORIZE
    for (size_t i = 0 ; i < (sizeof(r_.u32) / sizeof(r_.u32[0])) ; i++) {
      if (simde_math_isnanf(a_.f32[i]) ||
          a_.f32[i] < SIMDE_FLOAT32_C(0.0)) {
        r_.u32[i] = UINT32_C(0);
      } else if (a_.f32[i] > HEDLEY_STATIC_CAST(simde_float32, UINT32_MAX)) {
        r_.u32[i] = UINT32_MAX;
      } else {
        r_.u32[i] = HEDLEY_STATIC_CAST(uint32_t, a_.f32[i]);
      }
    }

    return simde_v128_from_private(r_);
  #endif
}
#if defined(SIMDE_WASM_SIMD128_ENABLE_NATIVE_ALIASES)
  #define wasm_u32x4_trunc_saturate_f32x4(a) simde_wasm_u32x4_trunc_saturate_f32x4((a))
#endif

SIMDE_END_DECLS_

HEDLEY_DIAGNOSTIC_POP

#endif /* !defined(SIMDE_WASM_SIMD128_H) */
