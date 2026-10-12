#ifndef _ERGM_DEBUG_H_
#define _ERGM_DEBUG_H_

#ifdef NDEBUG
#if defined(_MSC_VER)
#define ALWAYS_INLINE __forceinline
#elif defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER) ||    \
  defined(__GNUC__) || defined(__clang__)
#define ALWAYS_INLINE inline __attribute__((always_inline))
#else
#define ALWAYS_INLINE inline
#endif // Compiler version
#else
#define ALWAYS_INLINE
#endif // NDEBUG
#endif // _ERGM_DEBUG_H_
