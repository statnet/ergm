/*  File inst/include/ergm_kvec.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2024 Statnet Commons
 */
/* The MIT License

   Copyright (c) 2008, by Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/*
  An example:

#include "kvec.h"
int main() {
	kvec_t(int) array;
	kv_init(array);
	kv_push(int, array, 10); // append
	kv_a(int, array, 20) = 5; // dynamic --- avoid using to lookup values
	kv_A(array, 20) = 4; // static --- also use when reading
	kv_destroy(array);
	return 0;
}
*/

/*
  2024-09-27 (0.1.0.R.1):

    * kv_a() now works as an lvalue in C, at the cost of potentially causing problems in other situations. (Suggested by Jason Wang (https://github.com/wang-borong).)
    * New macros:
        * kv_blank can be used to initialise a kvec by assigning.
        * kv_del_plug(v, i) deletes the i'th element, moving the last element in the array to its position. O(1), but it doesn't preserve the order of the elements.
        * kv_del_shift(type, v, i) deletes the i'th element, shifting the remaining elements to fill the vacant space. O(n), but preserves the order of the elements.
        * kv_ins_shift(type, v, i, x) inserts an element at position i, shifting the elements to make room. O(n).

  2024-08-29 (0.1.0.R):

    * Adaptation to R memory management.
        * Using R_Free and R_Realloc instead of free and realloc.

  2008-09-22 (0.1.0):

    * The initial version.

*/

#ifndef _ERGM_KVEC_H_
#define _ERGM_KVEC_H_

#include <R.h>

#define kv_roundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

#define kvec_t(type) struct { size_t n, m; type *a; }
#define kv_init(v) ((v).n = (v).m = 0, (v).a = NULL) // Not necessary if the data structure is calloc()ed.
#define kv_blank {.n = 0, .m = 0, .a = NULL}
#define kv_destroy(v) {R_Free((v).a); (v).n = (v).m = 0;}
#define kv_A(v, i) ((v).a[(i)])
#define kv_pop(v) ((v).a[--(v).n])
#define kv_size(v) ((v).n)
#define kv_max(v) ((v).m)

#define kv_resize(type, v, s)  ((v).m = (s), (v).a = R_Realloc((v).a, (v).m, type))

#define kv_copy(type, v1, v0) do {							\
		if ((v1).m < (v0).n) kv_resize(type, v1, (v0).n);	\
		(v1).n = (v0).n;									\
		memcpy((v1).a, (v0).a, sizeof(type) * (v0).n);		\
	} while (0)												\

#define kv_push(type, v, x) do {									\
		if ((v).n == (v).m) {										\
			(v).m = (v).m? (v).m<<1 : 2;							\
			(v).a = R_Realloc((v).a, (v).m, type);	\
		}															\
		(v).a[(v).n++] = (x);										\
	} while (0)

#define kv_pushp(type, v) (((v).n == (v).m)?							\
						   ((v).m = ((v).m? (v).m<<1 : 2),				\
							(v).a = R_Realloc((v).a, (v).m, type), 0)	\
						   : 0), ((v).a + ((v).n++))

/* NB: Use primarily as lvalue. */
#define kv_a(type, v, i) ((v).m <= (size_t)(i)? \
						  ((v).m = (v).n = (i) + 1, kv_roundup32((v).m), \
						   (v).a = R_Realloc((v).a, (v).m, type), 0) \
						  : (v).n <= (size_t)(i)? (v).n = (i) + 1 \
						  : 0), (v).a[(i)]

#define kv_del_plug(v, i) if((i) != --(v).n) kv_A((v), (i)) = kv_A((v), (v).n)

#define kv_del_shift(type, v, i) memmove((v).a + (i), (v).a + (i) + 1, sizeof(type) * (--(v).n - (i)))

#define kv_ins_shift(type, v, i, x)  do {                               \
    if ((v).n == (v).m) {                                               \
      (v).m = (v).m? (v).m<<1 : 2;                                      \
      (v).a = R_Realloc((v).a, (v).m, type);                            \
    }                                                                   \
    memmove((v).a + (i) + 1, (v).a + (i), sizeof(type) * ((v).n++ - (i))); \
    (v).a[i] = (x);                                                     \
  } while (0)

#endif
