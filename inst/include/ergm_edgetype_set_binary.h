/*  File inst/include/ergm_edgetype_set_binary.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#define _ETYPE1(FUN) FUN
#define _ETYPE2(PREFIX, FUN) PREFIX ## FUN
#define ETYPE(...) _GET_OVERRIDE2(__VA_ARGS__, _ETYPE2, _ETYPE1, )(__VA_ARGS__)
#define IFEWT(...)
#define IFELSEEWT(yes, no) no
#define EWTTYPE Rboolean
#define EWTRTYPE INTEGER
