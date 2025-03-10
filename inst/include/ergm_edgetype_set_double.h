#define _ETYPE1(FUN) Wt ## FUN
#define _ETYPE2(PREFIX, FUN) PREFIX ## Wt ## FUN
#define ETYPE(...) _GET_OVERRIDE2(__VA_ARGS__, _ETYPE2, _ETYPE1, )(__VA_ARGS__)
#define IFEWT(...) __VA_ARGS__
#define IFELSEEWT(yes, no) yes
#define EWTTYPE double
#define EWTRTYPE REAL
