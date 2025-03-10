#define _ETYPE1(FUN) FUN
#define _ETYPE2(PREFIX, FUN) PREFIX ## FUN
#define ETYPE(...) _GET_OVERRIDE2(__VA_ARGS__, _ETYPE2, _ETYPE1, )(__VA_ARGS__)
#define IFEWT(...)
#define IFELSEEWT(yes, no) no
#define EWTTYPE Rboolean
#define EWTRTYPE INTEGER
