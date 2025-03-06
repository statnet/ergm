#define EDGETYPE(FUN) Wt ## FUN
#define IFEDGEWT(...) __VA_ARGS__
#define IFELSEEDGEWT(yes, no) yes
#define EDGEWTTYPE double
#define EDGEWTRTYPE REAL
