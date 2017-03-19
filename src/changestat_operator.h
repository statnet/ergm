#include "changestat.h"
#include "model.h"

typedef struct{void **aux_storage; Model *m;} StoreAuxAndModel;

unsigned char *unpack_strasdouble(double **x);
Model *unpack_Modelasdouble(double **x);
