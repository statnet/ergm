#include "ergm_Rutil.h"
#include "ergm_constants.h"
#include "ergm_kvec.h"
#include "ergm_khash.h"

unsigned int GetBuiltErgmAPIMajor(void){
  return ERGM_API_MAJOR;
}

unsigned int GetBuiltErgmAPIMinor(void){
  return ERGM_API_MINOR;
}

static kvec_t(khint_t) build_version_warned = kv_blank;

void warn_API_version(const char *sn){

  /* Here, we are saving the hash of the library name rather than the
     name itself, since keeping track of strings allocated from the
     heap is a hassle, and we don't want to assume that sn will always
     be around. Chances of collision are negligible. */
  khint_t snhash = kh_str_hash_func(sn);
  unsigned int i = 0;
  while(i < kv_size(build_version_warned) && kv_A(build_version_warned, i) != snhash) i++;
  if(i < kv_size(build_version_warned)) return;

  kv_push(khint_t, build_version_warned, snhash);

  unsigned int (*bv)() = (unsigned int (*)()) R_FindSymbol("GetBuiltErgmAPIMajor",sn,NULL);
  unsigned int bv_major = bv ? bv() : 0;
  if(bv_major){
    unsigned int bv_minor = ((unsigned int (*)())R_FindSymbol("GetBuiltErgmAPIMinor",sn,NULL))();
    if(bv_major != ERGM_API_MAJOR || bv_minor != ERGM_API_MINOR)
      warningcall_immediate(R_NilValue, "Package '%s' was compiled against 'ergm' with C API version\n    %d.%d, but it is being loaded by 'ergm' with C API version %d.%d. Inconsistent\n    versions may result in malfunctions ranging from incorrect results to R\n    crashing. Please rebuild the package against the current 'ergm' version.",
              sn, bv_major, bv_minor, ERGM_API_MAJOR, ERGM_API_MINOR);
  }else
    warningcall_immediate(R_NilValue, "Package '%s' was compiled against an older version of 'ergm'\n    that did not store '%s''s C API version information. Inconsistent\n    versions may cause malfunctions ranging from incorrect results to R\n    crashing. Please rebuild the package against the current 'ergm' version.",
                          sn, sn);
}
