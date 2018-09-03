#ifndef _KHASH_HELPERS_H_
#define _KHASH_HELPERS_H_

#include "khash.h"

#define kh_getval(name, h, k, def, into)				\
  {									\
    khiter_t _kh_getval_pos = kh_get(name, h, k);			\
    into = _kh_getval_pos==kh_end(h) ? def : kh_value(h, _kh_getval_pos); \
  }

#define kh_set(name, h, k, val)					\
  {								\
    /* Rprintf("(%d,%d)<-%u\n",k.tail,k.head,val);			 */\
    khiter_t _kh_set_pos = kh_get(name, h, k);			\
    int _kh_set_ret;						\
    if(_kh_set_pos==kh_end(h)){					\
      _kh_set_pos = kh_put(name, h, k, &_kh_set_ret);		\
    }								\
    kh_val(h, _kh_set_pos) = val;				\
  }

#define kh_unset(name, h, k)					\
  {								\
    /* Rprintf("(%d,%d)<-X\n",k.tail,k.head);			 */\
    khiter_t _kh_unset_pos = kh_get(name, h, k);		   \
    if(_kh_unset_pos!=kh_end(h)){				\
      kh_del(name, h, _kh_unset_pos);				\
    }								\
  }

#endif
