#ifndef _ERGM_CLI_H_
#define _ERGM_CLI_H_
#include <cli/progress.h>

#define VERBOSE_ENCODE_CLI_BAR(verbose) -((verbose) + 1)


#define VERBOSE_CLI_BAR(verbose) ((verbose) < 0)


#define VERBOSE_DECODE_CLI_BAR(verbose) verbose = -verbose - 1;


#define CLI_BAR(name, max, text)                          \
  SEXP name = PROTECT(cli_progress_bar(max, NULL));       \
  cli_progress_set_name(name, text);


#define CLI_BAR_IF_VERBOSE(name, max, text)             \
  SEXP name = PROTECT(VERBOSE_CLI_BAR(verbose) ?        \
                      cli_progress_bar(max, NULL) :     \
                      R_NilValue);                      \
  if (name != R_NilValue) {                             \
    VERBOSE_DECODE_CLI_BAR(verbose);                    \
    cli_progress_set_name(name, text);                  \
  }


#define CLI_BAR_SET(name, value)                                        \
  if (name != R_NilValue && CLI_SHOULD_TICK) cli_progress_set(name, value);


#define CLI_BAR_FINISH(name)                            \
  if (name != R_NilValue) cli_progress_done(name);      \
  UNPROTECT(1);

#endif // _ERGM_CLI_H_
